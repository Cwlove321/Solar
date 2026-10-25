"""Strict preprocessed-input checks for portable, isolated cases."""
from __future__ import annotations

import json
import math
import struct
from pathlib import Path

import numpy as np
from astropy.io import fits


def number(value, name: str, *, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be a finite number")
    if positive and value <= 0:
        raise ValueError(f"{name} must be positive")
    return float(value)


def integer(value, name: str, *, minimum: int = 1) -> int:
    result = number(value, name)
    if not result.is_integer() or result < minimum:
        raise ValueError(f"{name} must be an integer >= {minimum}")
    return int(result)


def audit_case_inputs(root: Path) -> dict:
    """Recheck actual arrays in an existing case; raise on any inconsistency.

    Returns the same derived-input dictionary as inspect_inputs. This is a
    deep geometry/data check, independent of file-manifest integrity checks.
    """
    root = Path(root).resolve()
    resolved = json.loads((root/"case.resolved.json").read_text(encoding="utf-8"))
    locations = []
    for key,default in (("boundary_dir","inputs/boundary"),("path_dat","inputs/axis/path.dat")):
        path = (root/resolved.get(key,default)).resolve()
        if not path.is_relative_to(root):
            raise ValueError(f"case {key} must refer to an input snapshot inside the case directory")
        locations.append(path)
    checked = inspect_inputs(*locations,resolved["grid"],resolved["rbsl"],allow_reference_time_mismatch=resolved.get("allow_reference_time_mismatch",False))
    for key,value in checked["global_pfss"].items():
        if resolved["global_pfss"].get(key) != value:
            raise ValueError(f"case global_pfss.{key} disagrees with actual input geometry")
    return checked


def inspect_inputs(boundary: Path, axis: Path, grid_settings: dict, rbsl: dict, *, allow_reference_time_mismatch: bool = False) -> dict:
    """Check actual binary/FITS arrays and derive the common local/global grid.

    Boundary coordinates follow prepare_rbsl_boundary.py: lon is Stonyhurst,
    global phi=lon+180 degrees, theta=90-lat degrees. path.dat contains those
    phi degrees, latitude degrees, and radius in Rsun as little-endian triples.
    """
    metadata = json.loads((boundary / "parameters.json").read_text(encoding="utf-8"))
    bg = metadata["grid"]
    nlat = integer(bg["nlat"], "boundary.grid.nlat")
    nlon = integer(bg["nlon"], "boundary.grid.nlon")
    resolution = number(bg["resolution_degree"], "boundary.grid.resolution_degree", positive=True)
    lon0, lon1 = [number(v, "boundary longitude") for v in bg["lon_range"]]
    lat0, lat1 = [number(v, "boundary latitude") for v in bg["lat_range"]]
    if not (-180 <= lon0 < lon1 <= 180 and -90 < lat0 < lat1 < 90):
        raise ValueError("boundary must be a single non-wrapping longitude wedge in [-180,180] and avoid poles")
    if not np.allclose([(lon1-lon0)/resolution, (lat1-lat0)/resolution], [nlon, nlat], atol=1e-7, rtol=0):
        raise ValueError("boundary angular ranges/resolution disagree with nlon/nlat")
    x0 = integer(bg["x0"], "boundary.grid.x0", minimum=0)
    y0 = integer(bg["y0"], "boundary.grid.y0", minimum=0)
    if not np.allclose([(lon0+180)/resolution, (lat0+90)/resolution], [x0,y0], atol=1e-7, rtol=0):
        raise ValueError("boundary patch indices disagree with longitude/latitude ranges; re-run preprocessing")

    global_path = boundary / "pfss_boundary.dat"
    with global_path.open("rb") as stream:
        header = stream.read(8)
        if len(header) != 8:
            raise ValueError("pfss_boundary.dat has an incomplete header")
        nphi, ntheta = struct.unpack("<ii", header)
        if nphi < 1 or ntheta < 1 or nphi*ntheta > 100_000_000:
            raise ValueError("invalid global boundary shape")
        expected_bytes = 8 + 8*(nphi+ntheta+nphi*ntheta)
        if global_path.stat().st_size != expected_bytes:
            raise ValueError("pfss_boundary.dat byte count disagrees with its header")
        theta = np.fromfile(stream, "<f8", count=ntheta)
        phi = np.fromfile(stream, "<f8", count=nphi)
        global_br = np.fromfile(stream, "<f8").reshape(ntheta,nphi)
    if not all(np.all(np.isfinite(a)) for a in (theta,phi,global_br)):
        raise ValueError("global boundary contains NaN/Inf")
    if not np.allclose([nphi*resolution,ntheta*resolution], [360,180], atol=1e-7, rtol=0):
        raise ValueError("global boundary dimensions disagree with resolution_degree")
    expected_theta = np.radians(180-(np.arange(ntheta)+.5)*resolution)
    expected_phi = np.radians((np.arange(nphi)+.5)*resolution)
    if not np.allclose(theta,expected_theta,atol=1e-10,rtol=0) or not np.allclose(phi,expected_phi,atol=1e-10,rtol=0):
        raise ValueError("global boundary coordinate orientation/centres do not match the preprocessor convention")
    if x0+nlon > nphi or y0+nlat > ntheta:
        raise ValueError("local patch is outside global boundary")

    brll_path = boundary / "Brll0000.fits"
    with fits.open(brll_path, memmap=False) as hdus:
        br = np.asarray(hdus[0].data)
        fh = hdus[0].header
        if str(fh.get("BUNIT","")).strip().lower() not in ("g","gauss"):
            raise ValueError("Brll0000.fits BUNIT must explicitly be G/gauss")
        if br.shape != (nlon,nlat) or int(fh.get("NGHOST",-1)) != 0:
            raise ValueError(f"Brll0000.fits must have no ghosts and Python shape {(nlon,nlat)}")
        for key,expected in (("NTHETA",nlat),("NPHI",nlon)):
            if key in fh and int(fh[key]) != expected:
                raise ValueError(f"FITS {key} disagrees with actual physical dimensions")
        if not np.all(np.isfinite(br)):
            raise ValueError("Brll0000.fits contains NaN/Inf")
        for key, expected in zip(("LONMIN","LONMAX","LATMIN","LATMAX"),(lon0,lon1,lat0,lat1)):
            if key not in fh or not math.isclose(float(fh[key]),expected,abs_tol=1e-8,rel_tol=0):
                raise ValueError(f"FITS {key} disagrees with parameters.json")
        if not np.allclose(global_br[y0:y0+nlat,x0:x0+nlon],br.T[::-1,:],atol=1e-6,rtol=1e-6):
            raise ValueError("Br FITS orientation/values disagree with the global PFSS patch")

    level = integer(grid_settings["refine_max_level"], "grid.refine_max_level")
    if level > 20:
        raise ValueError("refine_max_level must not exceed the 20-entry method arrays")
    divisor = 2**(level-1)
    derived = dict(grid_settings)
    for dimension, count in ((2,nlat),(3,nlon)):
        base, remainder = divmod(count,divisor)
        if remainder or base % 12:
            raise ValueError(f"boundary dimension {count} must equal domain_nx{dimension} * 2^(Lmax-1), with domain_nx{dimension} divisible by block size 12")
        if f"domain_nx{dimension}" in grid_settings and grid_settings[f"domain_nx{dimension}"] != base:
            raise ValueError(f"grid.domain_nx{dimension} conflicts with the preprocessed boundary; remove it to derive {base}")
        derived[f"domain_nx{dimension}"] = base
    derived["domain_nx1"] = integer(grid_settings["domain_nx1"], "grid.domain_nx1")
    if derived["domain_nx1"] % 12:
        raise ValueError("grid.domain_nx1 must be divisible by block size 12")
    rmin = number(grid_settings["xprobmin1"], "grid.xprobmin1",positive=True)
    rmax = number(grid_settings["xprobmax1"], "grid.xprobmax1",positive=True)
    if rmin >= rmax:
        raise ValueError("grid.xprobmin1 must be below xprobmax1")
    if not math.isclose(rmin,1.,rel_tol=0,abs_tol=1e-12):
        raise ValueError("grid.xprobmin1 must be 1: this spherical atmosphere/PFSS implementation is anchored to the solar surface")
    if rmax > 2.+1e-12 or not math.isclose(1/(rmax-rmin),round(1/(rmax-rmin)),rel_tol=0,abs_tol=1e-10):
        raise ValueError("the local PDFI_SS radial sampling requires xprobmax1<=2 and 1/(xprobmax1-1) to be an integer")
    q = number(grid_settings["qstretch_baselevel1"], "grid.qstretch_baselevel1",positive=True)
    if math.isclose(q,1.,rel_tol=0,abs_tol=1e-12):
        raise ValueError("qstretch_baselevel1=1 is singular with stretch_dim='uni'; choose a non-unit geometric stretching ratio")
    angular = {"xprobmin2":(90-lat1)/360,"xprobmax2":(90-lat0)/360,"xprobmin3":(lon0+180)/360,"xprobmax3":(lon1+180)/360}
    for key,value in angular.items():
        if key in grid_settings and not math.isclose(number(grid_settings[key],f"grid.{key}"),value,abs_tol=1e-10,rel_tol=0):
            raise ValueError(f"grid.{key} conflicts with boundary metadata; remove this old event setting")
    derived.update(angular)

    size = axis.stat().st_size
    if size % 24 or size < 48:
        raise ValueError("path.dat must contain at least two interleaved little-endian float64 triples (24 bytes each)")
    npoints = size//24
    summary_path = axis.parent / "path_summary.json"
    if summary_path.exists():
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        if integer(summary["npoints"],"path_summary.npoints",minimum=2) != npoints:
            raise ValueError("path_summary.npoints disagrees with actual path.dat bytes")
        if metadata.get("reference_time") and summary.get("reference_time"):
            from astropy.time import Time
            try:
                offset = abs((Time(metadata["reference_time"])-Time(summary["reference_time"])).to_value("s"))
            except (ValueError,TypeError) as exc:
                raise ValueError("boundary/axis reference_time is not a valid observation time") from exc
            if offset > 1e-3 and not allow_reference_time_mismatch:
                raise ValueError("boundary and axis reference_time differ; align them or explicitly set allow_reference_time_mismatch=true")
    if "npoints" in rbsl and integer(rbsl["npoints"],"rbsl.npoints",minimum=2) != npoints:
        raise ValueError("rbsl.npoints disagrees with actual path.dat bytes; remove it to derive the count")
    points = np.fromfile(axis,"<f8").reshape(-1,3)
    if not np.all(np.isfinite(points)) or np.any(points[:,2] <= 0):
        raise ValueError("path.dat contains NaN/Inf or nonpositive radius")
    if np.any(points[:,0] < lon0+180-1e-8) or np.any(points[:,0] > lon1+180+1e-8) or np.any(points[:,1] < lat0-1e-8) or np.any(points[:,1] > lat1+1e-8):
        raise ValueError("path.dat lies outside the boundary wedge; expected columns (longitude+180 degrees, latitude degrees, radius/Rsun)")
    if np.any(points[:,2] < rmin-1e-8) or np.any(points[:,2] > rmax+1e-8):
        raise ValueError("path.dat lies outside the configured radial domain (photospheric footpoints r=1 require xprobmin1<=1)")
    if np.any(np.linalg.norm(np.diff(points,axis=0),axis=1) == 0):
        raise ValueError("path.dat contains repeated consecutive axis points")
    return {"grid":derived,"npoints":npoints,"boundary_metadata":metadata,"global_pfss":{"patch_x0":x0,"patch_y0":y0,"expected_nphi":nphi,"expected_ntheta":ntheta}}

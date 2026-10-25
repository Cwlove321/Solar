#!/usr/bin/env python3
"""Run solar_code AIA preparation, project to SHARP CEA, and write a VTK plane."""

from __future__ import annotations

import argparse
import json
import warnings
from pathlib import Path

import astropy.io.ascii as astropy_ascii
import astropy.units as u
import numpy as np
import sunpy.map
from astropy.table import QTable
from astropy.time import Time
from erfa import ErfaWarning
from sunpy.time import parse_time

from solarutils.obs.maps import aia_prep
from aia_npy_to_bottom_boundary_vtk import parse_amrvac_par, write_plane
from project_aia_to_sharp import reproject_to_sharp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aia", type=Path, required=True)
    parser.add_argument("--sharp", type=Path, required=True)
    parser.add_argument("--par", type=Path, required=True)
    parser.add_argument("--pointing-table", type=Path, required=True)
    parser.add_argument("--psf", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--x0", type=int, default=230)
    parser.add_argument("--x1", type=int, default=950)
    parser.add_argument("--y0", type=int, default=80)
    parser.add_argument("--y1", type=int, default=400)
    parser.add_argument("--deconvolve-iterations", type=int, default=25)
    parser.add_argument("--rotation-model", default="howard")
    return parser.parse_args()


def load_pointing_table(path: Path) -> QTable:
    """Load the LMSAL AIA pointing CSV with the units expected by aiapy."""

    table = QTable(astropy_ascii.read(path))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=ErfaWarning)
        table["T_START"] = parse_time(table["T_START"])
        table["T_STOP"] = parse_time(table["T_STOP"])
    for column in table.colnames:
        if "X0" in column or "Y0" in column:
            table[column].unit = "pixel"
        if "IMSCALE" in column:
            table[column].unit = "arcsecond / pixel"
        if "INSTROT" in column:
            table[column].unit = "degree"
        if any(name in column for name in ("X0", "Y0", "IMSCALE", "INSTROT")) and hasattr(table[column], "mask"):
            table[column] = table[column].filled(np.nan)
    return table


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    aia = sunpy.map.Map(args.aia)
    pointing = load_pointing_table(args.pointing_table)
    prepared = aia_prep(
        aia,
        pointing_table=pointing,
        deconvolve_psf=True,
        psf=args.psf,
        deconvolve_iterations=args.deconvolve_iterations,
        register_map=True,
    )

    timestamp = Time(aia.date).strftime("%Y%m%dT%H%M%S")
    wavelength = int(round(aia.wavelength.to_value(u.angstrom)))
    tag = f"aia_{timestamp}_{wavelength:03d}_prep_psf{args.deconvolve_iterations}"
    prepared_fits = args.out / f"{tag}.fits"
    prepared.save(prepared_fits, overwrite=True)

    sharp = sunpy.map.Map(args.sharp)
    target = sharp.submap(
        (args.x0, args.y0) * u.pix,
        top_right=(args.x1 - 1, args.y1 - 1) * u.pix,
    )
    projected = reproject_to_sharp(
        prepared,
        target,
        differential_rotation=True,
        rotation_model=args.rotation_model,
    )
    data = np.asarray(projected.data, dtype=np.float32)
    expected_shape = (args.y1 - args.y0, args.x1 - args.x0)
    if data.shape != expected_shape:
        raise RuntimeError(f"projected shape {data.shape} != requested {expected_shape}")
    if not np.isfinite(data).all():
        raise RuntimeError("projected array contains non-finite values")

    projected_npy = args.out / f"{tag}_sharp_cea_720x320.npy"
    np.save(projected_npy, data)

    params = parse_amrvac_par(args.par)
    grid_scale_x = data.shape[1] // int(params["domain_nx1"])
    grid_scale_y = data.shape[0] // int(params["domain_nx2"])
    if grid_scale_x != grid_scale_y or grid_scale_x < 1:
        raise RuntimeError("projected array is not an integer isotropic multiple of the par grid")

    log_data = np.log10(np.maximum(data, 0.0) + 1.0).astype(np.float32)
    vtk_path = args.out / f"{tag}_bottom_x3min_720x320.vtk"
    vtk_info = write_plane(
        log_data,
        vtk_path,
        params,
        wavelength,
        aia.date.isot,
        aia.date.isot,
        args.aia.name,
        "log10p1",
        grid_scale_x,
    )

    metadata = {
        "source_aia": str(args.aia),
        "source_sharp": str(args.sharp),
        "amrvac_par": str(args.par),
        "preprocessing": {
            "implementation": "solarutils.obs.maps.aia_prep",
            "steps": ["update_pointing", "psf_deconvolve", "register"],
            "pointing_table": str(args.pointing_table),
            "psf": str(args.psf),
            "deconvolve_iterations": args.deconvolve_iterations,
        },
        "projection": {
            "target": "SHARP CEA",
            "sharp_pixel_slice_python": {"x": [args.x0, args.x1], "y": [args.y0, args.y1]},
            "shape_yx": list(data.shape),
            "differential_rotation": True,
            "rotation_model": args.rotation_model,
            "sharp_reference_time": Time(target.coordinate_frame.obstime).utc.isot,
        },
        "transform": "log10(max(I, 0) + 1)",
        "prepared_fits": prepared_fits.name,
        "projected_npy": projected_npy.name,
        "vtk": vtk_info,
    }
    (args.out / f"{tag}_bottom_x3min_720x320.metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()

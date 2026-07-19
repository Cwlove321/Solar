#!/usr/bin/env python
"""Make a SHARP-like CEA vector magnetic patch from full-disk HMI vector data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.ndimage import map_coordinates
from sunpy.coordinates import frames, sun
from sunpy.map import Map

try:
    from .nlfff_core.cea import (
        gaussian_block_reduce,
        heliographic_to_local_cartesian_boundary,
        hmi_native_vector_components,
        make_cea_patch_grid,
        native_to_heliographic_components,
    )
    from .nlfff_core.hmi import hmi_disambiguate_azimuth
except ImportError:  # Direct execution: python Preprocess/make_cea_patch.py
    from nlfff_core.cea import (
        gaussian_block_reduce,
        heliographic_to_local_cartesian_boundary,
        hmi_native_vector_components,
        make_cea_patch_grid,
        native_to_heliographic_components,
    )
    from nlfff_core.hmi import hmi_disambiguate_azimuth


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("Loading HMI vector segments...")
    field_map = Map(args.field)
    inclination_map = Map(args.inclination)
    azimuth_map = Map(args.azimuth)
    disambig_map = Map(args.disambig)

    print("Applying HMI disambiguation and forming native image-plane components...")
    corrected_azimuth = hmi_disambiguate_azimuth(
        azimuth_map.data,
        disambig_map.data,
        method=args.disambig_method,
    )
    bxi, beta, bzeta = hmi_native_vector_components(field_map.data, inclination_map.data, corrected_azimuth)

    print("Building target CEA grid...")
    grid = make_cea_patch_grid(
        args.center_lon,
        args.center_lat,
        args.width_degree,
        args.height_degree,
        args.resolution_degree,
    )
    print("Sampling native components on the CEA grid...")
    sampled_bxi, sampled_beta, sampled_bzeta, sampling_metadata = sample_native_components(
        field_map,
        bxi,
        beta,
        bzeta,
        grid,
        args,
    )
    sampled_disambig = None if args.components_only else sample_nearest_segment(disambig_map, grid)

    disk_lat = args.disk_latitude
    if disk_lat is None:
        disk_lat = float(field_map.observer_coordinate.lat.to_value(u.deg))
    p_angle = args.p_angle
    if p_angle is None:
        p_angle = infer_p_angle_degree(field_map)

    print("Transforming sampled vectors to SHARP-like Br/Btheta/Bphi...")
    br, bt, bp = native_to_heliographic_components(
        sampled_bxi,
        sampled_beta,
        sampled_bzeta,
        grid.lon,
        grid.lat,
        disk_lon_degree=args.disk_lon,
        disk_lat_degree=disk_lat,
        p_angle_degree=p_angle,
    )
    bx, by, bz = heliographic_to_local_cartesian_boundary(br, bt, bp)

    print("Writing outputs...")
    write_component_fits(output_dir / "Br.fits", br, field_map, grid, "Br", disk_lat, p_angle, args)
    write_component_fits(output_dir / "Bt.fits", bt, field_map, grid, "Bt", disk_lat, p_angle, args)
    write_component_fits(output_dir / "Bp.fits", bp, field_map, grid, "Bp", disk_lat, p_angle, args)
    if not args.components_only:
        write_component_fits(output_dir / "disambig.fits", sampled_disambig, field_map, grid, "disambig", disk_lat, p_angle, args)
        np.savez_compressed(
            output_dir / "cea_patch.npz",
            x=grid.x,
            y=grid.y,
            lon=grid.lon,
            lat=grid.lat,
            bxi=sampled_bxi,
            beta=sampled_beta,
            bzeta=sampled_bzeta,
            br=br,
            bt=bt,
            bp=bp,
            bx=bx,
            by=by,
            bz=bz,
            disambig=sampled_disambig,
        )
        np.savez_compressed(output_dir / "cartesian_boundary.npz", bx=bx, by=by, bz=bz)
        write_cartesian_boundary_ascii(output_dir / "cartesian_boundary.dat", bx, by, bz)

    metadata = {
        "inputs": {
            "field": str(args.field),
            "inclination": str(args.inclination),
            "azimuth": str(args.azimuth),
            "disambig": str(args.disambig),
        },
        "grid": {
            "center_lon": args.center_lon,
            "center_lat": args.center_lat,
            "width_degree": args.width_degree,
            "height_degree": args.height_degree,
            "resolution_degree": args.resolution_degree,
            "shape": list(br.shape),
        },
        "geometry": {
            "disk_lon": args.disk_lon,
            "disk_lat": disk_lat,
            "p_angle": p_angle,
            "longitude_frame": "HeliographicStonyhurst",
        },
        "sampling": sampling_metadata
        | {
            "interpolation_order": args.interpolation_order,
            "fill_value": args.fill_value,
            "disambig_method": args.disambig_method,
            "disambig_remap": "nearest_neighbor",
        },
        "cartesian_boundary": {
            "Bx": "Bp, positive along increasing longitude / CEA x",
            "By": "-Bt, positive northward / CEA y",
            "Bz": "Br, positive radially outward",
        },
        "components_only": bool(args.components_only),
    }
    (output_dir / "parameters.json").write_text(json.dumps(json_safe(metadata), indent=2, allow_nan=False), encoding="utf-8")

    if args.quicklook:
        plot_quicklook(output_dir / "quicklook.png", br, bt, bp, bx, by, bz, vmax=args.vmax)

    print(f"Saved CEA patch products in {output_dir}")
    return 0


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field", required=True, type=Path, help="HMI vector field FITS segment.")
    parser.add_argument("--inclination", required=True, type=Path, help="HMI inclination FITS segment.")
    parser.add_argument("--azimuth", required=True, type=Path, help="HMI azimuth FITS segment.")
    parser.add_argument("--disambig", required=True, type=Path, help="HMI disambiguation FITS segment.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for CEA patch products.")
    parser.add_argument("--center-lon", required=True, type=float, help="Patch center Stonyhurst longitude in deg.")
    parser.add_argument("--center-lat", required=True, type=float, help="Patch center Stonyhurst latitude in deg.")
    parser.add_argument("--width-degree", required=True, type=float, help="CEA patch width in degrees.")
    parser.add_argument("--height-degree", required=True, type=float, help="CEA patch height in degrees.")
    parser.add_argument("--resolution-degree", type=float, default=0.03, help="CEA pixel size in degrees.")
    parser.add_argument("--disambig-method", type=int, default=2, help="HMI disambiguation bit method, usually 2.")
    parser.add_argument("--interpolation-order", type=int, default=3, choices=range(0, 6), help="Spline order.")
    parser.add_argument(
        "--sampling-mode",
        choices=("direct", "sharp"),
        default="direct",
        help="direct samples the target grid once; sharp oversamples and Gaussian-resizes like SHARP CEA.",
    )
    parser.add_argument(
        "--oversample-resolution-degree",
        type=float,
        default=0.01,
        help="Requested high-resolution CEA pixel size for --sampling-mode sharp.",
    )
    parser.add_argument(
        "--smooth-sigma-degree",
        type=float,
        default=0.01,
        help="Gaussian sigma in degrees for --sampling-mode sharp.",
    )
    parser.add_argument(
        "--smooth-truncate",
        type=float,
        default=2.0,
        help="Gaussian kernel truncation in sigma units for --sampling-mode sharp.",
    )
    parser.add_argument("--fill-value", type=float, default=np.nan, help="Value for off-disk/out-of-frame samples.")
    parser.add_argument("--disk-lon", type=float, default=0.0, help="Disk-center longitude in the chosen lon frame.")
    parser.add_argument("--disk-latitude", type=float, default=None, help="Override disk-center latitude in deg.")
    parser.add_argument("--p-angle", type=float, default=None, help="Override solar P angle in deg.")
    parser.add_argument("--quicklook", action="store_true", default=True, help="Write quicklook.png.")
    parser.add_argument("--no-quicklook", dest="quicklook", action="store_false")
    parser.add_argument("--vmax", type=float, default=1500.0, help="Quicklook symmetric color limit in G.")
    parser.add_argument(
        "--components-only",
        action="store_true",
        help="Write only Br/Bt/Bp FITS and parameters.json; omit large NPZ/ASCII intermediate products.",
    )
    return parser.parse_args(argv)


def sample_native_components(source_map, bxi, beta, bzeta, grid, args):
    if args.sampling_mode == "direct":
        pixel_yx = lonlat_to_source_pixels(source_map, grid.lon, grid.lat)
        return (
            sample_array(bxi, pixel_yx, grid.shape, order=args.interpolation_order, cval=args.fill_value),
            sample_array(beta, pixel_yx, grid.shape, order=args.interpolation_order, cval=args.fill_value),
            sample_array(bzeta, pixel_yx, grid.shape, order=args.interpolation_order, cval=args.fill_value),
            {"mode": "direct", "oversample_factor": 1},
        )

    factor = int(np.round(grid.resolution_degree / args.oversample_resolution_degree))
    if factor < 1:
        raise ValueError("--oversample-resolution-degree must be no larger than --resolution-degree")
    high_resolution = grid.resolution_degree / factor
    high_grid = make_cea_patch_grid(
        grid.center_lon,
        grid.center_lat,
        grid.shape[1] * grid.resolution_degree,
        grid.shape[0] * grid.resolution_degree,
        high_resolution,
    )
    pixel_yx = lonlat_to_source_pixels(source_map, high_grid.lon, high_grid.lat)
    high_bxi = sample_array(bxi, pixel_yx, high_grid.shape, order=args.interpolation_order, cval=args.fill_value)
    high_beta = sample_array(beta, pixel_yx, high_grid.shape, order=args.interpolation_order, cval=args.fill_value)
    high_bzeta = sample_array(bzeta, pixel_yx, high_grid.shape, order=args.interpolation_order, cval=args.fill_value)

    sigma_pixel = args.smooth_sigma_degree / high_resolution
    sampled_bxi = gaussian_block_reduce(
        high_bxi,
        factor,
        factor,
        sigma_y=sigma_pixel,
        sigma_x=sigma_pixel,
        truncate=args.smooth_truncate,
        cval=args.fill_value,
    )
    sampled_beta = gaussian_block_reduce(
        high_beta,
        factor,
        factor,
        sigma_y=sigma_pixel,
        sigma_x=sigma_pixel,
        truncate=args.smooth_truncate,
        cval=args.fill_value,
    )
    sampled_bzeta = gaussian_block_reduce(
        high_bzeta,
        factor,
        factor,
        sigma_y=sigma_pixel,
        sigma_x=sigma_pixel,
        truncate=args.smooth_truncate,
        cval=args.fill_value,
    )
    metadata = {
        "mode": "sharp",
        "oversample_factor": factor,
        "requested_oversample_resolution_degree": args.oversample_resolution_degree,
        "actual_oversample_resolution_degree": high_resolution,
        "smooth_sigma_degree": args.smooth_sigma_degree,
        "smooth_sigma_pixel": sigma_pixel,
        "smooth_truncate": args.smooth_truncate,
        "note": "Approximation of SHARP: sample native components on a fine CEA grid, Gaussian smooth, block-reduce.",
    }
    return sampled_bxi, sampled_beta, sampled_bzeta, metadata


def sample_nearest_segment(segment_map, grid):
    pixel_yx = lonlat_to_source_pixels(segment_map, grid.lon, grid.lat)
    return sample_array(segment_map.data, pixel_yx, grid.shape, order=0, cval=np.nan)


def lonlat_to_source_pixels(smap, lon, lat):
    coords = SkyCoord(
        np.asarray(lon).ravel() * u.deg,
        np.asarray(lat).ravel() * u.deg,
        frame=frames.HeliographicStonyhurst,
        obstime=smap.date,
        rsun=smap.rsun_meters,
        observer="earth",
    )
    hpc = coords.transform_to(smap.coordinate_frame)
    pixel_x, pixel_y = smap.world_to_pixel(hpc)
    shape = np.asarray(lon).shape
    return np.vstack([pixel_y.value.reshape(shape).ravel(), pixel_x.value.reshape(shape).ravel()])


def sample_array(data, pixel_yx, shape, *, order: int, cval: float):
    coords = np.asarray(pixel_yx, dtype=float).copy()
    invalid = ~np.isfinite(coords).all(axis=0)
    coords[:, invalid] = -1.0
    data = np.asarray(data, dtype=float)
    finite_data = np.isfinite(data)
    clean_data = np.nan_to_num(data, nan=0.0, posinf=0.0, neginf=0.0)
    sampled = map_coordinates(
        clean_data,
        coords,
        order=order,
        mode="constant",
        cval=0.0,
    )
    sampled_weight = map_coordinates(
        finite_data.astype(float),
        coords,
        order=0,
        mode="constant",
        cval=0.0,
    )
    sampled[invalid | (sampled_weight < 0.5)] = cval
    return sampled.reshape(shape)


def infer_p_angle_degree(smap):
    try:
        return float(-sun.P(smap.date).to_value(u.deg))
    except Exception:
        return float(-smap.meta.get("crota2", smap.meta.get("CROTA2", 0.0)))


def component_header(source_map, grid, component: str, disk_lat: float, p_angle: float, args):
    header = fits.Header()
    header["SIMPLE"] = True
    header["BITPIX"] = -64
    header["NAXIS"] = 2
    header["NAXIS1"] = grid.shape[1]
    header["NAXIS2"] = grid.shape[0]
    header["CTYPE1"] = "HGLN-CEA"
    header["CTYPE2"] = "HGLT-CEA"
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["CRPIX1"] = (grid.shape[1] + 1.0) / 2.0
    header["CRPIX2"] = (grid.shape[0] + 1.0) / 2.0
    header["CRVAL1"] = grid.center_lon
    header["CRVAL2"] = grid.center_lat
    header["CDELT1"] = grid.resolution_degree
    header["CDELT2"] = grid.resolution_degree
    header["CROTA2"] = 0.0
    header["BUNIT"] = "G"
    header["CONTENT"] = f"CEA {component}"
    header["COMP"] = component
    header["DATE-OBS"] = str(source_map.date.isot)
    header["T_OBS"] = source_map.date.strftime("%Y.%m.%d_%H:%M:%S")
    header["T_REC"] = str(source_map.meta.get("t_rec", header["T_OBS"]))
    header["HGLN_OBS"] = float(args.disk_lon)
    header["HGLT_OBS"] = float(disk_lat)
    header["P_ANGLE"] = float(p_angle)
    header["DSUN_OBS"] = float(source_map.dsun.to_value(u.m))
    header["RSUN_OBS"] = float(source_map.rsun_obs.to_value(u.arcsec))
    header["RSUN_REF"] = float(source_map.rsun_meters.to_value(u.m))
    header["HISTORY"] = (
        "Created by Cartesian_TMF/Preprocess/"
        "make_cea_patch.py"
    )
    header["HISTORY"] = "SHARP-like CEA remap from full-disk HMI vector segments."
    return header


def write_component_fits(path, data, source_map, grid, component: str, disk_lat: float, p_angle: float, args):
    header = component_header(source_map, grid, component, disk_lat, p_angle, args)
    fits.PrimaryHDU(np.asarray(data, dtype=np.float64), header=header).writeto(path, overwrite=True)


def json_safe(value):
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


def write_cartesian_boundary_ascii(path, bx, by, bz):
    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    bz = np.asarray(bz, dtype=float)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("bx, by, and bz must have the same shape")
    with Path(path).open("w", encoding="utf-8") as handle:
        handle.write(f"# nx {bx.shape[1]}\n")
        handle.write(f"# ny {bx.shape[0]}\n")
        handle.write("# columns: Bx By Bz\n")
        for iy in range(bx.shape[0]):
            for ix in range(bx.shape[1]):
                handle.write(f"{bx[iy, ix]: .9e} {by[iy, ix]: .9e} {bz[iy, ix]: .9e}\n")


def plot_quicklook(path, br, bt, bp, bx, by, bz, *, vmax: float):
    fig, axes = plt.subplots(2, 3, figsize=(12, 7), dpi=160, constrained_layout=True)
    panels = [
        (br, "Br"),
        (bt, "Btheta"),
        (bp, "Bphi"),
        (bx, "Bx = Bphi"),
        (by, "By = -Btheta"),
        (bz, "Bz = Br"),
    ]
    for ax, (data, title) in zip(axes.ravel(), panels, strict=True):
        im = ax.imshow(data, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(title)
        ax.set_axis_off()
        fig.colorbar(im, ax=ax, shrink=0.72)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

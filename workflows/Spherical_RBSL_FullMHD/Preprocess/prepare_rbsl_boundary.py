#!/usr/bin/env python
"""Prepare spherical RBSL/PFSS boundary products from HMI vector data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from scipy.ndimage import map_coordinates
from sunpy.coordinates import frames
from sunpy.map import Map

try:
    from .rbsl_core.maps import get_submap
    from .rbsl_core.vectors import plot_signed_vectors
    from .rbsl_core.spherical import (
        hmi_disambiguate_azimuth,
        hmi_vector_components,
        image_magnetic_to_spherical_components,
        lon_lat_centers,
        match_template_location,
        patch_bounds_from_lon_lat,
        pfss_amrvac_parameters,
        pfss_mod_usr_parameters,
        preprocess_spherical_field,
        refined_ghost_cells,
        replace_synoptic_patch,
        resample_mrdailysynframe,
        spherical_flux_balance,
        write_cutout,
        write_pfss_boundary,
    )
except ImportError:  # Direct execution: python Preprocess/prepare_rbsl_boundary.py
    from rbsl_core.maps import get_submap
    from rbsl_core.vectors import plot_signed_vectors
    from rbsl_core.spherical import (
        hmi_disambiguate_azimuth,
        hmi_vector_components,
        image_magnetic_to_spherical_components,
        lon_lat_centers,
        match_template_location,
        patch_bounds_from_lon_lat,
        pfss_amrvac_parameters,
        pfss_mod_usr_parameters,
        preprocess_spherical_field,
        refined_ghost_cells,
        replace_synoptic_patch,
        resample_mrdailysynframe,
        spherical_flux_balance,
        write_cutout,
        write_pfss_boundary,
    )


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    nghost = refined_ghost_cells(args.nghost, args.refine_level)
    print(
        f"Ghost cells: base={args.nghost}, refinement level={args.refine_level}, "
        f"effective={nghost}"
    )

    print("Loading HMI vector maps...")
    map_bx, map_by, map_bz = load_hmi_vector_maps(
        args.field,
        args.inclination,
        args.azimuth,
        args.disambig,
        disambig_method=args.disambig_method,
    )
    context_extent = compute_context_extent(
        map_bz,
        extent_arcsec=args.extent_arcsec,
        lon_range=args.lon_range,
        lat_range=args.lat_range,
    )
    print("Building and saving HMI Br context submap FITS...")
    map_br = build_hmi_br_map(map_bx, map_by, map_bz, extent_arcsec=context_extent)
    write_map_fits(map_br, output_dir / "hmi_br_submap.fits")

    print("Plotting context and simulation region...")
    plot_region_context(
        map_br,
        lon_range=args.lon_range,
        lat_range=args.lat_range,
        output=output_dir / "cut.png",
        vmax=args.context_vmax,
    )

    print("Extracting Br/Btheta/Bphi on the spherical grid...")
    extracted = extract_spherical_patch(
        map_bx,
        map_by,
        map_bz,
        lon_range=args.lon_range,
        lat_range=args.lat_range,
        resolution_degree=args.resolution_degree,
        nghost=nghost,
    )
    br_g, bt_g, bp_g = extracted["br_g"], extracted["bt_g"], extracted["bp_g"]
    br, bt, bp = extracted["br"], extracted["bt"], extracted["bp"]
    lon_g, lat_g = extracted["lon_g"], extracted["lat_g"]
    lon, lat = extracted["lon"], extracted["lat"]

    flux = spherical_flux_balance(br, lat)
    print_flux(flux)
    print(f"nlon: {br.shape[1]}  nlat: {br.shape[0]}")
    plot_projection(
        br,
        bp,
        -bt,
        output=output_dir / "projection.png",
        interval=args.vector_interval,
        scale=args.vector_scale,
        vmax=args.vector_vmax,
        threshold=args.vector_threshold,
    )

    if args.preprocess:
        print("Preprocessing vector field to reduce force/torque and smooth the boundary...")
        theta_g = np.radians(90.0 - lat_g)
        phi_g = np.radians(lon_g + 180.0)
        brn_g, btn_g, bpn_g = preprocess_spherical_field(
            br_g,
            bt_g,
            bp_g,
            theta_g,
            phi_g,
            mu3=args.preprocess_mu3,
            mu4=args.preprocess_mu4,
            max_iter=args.preprocess_max_iter,
            dt=args.preprocess_report_interval,
            tol=args.preprocess_tol,
            nghost=nghost,
            progress=print_preprocess_progress,
        )
    else:
        print("Skipping preprocessing; using extracted field directly.")
        brn_g, btn_g, bpn_g = br_g, bt_g, bp_g

    if nghost:
        core = (slice(nghost, -nghost), slice(nghost, -nghost))
        brn, btn, bpn = brn_g[core], btn_g[core], bpn_g[core]
    else:
        brn, btn, bpn = brn_g, btn_g, bpn_g

    print("Building pfss_boundary.dat from mrdailysynframe with the Br patch inserted...")
    synoptic = read_synoptic_data(args.synoptic)
    global_br, global_lon, global_lat = resample_mrdailysynframe(synoptic, args.resolution_degree)
    x0, y0, nx, ny = patch_bounds_from_lon_lat(args.lon_range, args.lat_range, args.resolution_degree)
    if (ny, nx) != brn.shape:
        raise ValueError(f"patch shape mismatch: expected {(ny, nx)} from ranges, got {brn.shape}")
    match_x, match_y, match_score = match_template_location(global_br, br)
    roll_shift = x0 - match_x
    print(
        "synoptic/template match: "
        f"x={match_x}, y={match_y}, score={match_score:.5f}, "
        f"roll_shift={roll_shift}"
    )
    matched_global_br = np.roll(global_br, shift=roll_shift, axis=1)
    pfss_br = replace_synoptic_patch(matched_global_br, brn, x0, y0)
    write_pfss_boundary(output_dir / "pfss_boundary.dat", pfss_br, global_lon, global_lat, scale=args.b_scale)
    plot_synoptic_patch(
        matched_global_br,
        pfss_br,
        x0=x0,
        y0=y0,
        nx=nx,
        ny=ny,
        output=output_dir / "synoptic.png",
    )

    print("Writing local cutout.dat...")
    write_cutout(output_dir / args.cutout_name, brn_g, btn_g, bpn_g, scale=args.b_scale)

    mod_usr = pfss_mod_usr_parameters(pfss_br.shape, x0=x0, y0=y0, nx=nx, ny=ny)
    amrvac = pfss_amrvac_parameters(
        pfss_br.shape,
        x0=x0,
        y0=y0,
        nx=nx,
        ny=ny,
        nr=args.nr,
        length_r=args.length_r,
        refine_level=args.refine_level,
    )
    parameters = {
        "schema_version": 1,
        "workflow": "Spherical_RBSL_FullMHD",
        "coordinate_convention": "Stonyhurst longitude plus 180 degrees in simulation phi; theta is colatitude",
        "reference_time": map_bz.date.isot,
        "inputs": {
            "field": str(Path(args.field).expanduser().resolve()),
            "inclination": str(Path(args.inclination).expanduser().resolve()),
            "azimuth": str(Path(args.azimuth).expanduser().resolve()),
            "disambig": str(Path(args.disambig).expanduser().resolve()),
            "synoptic": str(Path(args.synoptic).expanduser().resolve()),
        },
        "pfss_mod_usr.t": mod_usr,
        "pfss_amrvac.par": amrvac,
        "flux": flux,
        "grid": {
            "resolution_degree": args.resolution_degree,
            "base_nghost": args.nghost,
            "refine_level": args.refine_level,
            "effective_nghost": nghost,
            "lon_range": args.lon_range,
            "lat_range": args.lat_range,
            "context_extent_arcsec": list(context_extent),
            "nlon": nx,
            "nlat": ny,
            "x0": x0,
            "y0": y0,
            "match_x": match_x,
            "match_y": match_y,
            "match_score": match_score,
            "roll_shift": roll_shift,
        },
        "outputs": {
            "hmi_br_submap_fits": str(output_dir / "hmi_br_submap.fits"),
            "pfss_boundary": str(output_dir / "pfss_boundary.dat"),
            "cutout": str(output_dir / args.cutout_name),
        },
    }
    (output_dir / "parameters.json").write_text(json.dumps(parameters, indent=2), encoding="utf-8")
    (output_dir / "parameters.txt").write_text(format_parameters(mod_usr, amrvac), encoding="utf-8")
    print(format_parameters(mod_usr, amrvac))
    print(f"Saved outputs in {output_dir}")
    return 0


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--field", required=True, help="HMI vector field FITS segment.")
    parser.add_argument("--inclination", required=True, help="HMI inclination FITS segment.")
    parser.add_argument("--azimuth", required=True, help="HMI azimuth FITS segment.")
    parser.add_argument("--disambig", required=True, help="HMI disambiguation FITS segment.")
    parser.add_argument("--synoptic", required=True, help="HMI mrdailysynframe_polfil Mr_polfil FITS file.")
    parser.add_argument("--output-dir", required=True, help="Directory for figures and boundary files.")
    parser.add_argument("--extent-arcsec", nargs=4, type=float, metavar=("X0", "X1", "Y0", "Y1"), required=True)
    parser.add_argument("--lon-range", nargs=2, type=float, metavar=("LON0", "LON1"), required=True)
    parser.add_argument("--lat-range", nargs=2, type=float, metavar=("LAT0", "LAT1"), required=True)
    parser.add_argument("--resolution-degree", type=float, default=0.05)
    parser.add_argument(
        "--nghost",
        type=int,
        default=2,
        help="Base ghost-cell count; effective width is nghost*2**(refine_level-1), matching IDL dixB.",
    )
    parser.add_argument("--nr", type=int, default=400)
    parser.add_argument("--length-r", type=float, default=0.8)
    parser.add_argument("--refine-level", type=int, default=1)
    parser.add_argument("--b-scale", type=float, default=1.0)
    parser.add_argument("--disambig-method", type=int, default=1)
    parser.add_argument("--cutout-name", default="cutout.dat")
    parser.add_argument("--context-vmax", type=float, default=300.0)
    parser.add_argument("--vector-interval", type=int, default=4)
    parser.add_argument("--vector-scale", type=float, default=50000.0)
    parser.add_argument("--vector-vmax", type=float, default=1000.0)
    parser.add_argument("--vector-threshold", type=float, default=100.0)
    parser.add_argument("--preprocess", dest="preprocess", action="store_true", default=True)
    parser.add_argument("--no-preprocess", dest="preprocess", action="store_false")
    parser.add_argument("--preprocess-mu3", type=float, default=0.1)
    parser.add_argument("--preprocess-mu4", type=float, default=1.0)
    parser.add_argument("--preprocess-max-iter", type=int, default=5000)
    parser.add_argument("--preprocess-report-interval", type=int, default=200)
    parser.add_argument("--preprocess-tol", type=float, default=1e-4)
    return parser.parse_args(argv)


def load_hmi_vector_maps(field_path, inclination_path, azimuth_path, disambig_path, *, disambig_method: int):
    field = Map(field_path).rotate()
    inclination = Map(inclination_path).rotate()
    azimuth = Map(azimuth_path).rotate()
    disambig = Map(disambig_path).rotate(missing=0)
    corrected_azimuth = hmi_disambiguate_azimuth(azimuth.data, disambig.data, method=disambig_method)
    bx, by, bz = hmi_vector_components(field.data, inclination.data, corrected_azimuth)
    return Map(bx, field.meta), Map(by, field.meta), Map(bz, field.meta)


def build_hmi_br_map(map_bx, map_by, map_bz, *, extent_arcsec=None):
    if extent_arcsec is not None:
        map_bx = get_submap(map_bx, tuple(extent_arcsec))
        map_by = get_submap(map_by, tuple(extent_arcsec))
        map_bz = get_submap(map_bz, tuple(extent_arcsec))
    ypix, xpix = np.indices(map_bz.data.shape, dtype=float)
    coords_hpc = map_bz.pixel_to_world(xpix * u.pix, ypix * u.pix)
    _, _, _, br, _, _ = image_magnetic_to_spherical_components(
        coords_hpc.Tx.value.ravel(),
        coords_hpc.Ty.value.ravel(),
        map_bx.data.ravel(),
        map_by.data.ravel(),
        map_bz.data.ravel(),
        map_bz.observer_coordinate.lat.value - 90.0,
        map_bz.rsun_obs.value,
    )
    meta = map_bz.meta.copy()
    meta["content"] = "BR"
    meta["bunit"] = "G"
    meta["SOLBRMAP"] = 1
    meta.pop("BLANK", None)
    return Map(br.reshape(map_bz.data.shape).astype(np.float32), meta)


def write_map_fits(smap, output):
    smap.save(output, overwrite=True)


def compute_context_extent(smap, *, extent_arcsec, lon_range, lat_range):
    region_corners = [
        stonyhurst_coord(lon_range[0], lat_range[0], smap),
        stonyhurst_coord(lon_range[0], lat_range[1], smap),
        stonyhurst_coord(lon_range[1], lat_range[0], smap),
        stonyhurst_coord(lon_range[1], lat_range[1], smap),
    ]
    return expand_extent_for_coords(tuple(extent_arcsec), region_corners, smap)


def plot_region_context(smap, *, lon_range, lat_range, output, vmax: float):
    region_corners = [
        stonyhurst_coord(lon_range[0], lat_range[0], smap),
        stonyhurst_coord(lon_range[0], lat_range[1], smap),
        stonyhurst_coord(lon_range[1], lat_range[0], smap),
        stonyhurst_coord(lon_range[1], lat_range[1], smap),
    ]
    bottom_left = region_corners[0]
    top_right = region_corners[-1]
    fig = plt.figure(figsize=(7.0, 5.5), dpi=180)
    ax = fig.add_subplot(projection=smap)
    smap.plot(axes=ax, norm=Normalize(-vmax, vmax), cmap="gray")
    ax.coords.grid(False)
    title = ax.get_title()
    ax.set_title("")
    ax.coords[0].set_auto_axislabel(False)
    ax.coords[1].set_auto_axislabel(False)
    ax.set_xlabel("Solar-X")
    ax.set_ylabel("Solar-Y")
    smap.draw_grid(
        axes=ax,
        grid_spacing=(10, 10) * u.deg,
        annotate=True,
        color="yellow",
        alpha=0.45,
        linewidth=0.8,
    )
    if title:
        fig.suptitle(title, y=0.99)
        fig.subplots_adjust(top=0.90)
    smap.draw_quadrangle(bottom_left, top_right=top_right, axes=ax, edgecolor="cyan", linewidth=1.5)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def expand_extent_for_coords(extent_arcsec, coords, smap, *, margin_arcsec: float | None = None):
    x_min, x_max, y_min, y_max = map(float, extent_arcsec)
    if margin_arcsec is None:
        margin_arcsec = max(20.0, 0.05 * max(x_max - x_min, y_max - y_min))

    xs = []
    ys = []
    for coord in coords:
        projected = coord.transform_to(smap.coordinate_frame)
        x = projected.Tx.to_value(u.arcsec)
        y = projected.Ty.to_value(u.arcsec)
        if np.isfinite(x) and np.isfinite(y):
            xs.append(x)
            ys.append(y)

    if not xs or not ys:
        return (x_min, x_max, y_min, y_max)

    map_x_min = smap.bottom_left_coord.Tx.to_value(u.arcsec)
    map_x_max = smap.top_right_coord.Tx.to_value(u.arcsec)
    map_y_min = smap.bottom_left_coord.Ty.to_value(u.arcsec)
    map_y_max = smap.top_right_coord.Ty.to_value(u.arcsec)

    expanded = (
        max(map_x_min, min(x_min, min(xs) - margin_arcsec)),
        min(map_x_max, max(x_max, max(xs) + margin_arcsec)),
        max(map_y_min, min(y_min, min(ys) - margin_arcsec)),
        min(map_y_max, max(y_max, max(ys) + margin_arcsec)),
    )
    return expanded


def extract_spherical_patch(map_bx, map_by, map_bz, *, lon_range, lat_range, resolution_degree, nghost):
    lon_g, lat_g = lon_lat_centers(tuple(lon_range), tuple(lat_range), resolution_degree, nghost=nghost)
    lon2, lat2 = np.meshgrid(lon_g, lat_g)
    coords_car = SkyCoord(
        lon2.ravel() * u.degree,
        lat2.ravel() * u.degree,
        frame=frames.HeliographicStonyhurst,
        obstime=map_bz.date,
        rsun=map_bz.rsun_meters,
        observer="earth",
    )
    coords_arc = coords_car.transform_to(frame=map_bz.coordinate_frame)
    coords_pixel = map_bz.world_to_pixel(coords_arc)
    coords_pixel = np.array([coords_pixel[1].value, coords_pixel[0].value])
    bz_select = map_coordinates(map_bz.data, coords_pixel, order=1)
    by_select = map_coordinates(map_by.data, coords_pixel, order=1)
    bx_select = map_coordinates(map_bx.data, coords_pixel, order=1)
    _, _, _, br_g, bt_g, bp_g = image_magnetic_to_spherical_components(
        coords_arc.Tx.value,
        coords_arc.Ty.value,
        bx_select,
        by_select,
        bz_select,
        map_bz.observer_coordinate.lat.value - 90.0,
        map_bz.rsun_obs.value,
    )
    shape = lon2.shape
    br_g = br_g.reshape(shape)
    bt_g = bt_g.reshape(shape)
    bp_g = bp_g.reshape(shape)
    if nghost:
        core = (slice(nghost, -nghost), slice(nghost, -nghost))
        lon = lon_g[nghost:-nghost]
        lat = lat_g[nghost:-nghost]
    else:
        core = (slice(None), slice(None))
        lon = lon_g
        lat = lat_g
    return {
        "lon_g": lon_g,
        "lat_g": lat_g,
        "br_g": br_g,
        "bt_g": bt_g,
        "bp_g": bp_g,
        "lon": lon,
        "lat": lat,
        "br": br_g[core],
        "bt": bt_g[core],
        "bp": bp_g[core],
    }


def plot_projection(br, vx, vy, *, output, interval, scale, vmax, threshold):
    fig, ax = plt.subplots(figsize=(7.0, 5.5), dpi=180)
    plot_signed_vectors(br, vx, vy, interval=interval, scale=scale, ax=ax, bzmax=vmax, threshold=threshold)
    ax.set_axis_off()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def read_synoptic_data(path):
    with fits.open(path) as hdul:
        if len(hdul) > 1 and hdul[1].data is not None:
            return np.asarray(hdul[1].data, dtype=float)
        return np.asarray(hdul[0].data, dtype=float)


def plot_synoptic_patch(before, after, *, x0, y0, nx, ny, output):
    pad = 50
    y_min = max(0, y0 - pad)
    y_max = min(before.shape[0], y0 + ny + pad)
    x_min = max(0, x0 - pad)
    x_max = min(before.shape[1], x0 + nx + pad)
    fig, axes = plt.subplots(2, 1, figsize=(10, 9), dpi=160)
    for ax, data, title in zip(axes, [before, after], ["Synoptic before replacement", "PFSS boundary"], strict=True):
        ax.imshow(data[y_min:y_max, x_min:x_max], origin="lower", cmap="gray", vmax=300, vmin=-300)
        ax.set_title(title)
        ax.set_axis_off()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def stonyhurst_coord(lon, lat, smap):
    return SkyCoord(
        lon * u.degree,
        lat * u.degree,
        frame=frames.HeliographicStonyhurst,
        obstime=smap.date,
        rsun=smap.rsun_meters,
        observer="earth",
    )


def print_flux(flux):
    print(
        "flux balance: "
        f"signed={flux['signed']:.6g}, unsigned={flux['unsigned']:.6g}, "
        f"positive={flux['positive']:.6g}, negative={flux['negative']:.6g}, "
        f"balance={flux['balance']:.6g}"
    )


def print_preprocess_progress(values):
    print(
        f"it={values['iteration']:5d}  dL={values['dloss']:.3e}  "
        f"L12={values['l12']:.3e}  L3={values['l3']:.3e}  "
        f"L4={values['l4']:.3e}  L={values['loss']:.3e}  "
        f"eps_f={values['eps_f']:.2e}  eps_t={values['eps_t']:.2e}"
    )


def format_parameters(mod_usr, amrvac):
    lines = ["parameters should be changed in pfss_mod_usr.t"]
    lines.append(f"lat 0: {mod_usr['lat0']} lat1: {mod_usr['lat1']}")
    lines.append(f"lon0: {mod_usr['lon0']} lon1: {mod_usr['lon1']}")
    lines.append("")
    lines.append("parameters should be changed in pfss_amrvac.par")
    for key, value in amrvac.items():
        if key.startswith("xprob"):
            lines.append(f"        {key}= {value}d0")
        else:
            lines.append(f"        {key}= {value}")
    return "\n".join(lines) + "\n"


if __name__ == "__main__":
    raise SystemExit(main())

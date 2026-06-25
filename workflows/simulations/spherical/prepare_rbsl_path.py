#!/usr/bin/env python
"""Interactively pick an RBSL axis and write `path.dat`."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import plotly.graph_objects as go
from astropy.coordinates import SkyCoord
from matplotlib.patches import Circle
from scipy.interpolate import splev, splprep
from scipy.ndimage import map_coordinates
from sunpy.coordinates import frames
from sunpy.map import Map

from solarutils.obs.maps import get_submap
from solarutils.sim.spherical import (
    flux_rope_axis_from_points,
    flux_rope_boundary_fieldline,
    hmi_disambiguate_azimuth,
    hmi_vector_components,
    image_magnetic_to_spherical_components,
    lon_lat_centers,
    read_path_dat,
    write_path_dat,
)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    image_clip_percent = resolve_image_clip_percent(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    aia_map = Map(args.aia)
    submap = get_submap(aia_map, tuple(args.submap_arcsec))
    reference_path_coords = None
    if args.reference_points:
        reference_points = np.load(args.reference_points)
        reference_path_coords = [stonyhurst_to_hpc(lon, lat, aia_map) for lon, lat in reference_points]

    if args.points:
        points = np.load(args.points)
    elif args.path_points_stonyhurst:
        points = np.asarray(args.path_points_stonyhurst, dtype=float).reshape(-1, 2)
    else:
        points, _ = pick_path_points_on_submap(
            submap,
            output_dir / "picked_points.png",
            circle_radius_arcsec=args.footpoint_circle_arcsec,
            smooth_points=args.smooth_path_points,
            smoothing=args.smooth_path_smoothing,
            reference_path_coords=reference_path_coords,
        )

    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or points.shape[0] < 2:
        raise ValueError("path points must have shape (N, 2) with at least two points")
    footpoints = points[[0, -1]]
    footpoints_hpc = [stonyhurst_to_hpc(lon, lat, aia_map) for lon, lat in footpoints]
    plot_footpoints_on_submap(
        submap,
        footpoints_hpc,
        output_dir / "selected_footpoints.png",
        circle_radius_arcsec=args.footpoint_circle_arcsec,
    )

    np.save(output_dir / "pts.npy", points)
    if args.points or args.path_points_stonyhurst:
        plot_path_on_submap(
            submap,
            [stonyhurst_to_hpc(lon, lat, aia_map) for lon, lat in points],
            output_dir / "picked_points.png",
            circle_radius_arcsec=args.footpoint_circle_arcsec,
            reference_path_coords=reference_path_coords,
        )

    lon_range = tuple(args.lon_range) if args.lon_range is not None else None
    lat_range = tuple(args.lat_range) if args.lat_range is not None else None
    lonlat_range_inferred = False
    if args.write_3d and lon_range is None and lat_range is None:
        lon_range, lat_range = infer_lonlat_range(points)
        lonlat_range_inferred = True
    needs_lonlat_grid = lon_range is not None or lat_range is not None or args.hmi_field or args.write_3d
    if needs_lonlat_grid and (lon_range is None or lat_range is None):
        raise ValueError("--lon-range and --lat-range must be provided together; 3D-only runs can omit both")
    lon = lat = base_data = extent = None
    if needs_lonlat_grid:
        lon, lat = lon_lat_centers(lon_range, lat_range, args.resolution_degree, nghost=args.nghost)
        base_path = args.base_map or args.base_aia
        base_data = sample_map_on_stonyhurst(aia_map if base_path is None else Map(base_path), lon, lat)
        extent = [lon[0], lon[-1], lat[0], lat[-1]]
        plot_lonlat_base(
            base_data,
            extent,
            footpoints,
            output_dir / "base_footpoints.png",
            circle_radius_degree=args.footpoint_circle_degree,
            cmap=args.cmap,
            image_scale=args.image_scale,
            image_clip_percent=image_clip_percent,
            image_power=args.image_power,
            image_vmin=args.image_vmin,
            image_vmax=args.image_vmax,
        )
        plot_lonlat_base(
            base_data,
            extent,
            footpoints,
            output_dir / "picked_points_lonlat.png",
            circle_radius_degree=args.footpoint_circle_degree,
            cmap=args.cmap,
            path_points=points,
            image_scale=args.image_scale,
            image_clip_percent=image_clip_percent,
            image_power=args.image_power,
            image_vmin=args.image_vmin,
            image_vmax=args.image_vmax,
        )

    hmi_br = None
    if args.hmi_field:
        hmi_br = load_hmi_br_on_stonyhurst(
            args.hmi_field,
            args.hmi_inclination,
            args.hmi_azimuth,
            args.hmi_disambig,
            lon,
            lat,
            disambig_method=args.hmi_disambig_method,
        )
        plot_hmi_lonlat_base(
            hmi_br,
            extent,
            footpoints,
            output_dir / "hmi_footpoints.png",
            circle_radius_degree=args.footpoint_circle_degree,
            path_points=None,
            vmin=args.hmi_vmin,
            vmax=args.hmi_vmax,
        )
        plot_hmi_lonlat_base(
            hmi_br,
            extent,
            footpoints,
            output_dir / "hmi_picked_points.png",
            circle_radius_degree=args.footpoint_circle_degree,
            path_points=points,
            vmin=args.hmi_vmin,
            vmax=args.hmi_vmax,
        )

    solar_radius_mm = args.solar_radius_mm
    if solar_radius_mm is None:
        solar_radius_mm = aia_map.rsun_meters.to_value(u.Mm)
    result = flux_rope_axis_from_points(points, args.height_mm, npoints=args.npoints, solar_radius_mm=solar_radius_mm)
    boundary_path = None
    if args.rope_radius_mm and args.rope_radius_mm > 0.0 and args.rope_twist:
        boundary_path = flux_rope_boundary_fieldline(
            result["path"],
            args.rope_radius_mm,
            args.rope_twist,
            solar_radius_mm=solar_radius_mm,
        )
    path = write_path_dat(output_dir / "path.dat", result["path"], add_lon_180=True)
    np.savetxt(output_dir / "path_lonlat_radius.txt", result["path"], fmt="%25.16e")
    if boundary_path is not None:
        np.savetxt(output_dir / "boundary_lonlat_radius.txt", boundary_path, fmt="%25.16e")
    else:
        (output_dir / "boundary_lonlat_radius.txt").unlink(missing_ok=True)
    plot_projected_path_on_aia(
        submap,
        result["path"],
        output_dir / "projected_path_aia.png",
        circle_radius_arcsec=args.footpoint_circle_arcsec,
        boundary_lonlatr=boundary_path,
    )
    rope_3d_html = None
    rope_3d_hmi_html = None
    if args.write_3d:
        surface_data, surface_vmin, surface_vmax, _ = image_display_data(
            base_data,
            image_scale=args.image_scale,
            clip_percent=image_clip_percent,
            image_power=args.image_power,
            image_vmin=args.image_vmin,
            image_vmax=args.image_vmax,
        )
        rope_3d_html = output_dir / "rope_3d.html"
        plot_3d_rope_surface(
            result["path"],
            lon,
            lat,
            surface_data,
            rope_3d_html,
            solar_radius_mm=solar_radius_mm,
            surface_range=(surface_vmin, surface_vmax),
            surface_colorscale=plotly_colorscale_for_cmap(args.cmap),
            title="3D flux rope on AIA lon/lat map",
            rope_width=args.rope_3d_width,
        )
        if hmi_br is not None:
            rope_3d_hmi_html = output_dir / "rope_3d_hmi.html"
            plot_3d_rope_surface(
                result["path"],
                lon,
                lat,
                hmi_br,
                rope_3d_hmi_html,
                solar_radius_mm=solar_radius_mm,
                surface_range=(args.hmi_vmin, args.hmi_vmax),
                surface_colorscale="Gray",
                title="3D flux rope on HMI Br map",
                rope_width=args.rope_3d_width,
            )

    summary = {
        "height_mm": args.height_mm,
        "solar_radius_mm": solar_radius_mm,
        "npoints": args.npoints,
        "footpoints_hpc_arcsec": [[coord.Tx.to_value(u.arcsec), coord.Ty.to_value(u.arcsec)] for coord in footpoints_hpc],
        "footpoints_stonyhurst": footpoints.tolist(),
        "picked_points_count": int(points.shape[0]),
        "rope_radius_mm": args.rope_radius_mm,
        "rope_twist": args.rope_twist,
        "hmi_br_plotted": hmi_br is not None,
        "rope_3d_html": str(rope_3d_html) if rope_3d_html is not None else None,
        "rope_3d_hmi_html": str(rope_3d_hmi_html) if rope_3d_hmi_html is not None else None,
        "path_dat": str(path),
        "lon_range": list(lon_range) if lon_range is not None else None,
        "lat_range": list(lat_range) if lat_range is not None else None,
        "lonlat_range_inferred": lonlat_range_inferred,
        "image_scale": args.image_scale,
        "image_clip_percent": list(image_clip_percent),
        "depth_mm": result["depth_mm"],
        "major_radius_mm": result["major_radius_mm"],
        "endpoint_angle_rad": result["endpoint_angle_rad"],
        "projected_length_mm": result["projected_length_mm"],
    }
    (output_dir / "path_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

    print(f"footpoints lon/lat: {footpoints.tolist()}")
    print(f"picked points: {points.shape[0]}")
    print(f"depth = {result['depth_mm']:.3f} Mm")
    print(f"major radius = {result['major_radius_mm']:.3f} Mm")
    print(f"projected length = {result['projected_length_mm']:.3f} Mm")
    print(f"wrote: {path}")
    return 0


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aia", required=True, help="AIA FITS used for the path-picking submap.")
    parser.add_argument("--base-map", help="Optional AIA/HMI FITS used as lon/lat path-picking background.")
    parser.add_argument("--base-aia", help=argparse.SUPPRESS)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--submap-arcsec", nargs=4, type=float, metavar=("X0", "X1", "Y0", "Y1"), required=True)
    parser.add_argument("--lon-range", nargs=2, type=float, metavar=("LON0", "LON1"), help="Optional lon/lat map longitude range.")
    parser.add_argument("--lat-range", nargs=2, type=float, metavar=("LAT0", "LAT1"), help="Optional lon/lat map latitude range.")
    parser.add_argument("--resolution-degree", type=float, default=0.1)
    parser.add_argument("--nghost", type=int, default=2)
    parser.add_argument("--height-mm", type=float, required=True)
    parser.add_argument("--solar-radius-mm", type=float)
    parser.add_argument("--npoints", type=int, default=6000)
    parser.add_argument("--footpoints-arcsec", nargs=4, type=float, metavar=("X0", "Y0", "X1", "Y1"), help=argparse.SUPPRESS)
    parser.add_argument(
        "--footpoints-stonyhurst",
        nargs=4,
        type=float,
        metavar=("LON0", "LAT0", "LON1", "LAT1"),
        help=argparse.SUPPRESS,
    )
    parser.add_argument("--path-points-stonyhurst", nargs="+", type=float, help="Flat lon lat lon lat ... list.")
    parser.add_argument("--points", help="Existing pts.npy with lon/lat path points.")
    parser.add_argument("--reference-points", help="Optional pts.npy to overlay as a reference curve while picking.")
    parser.add_argument(
        "--smooth-path-points",
        type=int,
        default=0,
        help="Number of smoothed path points to generate from interactive clicks; 0 keeps raw clicks.",
    )
    parser.add_argument(
        "--smooth-path-smoothing",
        type=float,
        default=0.0,
        help="Spline smoothing factor in clicked pixel coordinates.",
    )
    parser.add_argument("--footpoint-circle-arcsec", type=float, default=20.0)
    parser.add_argument("--footpoint-circle-degree", type=float, default=20.0 / (960.0 * np.pi) * 180.0)
    parser.add_argument("--cmap", default="sdoaia304")
    parser.add_argument("--image-scale", choices=["log", "power"], default="log", help="Lon/lat AIA display scaling.")
    parser.add_argument("--image-clip-percent", nargs=2, type=float, default=[1.0, 99.5], metavar=("LOW", "HIGH"))
    parser.add_argument("--image-vmin-percent", type=float, help="Lower display percentile; overrides --image-clip-percent LOW.")
    parser.add_argument("--image-vmax-percent", type=float, help="Upper display percentile; overrides --image-clip-percent HIGH.")
    parser.add_argument("--image-power", type=float, help="Legacy power applied to abs(base data) before plotting.")
    parser.add_argument("--image-vmin", type=float, help="Legacy lower display limit after image-power.")
    parser.add_argument("--image-vmax", type=float, help="Legacy upper display limit after image-power.")
    parser.add_argument("--hmi-field", help="Optional HMI field FITS for Br marker plots.")
    parser.add_argument("--hmi-inclination", help="HMI inclination FITS, required with --hmi-field.")
    parser.add_argument("--hmi-azimuth", help="HMI azimuth FITS, required with --hmi-field.")
    parser.add_argument("--hmi-disambig", help="HMI disambiguation FITS, required with --hmi-field.")
    parser.add_argument("--hmi-disambig-method", type=int, default=2)
    parser.add_argument("--hmi-vmin", type=float, default=-300.0)
    parser.add_argument("--hmi-vmax", type=float, default=300.0)
    parser.add_argument("--write-3d", action="store_true", help="Write rope_3d HTML files; disabled by default to save disk space.")
    parser.add_argument("--rope-3d-width", type=float, default=10.0)
    parser.add_argument("--rope-radius-mm", type=float, default=0.0, help="Optional radius for a helical boundary field line in Mm.")
    parser.add_argument("--rope-twist", type=float, default=0.0, help="Number of turns for the optional boundary field line.")
    return parser.parse_args(argv)


def resolve_image_clip_percent(args):
    low, high = args.image_clip_percent
    if args.image_vmin_percent is not None:
        low = args.image_vmin_percent
    if args.image_vmax_percent is not None:
        high = args.image_vmax_percent
    if not (0.0 <= low < high <= 100.0):
        raise ValueError("image display percentiles must satisfy 0 <= vmin < vmax <= 100")
    return float(low), float(high)


def infer_lonlat_range(points, *, padding_fraction=0.15, min_padding_degree=5.0):
    points = np.asarray(points, dtype=float)
    lon_min, lat_min = np.nanmin(points, axis=0)
    lon_max, lat_max = np.nanmax(points, axis=0)
    lon_span = lon_max - lon_min
    lat_span = lat_max - lat_min
    lon_padding = max(lon_span * padding_fraction, min_padding_degree)
    lat_padding = max(lat_span * padding_fraction, min_padding_degree)
    return (float(lon_min - lon_padding), float(lon_max + lon_padding)), (
        float(lat_min - lat_padding),
        float(lat_max + lat_padding),
    )


def pick_two_hpc_points(smap, output):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection=smap)
    smap.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
    set_solar_xy_labels(ax)
    ax.set_title("Click two footpoints, then close the window")
    clicks = plt.ginput(2, timeout=0)
    if len(clicks) != 2:
        raise RuntimeError("expected exactly two clicked footpoints")
    coords = [smap.pixel_to_world(x * u.pix, y * u.pix) for x, y in clicks]
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return coords


def plot_footpoints_on_submap(smap, coords, output, *, circle_radius_arcsec):
    fig = plt.figure(figsize=(8, 8), dpi=160)
    ax = fig.add_subplot(projection=smap)
    smap.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
    set_solar_xy_labels(ax)
    colors = ["white", "black"]
    for coord, color in zip(coords, colors, strict=True):
        ax.plot_coord(
            coord,
            "o",
            markersize=12,
            markerfacecolor=color,
            markeredgecolor=color,
            markeredgewidth=2,
            alpha=0.45,
        )
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def sample_map_on_stonyhurst(smap, lon, lat):
    lon2, lat2 = np.meshgrid(lon, lat)
    coords_hgs = SkyCoord(
        lon2.ravel() * u.degree,
        lat2.ravel() * u.degree,
        frame=frames.HeliographicStonyhurst,
        obstime=smap.date,
        rsun=smap.rsun_meters,
        observer="earth",
    )
    coords_hpc = coords_hgs.transform_to(smap.coordinate_frame)
    pixels = smap.world_to_pixel(coords_hpc)
    pixel_coords = np.array([pixels[1].value, pixels[0].value])
    data = map_coordinates(np.asarray(smap.data, dtype=float), pixel_coords, order=1, mode="constant", cval=np.nan)
    return data.reshape(lat2.shape)


def load_hmi_br_on_stonyhurst(field_path, inclination_path, azimuth_path, disambig_path, lon, lat, *, disambig_method):
    missing = [
        name
        for name, value in [
            ("--hmi-inclination", inclination_path),
            ("--hmi-azimuth", azimuth_path),
            ("--hmi-disambig", disambig_path),
        ]
        if value is None
    ]
    if missing:
        raise ValueError(f"{', '.join(missing)} required with --hmi-field")

    field = Map(field_path).rotate()
    inclination = Map(inclination_path).rotate()
    azimuth = Map(azimuth_path).rotate()
    disambig = Map(disambig_path).rotate(missing=0)
    corrected_azimuth = hmi_disambiguate_azimuth(azimuth.data, disambig.data, method=disambig_method)
    bx, by, bz = hmi_vector_components(field.data, inclination.data, corrected_azimuth)
    map_bx = Map(bx, field.meta)
    map_by = Map(by, field.meta)
    map_bz = Map(bz, field.meta)

    lon2, lat2 = np.meshgrid(lon, lat)
    coords_hgs = SkyCoord(
        lon2.ravel() * u.degree,
        lat2.ravel() * u.degree,
        frame=frames.HeliographicStonyhurst,
        obstime=map_bz.date,
        rsun=map_bz.rsun_meters,
        observer="earth",
    )
    coords_hpc = coords_hgs.transform_to(map_bz.coordinate_frame)
    pixels = map_bz.world_to_pixel(coords_hpc)
    pixel_coords = np.array([pixels[1].value, pixels[0].value])
    bx_select = map_coordinates(map_bx.data, pixel_coords, order=1, mode="constant", cval=np.nan)
    by_select = map_coordinates(map_by.data, pixel_coords, order=1, mode="constant", cval=np.nan)
    bz_select = map_coordinates(map_bz.data, pixel_coords, order=1, mode="constant", cval=np.nan)
    _, _, _, br, _, _ = image_magnetic_to_spherical_components(
        coords_hpc.Tx.value,
        coords_hpc.Ty.value,
        bx_select,
        by_select,
        bz_select,
        map_bz.observer_coordinate.lat.value - 90.0,
        map_bz.rsun_obs.value,
    )
    return br.reshape(lon2.shape)


def plot_lonlat_base(
    data,
    extent,
    footpoints,
    output,
    *,
    circle_radius_degree,
    cmap,
    image_scale,
    image_clip_percent,
    image_power,
    image_vmin,
    image_vmax,
    path_points=None,
):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=160)
    image, vmin, vmax, _ = image_display_data(
        data,
        image_scale=image_scale,
        clip_percent=image_clip_percent,
        image_power=image_power,
        image_vmin=image_vmin,
        image_vmax=image_vmax,
    )
    ax.imshow(
        image,
        origin="lower",
        extent=extent,
        cmap=cmap if cmap in plt.colormaps() else "gray",
        aspect="auto",
        vmin=vmin,
        vmax=vmax,
    )
    colors = ["white", "black"]
    for (lon, lat), color in zip(footpoints, colors, strict=True):
        ax.add_patch(Circle((lon, lat), circle_radius_degree, facecolor=color, edgecolor=color, linewidth=1.5, alpha=0.35))
    if path_points is not None:
        ax.plot(path_points[:, 0], path_points[:, 1], "o-", color="yellow", markersize=3, linewidth=1)
    ax.set_xlabel("Stonyhurst longitude [deg]")
    ax.set_ylabel("Stonyhurst latitude [deg]")
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def image_display_data(data, *, image_scale, clip_percent, image_power, image_vmin, image_vmax):
    data = np.asarray(data, dtype=float)
    if image_scale == "power" or any(value is not None for value in (image_power, image_vmin, image_vmax)):
        power = 0.25 if image_power is None else image_power
        image = np.abs(data) ** power
        finite = image[np.isfinite(image)]
        if finite.size == 0:
            return image, image_vmin, image_vmax, "power"
        vmin = float(np.nanpercentile(finite, clip_percent[0])) if image_vmin is None else image_vmin
        vmax = float(np.nanpercentile(finite, clip_percent[1])) if image_vmax is None else image_vmax
        return image, vmin, vmax, "power"

    positive = np.abs(data)
    finite_positive = positive[np.isfinite(positive) & (positive > 0.0)]
    if finite_positive.size == 0:
        image = np.full_like(positive, np.nan, dtype=float)
        return image, None, None, "log"
    floor = float(np.nanpercentile(finite_positive, 0.1))
    floor = max(floor, np.finfo(float).tiny)
    image = np.log10(np.clip(positive, floor, None))
    finite = image[np.isfinite(image)]
    vmin, vmax = np.nanpercentile(finite, clip_percent)
    return image, float(vmin), float(vmax), "log"


def plot_hmi_lonlat_base(data, extent, footpoints, output, *, circle_radius_degree, path_points, vmin, vmax):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=160)
    ax.imshow(data, origin="lower", extent=extent, cmap="gray", aspect="auto", vmin=vmin, vmax=vmax)
    for (lon, lat), color in zip(footpoints, ["white", "black"], strict=True):
        ax.add_patch(Circle((lon, lat), circle_radius_degree, facecolor=color, edgecolor=color, linewidth=1.5, alpha=0.35))
    if path_points is not None:
        ax.plot(path_points[:, 0], path_points[:, 1], "o-", color="yellow", markersize=3, linewidth=1)
    ax.set_xlabel("Stonyhurst longitude [deg]")
    ax.set_ylabel("Stonyhurst latitude [deg]")
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_3d_rope_surface(
    path,
    lon,
    lat,
    surface_data,
    output,
    *,
    solar_radius_mm,
    surface_range,
    surface_colorscale,
    title,
    rope_width,
):
    lon_path = path[:, 0]
    lat_path = path[:, 1]
    radius_path = path[:, 2] * solar_radius_mm
    x_path, y_path, z_path = spherical_surface_xyz(lon_path, lat_path, radius_path)

    lon2, lat2 = np.meshgrid(lon, lat)
    x_surface, y_surface, z_surface = spherical_surface_xyz(lon2, lat2, solar_radius_mm)
    surface = go.Surface(
        x=x_surface,
        y=y_surface,
        z=z_surface,
        surfacecolor=surface_data,
        colorscale=surface_colorscale,
        cmin=surface_range[0],
        cmax=surface_range[1],
        showscale=False,
        opacity=1.0,
        name="Surface map",
    )
    rope = go.Scatter3d(
        x=x_path,
        y=y_path,
        z=z_path,
        mode="lines",
        line={"color": "rgba(0, 0, 255, 0.65)", "width": rope_width},
        name="Flux rope axis",
    )
    fig = go.Figure(data=[surface, rope])
    fig.update_layout(
        title=title,
        width=760,
        height=720,
        margin={"l": 0, "r": 0, "t": 45, "b": 0},
        scene={
            "aspectmode": "data",
            "camera": {"eye": {"x": 2.5, "y": 2.5, "z": 1.0}},
            "xaxis": {"title": "X [Mm]", "showgrid": False, "zeroline": False, "showticklabels": False},
            "yaxis": {"title": "Y [Mm]", "showgrid": False, "zeroline": False, "showticklabels": False},
            "zaxis": {"title": "Z [Mm]", "showgrid": False, "zeroline": False, "showticklabels": False},
        },
    )
    fig.write_html(output, include_plotlyjs="cdn", full_html=True)


def spherical_surface_xyz(lon_degree, lat_degree, radius):
    phi = np.radians(lon_degree)
    theta = np.radians(90.0 - lat_degree)
    x = radius * np.sin(theta) * np.cos(phi)
    y = radius * np.sin(theta) * np.sin(phi)
    z = radius * np.cos(theta)
    return x, y, z


def plotly_colorscale_for_cmap(cmap):
    if cmap not in plt.colormaps():
        return "Viridis"
    mpl_cmap = plt.get_cmap(cmap)
    samples = np.linspace(0.0, 1.0, 16)
    colorscale = []
    for value in samples:
        red, green, blue, _ = mpl_cmap(value)
        colorscale.append([float(value), f"rgb({int(red * 255)},{int(green * 255)},{int(blue * 255)})"])
    return colorscale


def pick_path_points(
    data,
    extent,
    footpoints,
    output,
    *,
    circle_radius_degree,
    cmap,
    image_power,
    image_vmin,
    image_vmax,
):
    fig, ax = plt.subplots(figsize=(8, 6), dpi=160)
    image = np.abs(data) ** image_power
    ax.imshow(
        image,
        origin="lower",
        extent=extent,
        cmap=cmap if cmap in plt.colormaps() else "gray",
        aspect="auto",
        vmin=image_vmin,
        vmax=image_vmax,
    )
    for (lon, lat), color in zip(footpoints, ["black", "white"], strict=True):
        ax.add_patch(Circle((lon, lat), circle_radius_degree, facecolor="none", edgecolor=color, linewidth=1.5))
    ax.set_title("Click path points in order, press Enter when done")
    points = plt.ginput(n=-1, timeout=0)
    if not points:
        raise RuntimeError("no path points were clicked")
    points = np.asarray(points, dtype=float)
    ax.plot(points[:, 0], points[:, 1], "o-", color="yellow", markersize=3, linewidth=1)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    return points


def pick_path_points_on_submap(
    smap,
    output,
    *,
    circle_radius_arcsec,
    smooth_points,
    smoothing,
    reference_path_coords=None,
):
    fig = plt.figure(figsize=(8, 8), dpi=160)
    ax = fig.add_subplot(projection=smap)
    smap.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
    set_solar_xy_labels(ax)
    if reference_path_coords is not None:
        ax.plot_coord(
            SkyCoord(reference_path_coords),
            color="deepskyblue",
            linestyle="-.",
            linewidth=1.4,
            marker=None,
            alpha=0.9,
        )
    raw_points = []
    finished = {"done": False}
    raw_line, = ax.plot([], [], "o-", color="tab:orange", markersize=4, linewidth=1.0, transform=ax.get_transform("pixel"))
    ax.set_title("Click path points; delete/backspace removes last; enter finishes")

    def redraw_raw_points():
        if raw_points:
            pixels = np.asarray(raw_points, dtype=float)
            raw_line.set_data(pixels[:, 0], pixels[:, 1])
        else:
            raw_line.set_data([], [])
        fig.canvas.draw_idle()

    def on_click(event):
        if event.inaxes != ax or event.button != 1 or event.xdata is None or event.ydata is None:
            return
        raw_points.append([float(event.xdata), float(event.ydata)])
        redraw_raw_points()

    def on_key(event):
        key = (event.key or "").lower()
        if key in {"delete", "backspace"}:
            if raw_points:
                raw_points.pop()
                redraw_raw_points()
        elif key in {"enter", "return"}:
            finished["done"] = True
            plt.close(fig)

    fig.canvas.mpl_connect("button_press_event", on_click)
    fig.canvas.mpl_connect("key_press_event", on_key)
    plt.show()

    if len(raw_points) < 2:
        raise RuntimeError("at least two path points are required")

    clicked_pixels = np.asarray(raw_points, dtype=float)
    path_pixels = smooth_pixel_path(clicked_pixels, smooth_points=smooth_points, smoothing=smoothing)
    np.save(output.with_name("clicked_points_raw.npy"), clicked_pixels)
    coords_hpc = [smap.pixel_to_world(x * u.pix, y * u.pix) for x, y in path_pixels]
    points = np.asarray([hpc_to_stonyhurst(coord) for coord in coords_hpc], dtype=float)
    plot_path_on_submap(
        smap,
        coords_hpc,
        output,
        circle_radius_arcsec=circle_radius_arcsec,
        reference_path_coords=reference_path_coords,
    )
    plt.close(fig)
    return points, coords_hpc


def smooth_pixel_path(pixels, *, smooth_points, smoothing):
    pixels = np.asarray(pixels, dtype=float)
    if pixels.ndim != 2 or pixels.shape[1] != 2:
        raise ValueError("pixels must have shape (N, 2)")
    if pixels.shape[0] > 1:
        keep = np.concatenate([[True], np.linalg.norm(np.diff(pixels, axis=0), axis=1) > 1e-6])
        pixels = pixels[keep]
    if smooth_points <= 0 or pixels.shape[0] < 3:
        return pixels

    distances = np.concatenate([[0.0], np.cumsum(np.linalg.norm(np.diff(pixels, axis=0), axis=1))])
    if distances[-1] == 0.0:
        return pixels
    u_values = distances / distances[-1]
    k = min(3, pixels.shape[0] - 1)
    u_new = np.linspace(0.0, 1.0, smooth_points)
    try:
        weights = np.ones(pixels.shape[0], dtype=float)
        if smoothing > 0.0:
            weights[[0, -1]] = 1e3
        tck, _ = splprep([pixels[:, 0], pixels[:, 1]], u=u_values, w=weights, s=smoothing, k=k)
        x_new, y_new = splev(u_new, tck)
    except ValueError:
        x_new = np.interp(u_new, u_values, pixels[:, 0])
        y_new = np.interp(u_new, u_values, pixels[:, 1])
    smoothed = np.column_stack([x_new, y_new])
    smoothed[0] = pixels[0]
    smoothed[-1] = pixels[-1]
    return smoothed


def plot_path_on_submap(
    smap,
    path_coords,
    output,
    *,
    circle_radius_arcsec,
    reference_path_coords=None,
    boundary_coords=None,
):
    fig = plt.figure(figsize=(8, 8), dpi=160)
    ax = fig.add_subplot(projection=smap)
    smap.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
    set_solar_xy_labels(ax)
    if reference_path_coords is not None:
        ax.plot_coord(
            SkyCoord(reference_path_coords),
            color="deepskyblue",
            linestyle="-.",
            linewidth=1.4,
            marker=None,
            alpha=0.9,
        )
    ax.plot_coord(SkyCoord(path_coords), color="lime", linestyle="--", linewidth=1.6, marker=None)
    if boundary_coords is not None:
        ax.plot_coord(SkyCoord(boundary_coords), color="cyan", linestyle="-", linewidth=1.2, marker=None, alpha=0.95)
    endpoint_coords = [path_coords[0], path_coords[-1]]
    colors = ["white", "black"]
    for coord, color in zip(endpoint_coords, colors, strict=True):
        ax.plot_coord(
            coord,
            "o",
            markersize=13,
            markerfacecolor=color,
            markeredgecolor=color,
            markeredgewidth=2,
            alpha=0.45,
        )
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_projected_path_on_aia(smap, path_lonlatr, output, *, circle_radius_arcsec, boundary_lonlatr=None):
    coords = lonlatr_to_hpc(path_lonlatr, smap)
    boundary_coords = lonlatr_to_hpc(boundary_lonlatr, smap) if boundary_lonlatr is not None else None
    plot_path_on_submap(
        smap,
        coords,
        output,
        circle_radius_arcsec=circle_radius_arcsec,
        boundary_coords=boundary_coords,
    )


def ensure_footpoint_endpoints(points, footpoints):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")
    if points.shape[0] == 0:
        return footpoints.copy()
    start_distance = np.linalg.norm(points[0] - footpoints[0])
    end_distance = np.linalg.norm(points[-1] - footpoints[1])
    if start_distance > 1e-8:
        points = np.vstack([footpoints[0], points])
    if end_distance > 1e-8:
        points = np.vstack([points, footpoints[1]])
    return points


def hpc_to_stonyhurst(coord):
    hgs = coord.transform_to(frames.HeliographicStonyhurst)
    return [hgs.lon.to_value(u.degree), hgs.lat.to_value(u.degree)]


def stonyhurst_to_hpc(lon, lat, smap):
    coord = SkyCoord(
        lon * u.degree,
        lat * u.degree,
        frame=frames.HeliographicStonyhurst,
        obstime=smap.date,
        rsun=smap.rsun_meters,
        observer="earth",
    )
    return coord.transform_to(smap.coordinate_frame)


def lonlatr_to_hpc(lon_lat_r, smap):
    lon_lat_r = np.asarray(lon_lat_r, dtype=float)
    coords = SkyCoord(
        lon_lat_r[:, 0] * u.degree,
        lon_lat_r[:, 1] * u.degree,
        lon_lat_r[:, 2] * smap.rsun_meters,
        frame=frames.HeliographicStonyhurst,
        obstime=smap.date,
    )
    return coords.transform_to(smap.coordinate_frame)


def set_solar_xy_labels(ax):
    ax.set_xlabel("Solar-X")
    ax.set_ylabel("Solar-Y")


if __name__ == "__main__":
    raise SystemExit(main())

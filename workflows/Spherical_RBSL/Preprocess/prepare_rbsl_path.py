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
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.patches import Circle, Polygon
from scipy.interpolate import splev, splprep
from scipy.ndimage import map_coordinates
from sunpy.coordinates import frames
from sunpy.map import Map

try:
    from .rbsl_core.maps import get_submap
    from .rbsl_core.spherical import (
        footpoint_fluxes_from_hmi_map,
        flux_rope_axis_from_points,
        flux_rope_boundary_fieldline,
        lon_lat_centers,
        sample_sunpy_map_on_stonyhurst,
        write_path_dat,
    )
except ImportError:  # Direct execution: python Preprocess/prepare_rbsl_path.py
    from rbsl_core.maps import get_submap
    from rbsl_core.spherical import (
        footpoint_fluxes_from_hmi_map,
        flux_rope_axis_from_points,
        flux_rope_boundary_fieldline,
        lon_lat_centers,
        sample_sunpy_map_on_stonyhurst,
        write_path_dat,
    )

LIGHT_TUBE_COLOR = "#8EEBFF"


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    image_clip_percent = resolve_image_clip_percent(args)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    image_output_dir = resolve_image_output_dir(output_dir, rope_radius_mm=args.rope_radius_mm)
    image_output_dir.mkdir(parents=True, exist_ok=True)

    pick_map = resolve_pick_map(args)
    submap = pick_map
    reference_path_coords = None
    if args.reference_points:
        reference_points = np.load(args.reference_points)
        reference_path_coords = [stonyhurst_to_hpc(lon, lat, pick_map) for lon, lat in reference_points]

    if args.points:
        points = np.load(args.points)
    elif args.path_points_stonyhurst:
        points = np.asarray(args.path_points_stonyhurst, dtype=float).reshape(-1, 2)
    else:
        points, _ = pick_path_points_on_submap(
            submap,
            image_output_dir / "picked_points.png",
            circle_radius_arcsec=args.footpoint_circle_arcsec,
            smooth_points=args.smooth_path_points,
            smoothing=args.smooth_path_smoothing,
            reference_path_coords=reference_path_coords,
        )

    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or points.shape[0] < 2:
        raise ValueError("path points must have shape (N, 2) with at least two points")
    footpoints = points[[0, -1]]
    footpoints_hpc = [stonyhurst_to_hpc(lon, lat, pick_map) for lon, lat in footpoints]
    (output_dir / "selected_footpoints.png").unlink(missing_ok=True)
    solar_radius_mm = args.solar_radius_mm
    if solar_radius_mm is None:
        solar_radius_mm = pick_map.rsun_meters.to_value(u.Mm)

    np.save(output_dir / "pts.npy", points)
    if args.points or args.path_points_stonyhurst:
        plot_path_on_submap(
            submap,
            [stonyhurst_to_hpc(lon, lat, pick_map) for lon, lat in points],
            image_output_dir / "picked_points.png",
            circle_radius_arcsec=args.footpoint_circle_arcsec,
            reference_path_coords=reference_path_coords,
        )

    lon_range = tuple(args.lon_range) if args.lon_range is not None else None
    lat_range = tuple(args.lat_range) if args.lat_range is not None else None
    lonlat_range_inferred = False
    if args.write_3d and lon_range is None and lat_range is None:
        lon_range, lat_range = infer_lonlat_range(points)
        lonlat_range_inferred = True
    needs_lonlat_grid = lon_range is not None or lat_range is not None or args.write_3d
    if needs_lonlat_grid and (lon_range is None or lat_range is None):
        raise ValueError("--lon-range and --lat-range must be provided together; 3D-only runs can omit both")
    lon = lat = base_data = extent = None
    if needs_lonlat_grid:
        lon, lat = lon_lat_centers(lon_range, lat_range, args.resolution_degree, nghost=args.nghost)
        base_map = resolve_base_map(args, pick_map)
        base_data = sample_map_on_stonyhurst(base_map, lon, lat)
        extent = [lon[0], lon[-1], lat[0], lat[-1]]
        (output_dir / "base_footpoints.png").unlink(missing_ok=True)
        (output_dir / "picked_points_lonlat.png").unlink(missing_ok=True)

    hmi_map = resolve_hmi_map(args)
    hmi_data = None
    hmi_map_type = None
    if hmi_map is not None:
        hmi_map_type = classify_hmi_map(hmi_map)
    if hmi_map is not None and needs_lonlat_grid:
        hmi_data, _ = sample_sunpy_map_on_stonyhurst(hmi_map, lon, lat, rotate=False)
        plot_hmi_lonlat_base(
            hmi_data,
            extent,
            footpoints,
            image_output_dir / "hmi_footpoints.png",
            circle_radius_degree=args.footpoint_circle_degree,
            path_points=None,
            vmin=args.hmi_vmin,
            vmax=args.hmi_vmax,
        )
        plot_hmi_lonlat_base(
            hmi_data,
            extent,
            footpoints,
            image_output_dir / "hmi_picked_points.png",
            circle_radius_degree=args.footpoint_circle_degree,
            path_points=points,
            vmin=args.hmi_vmin,
            vmax=args.hmi_vmax,
        )
    else:
        (image_output_dir / "hmi_footpoints.png").unlink(missing_ok=True)
        (image_output_dir / "hmi_picked_points.png").unlink(missing_ok=True)

    footpoint_fluxes = None
    footpoint_flux_note = None
    footpoint_flux_radius_mm = args.footpoint_flux_radius_mm if args.footpoint_flux_radius_mm is not None else 20.0
    if hmi_map is not None:
        footpoint_fluxes = footpoint_fluxes_from_hmi_map(
            footpoints,
            footpoint_flux_radius_mm,
            args.hmi,
            solar_radius_mm=solar_radius_mm,
            resolution_degree=args.footpoint_flux_resolution_degree,
            rotate=True,
        )
        if hmi_map_type != "br":
            footpoint_flux_note = (
                "Fluxes were computed from the provided HMI map values after rotation/sampling. "
                "Because the input is not an HMI Br map (for example HMI45 / LOS magnetogram), "
                "these are LOS-based approximate fluxes rather than strictly radial Br fluxes."
            )
    footpoint_value_labels = None
    if footpoint_fluxes is not None:
        footpoint_value_labels = [format_flux_label(item["signed_flux_Mx"]) for item in footpoint_fluxes["results"]]

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
        image_output_dir / "projected_path_aia.png",
        circle_radius_arcsec=args.footpoint_circle_arcsec,
        boundary_lonlatr=boundary_path,
        tube_radius_mm=args.rope_radius_mm if args.rope_radius_mm and args.rope_radius_mm > 0.0 else None,
        tube_alpha=0.25,
        solar_radius_mm=solar_radius_mm,
        show_boundary_line=not args.hide_boundary_line,
        image_vmin=args.hmi_vmin if is_hmi_map(submap) else None,
        image_vmax=args.hmi_vmax if is_hmi_map(submap) else None,
        footpoint_value_labels=footpoint_value_labels,
        show_length_label=not args.hide_length_label,
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
        if hmi_data is not None:
            rope_3d_hmi_html = output_dir / "rope_3d_hmi.html"
            plot_3d_rope_surface(
                result["path"],
                lon,
                lat,
                hmi_data,
                rope_3d_hmi_html,
                solar_radius_mm=solar_radius_mm,
                surface_range=(args.hmi_vmin, args.hmi_vmax),
                surface_colorscale="Gray",
                title=f"3D flux rope on HMI map ({hmi_map_type or 'unknown'})",
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
        "hmi_background_used": hmi_map is not None,
        "hmi_map_type": hmi_map_type,
        "rope_3d_html": str(rope_3d_html) if rope_3d_html is not None else None,
        "rope_3d_hmi_html": str(rope_3d_hmi_html) if rope_3d_hmi_html is not None else None,
        "path_dat": str(path),
        "lon_range": list(lon_range) if lon_range is not None else None,
        "lat_range": list(lat_range) if lat_range is not None else None,
        "lonlat_range_inferred": lonlat_range_inferred,
        "image_scale": args.image_scale,
        "image_clip_percent": list(image_clip_percent),
        "footpoint_fluxes": footpoint_fluxes,
        "footpoint_flux_note": footpoint_flux_note,
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
    if footpoint_fluxes is not None:
        if footpoint_flux_note:
            print(f"note: {footpoint_flux_note}")
        for item in footpoint_fluxes["results"]:
            print(
                f"footpoint {item['footpoint_index']} flux (r={item['radius_mm']:.3f} Mm): "
                f"signed={item['signed_flux_Mx']:.6e} Mx, unsigned={item['unsigned_flux_Mx']:.6e} Mx"
            )
    print(f"wrote: {path}")
    return 0


def format_flux_label(value_mx: float) -> str:
    value_mx = float(value_mx)
    if value_mx == 0.0:
        return r"$0.0 \times 10^{0}\ \mathrm{Mx}$"
    exponent = int(np.floor(np.log10(abs(value_mx))))
    mantissa = value_mx / (10.0**exponent)
    return rf"${mantissa:.1f} \times 10^{{{exponent}}}\ \mathrm{{Mx}}$"


def resolve_image_output_dir(output_dir: Path, *, rope_radius_mm: float | None):
    rope_radius_mm = 0.0 if rope_radius_mm is None else float(rope_radius_mm)
    if rope_radius_mm > 0.0 and "tube" not in output_dir.parts:
        return output_dir / "tube"
    return output_dir


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--aia",
        help="AIA FITS used for path picking and projected overlay.",
    )
    parser.add_argument("--hmi", help="Optional HMI FITS used for HMI overlays / fluxes, or as the picking background if --aia is omitted.")
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
    parser.add_argument("--hmi-vmin", type=float, default=-300.0)
    parser.add_argument("--hmi-vmax", type=float, default=300.0)
    parser.add_argument(
        "--footpoint-flux-radius-mm",
        type=float,
        help="Footpoint flux radius in Mm. If --hmi is provided and this is omitted, defaults to 20 Mm.",
    )
    parser.add_argument("--footpoint-flux-resolution-degree", type=float, default=0.02)
    parser.add_argument("--write-3d", action="store_true", help="Write rope_3d HTML files; disabled by default to save disk space.")
    parser.add_argument("--rope-3d-width", type=float, default=10.0)
    parser.add_argument("--rope-radius-mm", type=float, default=0.0, help="Optional radius for a helical boundary field line in Mm.")
    parser.add_argument("--rope-twist", type=float, default=0.0, help="Number of turns for the optional boundary field line.")
    parser.add_argument("--hide-boundary-line", action="store_true", help="Hide the 2D projected twist/boundary line while keeping tube/boundary files.")
    parser.add_argument("--hide-length-label", action="store_true", help="Hide the top-right projected rope-length label.")
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


def resolve_pick_map(args):
    extent = tuple(args.submap_arcsec)
    pick_path = args.aia or args.hmi
    if pick_path is None:
        raise ValueError("one of --aia or --hmi must be provided")
    return get_submap(prepare_map_for_plotting(Map(pick_path)), extent)


def resolve_base_map(args, default_map):
    base_path = args.base_map or args.base_aia
    if base_path is None:
        return default_map
    return prepare_map_for_plotting(Map(base_path))


def resolve_hmi_map(args):
    if args.hmi is None:
        return None
    return prepare_map_for_plotting(Map(args.hmi))


def is_hmi_map(smap):
    if bool(smap.meta.get("SOLBRMAP")):
        return True
    detector = str(getattr(smap, "detector", "") or "").upper()
    instrument = str(getattr(smap, "instrument", "") or "").upper()
    telescop = str(smap.meta.get("telescop", "") or "").upper()
    content = str(smap.meta.get("content", "") or "").upper()
    return detector == "HMI" or "HMI" in instrument or "HMI" in telescop or "MAGNETOGRAM" in content


def prepare_map_for_plotting(smap):
    return smap.rotate() if is_hmi_map(smap) else smap


def classify_hmi_map(smap):
    content = str(smap.meta.get("content", "") or "").upper()
    if bool(smap.meta.get("SOLBRMAP")) or content == "BR" or " RADIAL" in content or content.startswith("BR"):
        return "br"
    if "MAGNETOGRAM" in content:
        return "los"
    return "generic_hmi"


def plot_submap_image(smap, ax, *, vmin=None, vmax=None):
    artist = None
    if vmin is not None or vmax is not None:
        plot_kwargs = {"axes": ax}
        if vmin is not None:
            plot_kwargs["vmin"] = vmin
        if vmax is not None:
            plot_kwargs["vmax"] = vmax
        try:
            artist = smap.plot(**plot_kwargs)
            if hasattr(artist, "set_clim"):
                artist.set_clim(vmin=vmin, vmax=vmax)
            return artist
        except ValueError as exc:
            if "Cannot manually specify vmax" not in str(exc):
                raise
    try:
        artist = smap.plot(axes=ax, clip_interval=(1, 99.99) * u.percent)
    except ValueError as exc:
        if "Cannot manually specify vmax" not in str(exc):
            raise
        artist = smap.plot(axes=ax)
    if (vmin is not None or vmax is not None) and hasattr(artist, "set_clim"):
        artist.set_clim(vmin=vmin, vmax=vmax)
    return artist


def format_projection_title(smap):
    if is_hmi_map(smap):
        hmi_type = classify_hmi_map(smap)
        if hmi_type == "br":
            prefix = r"HMI $B_r$"
        elif hmi_type == "los":
            prefix = "HMI LOS"
        else:
            prefix = "HMI"
        return f"{prefix} {smap.date.isot[:19].replace('T', ' ')}"
    return None


def format_colorbar_label(smap):
    if is_hmi_map(smap):
        return r"$B_r$ [G]" if classify_hmi_map(smap) == "br" else "Magnetic field [G]"
    return None


def pick_two_hpc_points(smap, output):
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(projection=smap)
    plot_submap_image(smap, ax)
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
    plot_submap_image(smap, ax)
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
    plot_submap_image(smap, ax)
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
    output.with_name("clicked_points_raw.npy").unlink(missing_ok=True)
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
    footpoint_radius_arcsec=None,
    reference_path_coords=None,
    boundary_coords=None,
    tube_vertices_pixel=None,
    tube_color=LIGHT_TUBE_COLOR,
    tube_alpha=0.25,
    annotation=None,
    corner_label=None,
    image_vmin=None,
    image_vmax=None,
    footpoint_value_labels=None,
    title_override=None,
    colorbar_label=None,
    show_colorbar=False,
):
    fig = plt.figure(figsize=(8, 8), dpi=160)
    ax = fig.add_subplot(projection=smap)
    image_artist = plot_submap_image(smap, ax, vmin=image_vmin, vmax=image_vmax)
    ax.grid(alpha=0)
    set_solar_xy_labels(ax)
    if title_override is not None:
        ax.set_title(title_override)
    if show_colorbar and image_artist is not None:
        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4.5%", pad=0.08, axes_class=plt.Axes)
        cbar = fig.colorbar(image_artist, cax=cax)
        if colorbar_label is not None:
            cbar.set_label(colorbar_label)
    if reference_path_coords is not None:
        ax.plot_coord(
            SkyCoord(reference_path_coords),
            color="deepskyblue",
            linestyle="-.",
            linewidth=1.4,
            marker=None,
            alpha=0.9,
        )
    if tube_vertices_pixel is not None and len(tube_vertices_pixel) >= 3:
        ax.add_patch(
            Polygon(
                tube_vertices_pixel,
                closed=True,
                facecolor=tube_color,
                edgecolor="none",
                alpha=tube_alpha,
                transform=ax.get_transform("pixel"),
                zorder=2,
            )
        )
    ax.plot_coord(SkyCoord(path_coords), color="lime", linestyle="--", linewidth=1.6, marker=None)
    if boundary_coords is not None:
        ax.plot_coord(SkyCoord(boundary_coords), color=tube_color, linestyle="-", linewidth=1.2, marker=None, alpha=0.9)
    endpoint_coords = [path_coords[0], path_coords[-1]]
    colors = ["white", "black"]
    symbols = ["+", "-"]
    symbol_colors = ["black", "white"]
    radius_arcsec = circle_radius_arcsec if footpoint_radius_arcsec is None else footpoint_radius_arcsec
    for index, (coord, color, symbol, symbol_color) in enumerate(
        zip(endpoint_coords, colors, symbols, symbol_colors, strict=True)
    ):
        circle_vertices_pixel = hpc_circle_arcsec_to_pixel_vertices(
            smap,
            coord.Tx.to_value(u.arcsec),
            coord.Ty.to_value(u.arcsec),
            radius_arcsec,
        )
        ax.add_patch(
            Polygon(
                circle_vertices_pixel,
                closed=True,
                facecolor=color,
                edgecolor=color,
                linewidth=1.2,
                alpha=0.35,
                transform=ax.get_transform("pixel"),
                zorder=4,
            )
        )
        center_pixel = hpc_arcsec_polygon_to_pixel_vertices(
            smap,
            np.asarray([[coord.Tx.to_value(u.arcsec), coord.Ty.to_value(u.arcsec)]], dtype=float),
        )[0]
        ax.text(
            center_pixel[0],
            center_pixel[1],
            symbol,
            color=symbol_color,
            fontsize=14,
            fontweight="bold",
            ha="center",
            va="center",
            transform=ax.get_transform("pixel"),
            zorder=5,
        )
        if footpoint_value_labels is not None and index < len(footpoint_value_labels):
            pixel_radius = float(np.max(np.linalg.norm(circle_vertices_pixel - center_pixel, axis=1)))
            ax.text(
                center_pixel[0],
                center_pixel[1] - 1.9 * pixel_radius,
                footpoint_value_labels[index],
                color="white",
                fontsize=12,
                ha="center",
                va="top",
                transform=ax.get_transform("pixel"),
                zorder=5,
            )
    if annotation is not None:
        ax.annotate(
            annotation["text"],
            xy=annotation["xy"],
            xytext=annotation["xytext"],
            xycoords=ax.get_transform("pixel"),
            textcoords=ax.get_transform("pixel"),
            color=annotation.get("color", LIGHT_TUBE_COLOR),
            fontsize=annotation.get("fontsize", 12),
            ha=annotation.get("ha", "center"),
            va=annotation.get("va", "center"),
            arrowprops={
                "arrowstyle": "->",
                "color": annotation.get("color", LIGHT_TUBE_COLOR),
                "lw": annotation.get("lw", 1.6),
                "alpha": annotation.get("alpha", 0.95),
            },
            zorder=5,
        )
    if corner_label is not None:
        ax.text(
            0.97,
            0.97,
            corner_label["text"],
            transform=ax.transAxes,
            ha=corner_label.get("ha", "right"),
            va=corner_label.get("va", "top"),
            fontsize=corner_label.get("fontsize", 13),
            color=corner_label.get("color", "white"),
            zorder=6,
        )
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_projected_path_on_aia(
    smap,
    path_lonlatr,
    output,
    *,
    circle_radius_arcsec,
    boundary_lonlatr=None,
    tube_radius_mm=None,
    tube_alpha=0.25,
    solar_radius_mm,
    show_boundary_line=True,
    image_vmin=None,
    image_vmax=None,
    footpoint_value_labels=None,
    show_length_label=True,
):
    coords = lonlatr_to_hpc(path_lonlatr, smap)
    boundary_coords = lonlatr_to_hpc(boundary_lonlatr, smap) if (boundary_lonlatr is not None and show_boundary_line) else None
    tube_vertices_pixel = None
    footpoint_radius_arcsec = circle_radius_arcsec
    annotation = None
    corner_label = None
    if tube_radius_mm is not None and tube_radius_mm > 0.0:
        projected_xy = hpc_coords_to_xy_arcsec(coords)
        half_width_arcsec = tube_radius_mm / solar_radius_mm * smap.rsun_obs.to_value(u.arcsec)
        tube_polygon_arcsec = build_constant_width_tube_polygon(projected_xy, half_width_arcsec)
        tube_vertices_pixel = hpc_arcsec_polygon_to_pixel_vertices(smap, tube_polygon_arcsec)
        footpoint_radius_arcsec = half_width_arcsec
        annotation = build_tube_radius_annotation(smap, projected_xy, half_width_arcsec, tube_radius_mm)
    if show_length_label:
        projected_length_mm = estimate_path_length_mm(path_lonlatr, solar_radius_mm)
        corner_label = build_length_corner_label(projected_length_mm)
    plot_path_on_submap(
        smap,
        coords,
        output,
        circle_radius_arcsec=circle_radius_arcsec,
        footpoint_radius_arcsec=footpoint_radius_arcsec,
        boundary_coords=boundary_coords,
        tube_vertices_pixel=tube_vertices_pixel,
        tube_alpha=tube_alpha,
        annotation=annotation,
        corner_label=corner_label,
        image_vmin=image_vmin,
        image_vmax=image_vmax,
        footpoint_value_labels=footpoint_value_labels,
        title_override=format_projection_title(smap),
        colorbar_label=format_colorbar_label(smap),
        show_colorbar=is_hmi_map(smap),
    )


def hpc_coords_to_xy_arcsec(coords):
    return np.column_stack([coords.Tx.to_value(u.arcsec), coords.Ty.to_value(u.arcsec)])


def build_constant_width_tube_polygon(xy_arcsec, half_width_arcsec, *, cap_points=48):
    xy_arcsec = np.asarray(xy_arcsec, dtype=float)
    if xy_arcsec.ndim != 2 or xy_arcsec.shape[1] != 2 or xy_arcsec.shape[0] < 2:
        raise ValueError("xy_arcsec must have shape (N, 2) with at least two points")

    tangent = np.empty_like(xy_arcsec)
    tangent[0] = xy_arcsec[1] - xy_arcsec[0]
    tangent[-1] = xy_arcsec[-1] - xy_arcsec[-2]
    tangent[1:-1] = xy_arcsec[2:] - xy_arcsec[:-2]
    norm = np.linalg.norm(tangent, axis=1, keepdims=True)
    norm[norm == 0.0] = 1.0
    tangent /= norm
    normal = np.column_stack([-tangent[:, 1], tangent[:, 0]])
    left = xy_arcsec + half_width_arcsec * normal
    right = xy_arcsec - half_width_arcsec * normal
    end_cap = semicircle_cap_points(
        xy_arcsec[-1],
        tangent[-1],
        normal[-1],
        half_width_arcsec,
        theta_start=0.0,
        theta_end=np.pi,
        npoints=cap_points,
    )
    start_cap = semicircle_cap_points(
        xy_arcsec[0],
        tangent[0],
        normal[0],
        half_width_arcsec,
        theta_start=np.pi,
        theta_end=2.0 * np.pi,
        npoints=cap_points,
    )
    return np.vstack([left, end_cap[1:], right[-2::-1], start_cap[1:]])


def semicircle_cap_points(center, tangent, normal, radius, *, theta_start, theta_end, npoints):
    center = np.asarray(center, dtype=float)
    tangent = np.asarray(tangent, dtype=float)
    normal = np.asarray(normal, dtype=float)
    angles = np.linspace(theta_start, theta_end, npoints)
    basis_t = tangent[None, :]
    basis_n = normal[None, :]
    return center + radius * (np.cos(angles)[:, None] * basis_n + np.sin(angles)[:, None] * basis_t)


def hpc_arcsec_polygon_to_pixel_vertices(smap, polygon_arcsec):
    polygon_arcsec = np.asarray(polygon_arcsec, dtype=float)
    coords = SkyCoord(
        polygon_arcsec[:, 0] * u.arcsec,
        polygon_arcsec[:, 1] * u.arcsec,
        frame=smap.coordinate_frame,
    )
    pixel_x, pixel_y = smap.world_to_pixel(coords)
    return np.column_stack([pixel_x.to_value(u.pix), pixel_y.to_value(u.pix)])


def hpc_circle_arcsec_to_pixel_vertices(smap, x_arcsec, y_arcsec, radius_arcsec, *, npoints=96):
    angles = np.linspace(0.0, 2.0 * np.pi, npoints, endpoint=False)
    polygon_arcsec = np.column_stack(
        [
            x_arcsec + radius_arcsec * np.cos(angles),
            y_arcsec + radius_arcsec * np.sin(angles),
        ]
    )
    return hpc_arcsec_polygon_to_pixel_vertices(smap, polygon_arcsec)


def build_tube_radius_annotation(smap, projected_xy_arcsec, half_width_arcsec, tube_radius_mm):
    projected_xy_arcsec = np.asarray(projected_xy_arcsec, dtype=float)
    mid = projected_xy_arcsec.shape[0] // 2
    if projected_xy_arcsec.shape[0] < 3:
        return None
    prev_idx = max(mid - 1, 0)
    next_idx = min(mid + 1, projected_xy_arcsec.shape[0] - 1)
    tangent = projected_xy_arcsec[next_idx] - projected_xy_arcsec[prev_idx]
    tangent_norm = np.linalg.norm(tangent)
    if tangent_norm == 0.0:
        tangent = np.array([1.0, 0.0])
    else:
        tangent = tangent / tangent_norm
    normal = np.array([-tangent[1], tangent[0]])
    xy_arcsec = projected_xy_arcsec[mid]
    xytext_arcsec = xy_arcsec + normal * (2.2 * half_width_arcsec)
    xy_pix = hpc_arcsec_polygon_to_pixel_vertices(smap, np.asarray([xy_arcsec]))[0]
    xytext_pix = hpc_arcsec_polygon_to_pixel_vertices(smap, np.asarray([xytext_arcsec]))[0]
    return {
        "text": rf"$a_0 = {tube_radius_mm:.0f}\,\mathrm{{Mm}}$",
        "xy": xy_pix,
        "xytext": xytext_pix,
        "color": LIGHT_TUBE_COLOR,
        "fontsize": 12,
        "ha": "center",
        "va": "center",
        "lw": 1.6,
        "alpha": 0.95,
    }


def build_length_corner_label(projected_length_mm):
    return {
        "text": rf"$\mathrm{{Length\ of\ flux\ rope}}\ \ L = {projected_length_mm:.1f}\,\mathrm{{Mm}}$",
        "color": "white",
        "fontsize": 13,
        "ha": "right",
        "va": "top",
    }


def estimate_path_length_mm(path_lonlatr, solar_radius_mm):
    arr = np.asarray(path_lonlatr, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 3 or arr.shape[0] < 2:
        raise ValueError("path_lonlatr must have shape (N, 3) with at least two rows")
    lon = np.deg2rad(arr[:, 0])
    lat = np.deg2rad(arr[:, 1])
    radius = arr[:, 2] * solar_radius_mm
    xyz = np.column_stack(
        [
            radius * np.cos(lat) * np.cos(lon),
            radius * np.cos(lat) * np.sin(lon),
            radius * np.sin(lat),
        ]
    )
    return float(np.linalg.norm(np.diff(xyz, axis=0), axis=1).sum())


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

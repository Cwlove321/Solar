#!/usr/bin/env python
"""Overlay an RBSL path on Solar Orbiter FSI and STEREO/EUVI images."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib.patches import Polygon
from sunpy.map import Map

try:
    from .rbsl_core.maps import get_submap
    from .rbsl_core.spherical import read_path_dat
except ImportError:  # Direct execution: python Preprocess/plot_multiview_path.py
    from rbsl_core.maps import get_submap
    from rbsl_core.spherical import read_path_dat

LIGHT_TUBE_COLOR = "#8EEBFF"
EUVI_TUBE_COLOR = "#F6B26B"


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    path = read_path_dat(args.path_dat, subtract_lon_180=True)
    footpoints = load_footpoints(args, path)
    tube_radius_mm, tube_alpha = load_tube_style(args)
    output_dir = resolve_image_output_dir(Path(args.output_dir), tube_radius_mm=tube_radius_mm)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.fsi:
        fsi_output = output_dir / infer_fsi_output_name(args.fsi)
        plot_view(
            args.fsi,
            path,
            footpoints,
            fsi_output,
            submap_arcsec=args.fsi_submap_arcsec,
            clip_percent=tuple(args.fsi_clip_percent),
            title=args.fsi_title,
            path_color="lime",
            tube_radius_mm=tube_radius_mm,
            tube_alpha=tube_alpha,
            image_cmap=None,
            tube_color=LIGHT_TUBE_COLOR,
        )

    if args.euvi:
        plot_view(
            args.euvi,
            path,
            footpoints,
            output_dir / "stereo_euvi_path.png",
            submap_arcsec=args.euvi_submap_arcsec,
            clip_percent=tuple(args.euvi_clip_percent),
            title=args.euvi_title,
            path_color="magenta",
            tube_radius_mm=tube_radius_mm,
            tube_alpha=tube_alpha,
            image_cmap=None,
            tube_color=EUVI_TUBE_COLOR,
        )

    print(f"wrote: {output_dir}")
    return 0


def infer_fsi_output_name(map_path: str) -> str:
    smap = Map(map_path)
    wavelength = getattr(smap, "wavelength", None)
    if wavelength is not None:
        try:
            wave = int(round(wavelength.to_value(u.AA)))
            return f"solar_orbiter_fsi{wave}_path.png"
        except Exception:
            pass
    return "solar_orbiter_fsi_path.png"


def resolve_image_output_dir(output_dir: Path, *, tube_radius_mm: float | None):
    tube_radius_mm = None if tube_radius_mm is None else float(tube_radius_mm)
    if tube_radius_mm is not None and tube_radius_mm > 0.0 and "tube" not in output_dir.parts:
        return output_dir / "tube"
    return output_dir


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path-dat", required=True, help="Binary path.dat with lon+180, lat, r/Rsun.")
    parser.add_argument("--points", help="Optional pts.npy; first and last rows are used as footpoints.")
    parser.add_argument("--summary", help="Optional path_summary.json containing footpoints_stonyhurst.")
    parser.add_argument("--footpoints-stonyhurst", nargs=4, type=float, metavar=("LON0", "LAT0", "LON1", "LAT1"))
    parser.add_argument("--fsi", help="Solar Orbiter EUI/FSI FITS.")
    parser.add_argument("--euvi", help="STEREO EUVI FITS/FTS.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--fsi-submap-arcsec", nargs=4, type=float, help="Optional event-specific FSI crop; omitted uses the full image.")
    parser.add_argument("--euvi-submap-arcsec", nargs=4, type=float, help="Optional event-specific EUVI crop; omitted uses the full image.")
    parser.add_argument("--fsi-clip-percent", nargs=2, type=float, default=[10.0, 99.9])
    parser.add_argument("--euvi-clip-percent", nargs=2, type=float, default=[1.0, 99.9])
    parser.add_argument("--path-step", type=int, default=500)
    parser.add_argument("--fsi-title", default="Solar Orbiter FSI")
    parser.add_argument("--euvi-title", default="STEREO/EUVI")
    parser.add_argument("--tube-radius-mm", type=float, help="Optional projected tube half-width in Mm.")
    parser.add_argument("--tube-alpha", type=float, default=0.25, help="Tube fill alpha.")
    return parser.parse_args(argv)


def load_footpoints(args, path):
    if args.footpoints_stonyhurst:
        return np.asarray(args.footpoints_stonyhurst, dtype=float).reshape(2, 2)
    if args.summary:
        summary = json.loads(Path(args.summary).read_text(encoding="utf-8"))
        if "footpoints_stonyhurst" in summary:
            return np.asarray(summary["footpoints_stonyhurst"], dtype=float)
    if args.points:
        points = np.load(args.points)
        return np.column_stack([points[[0, -1]], np.ones(2)])
    return path[[0, -1], :2]


def load_tube_style(args):
    if args.tube_radius_mm is not None:
        return float(args.tube_radius_mm), float(args.tube_alpha)
    if args.summary:
        summary = json.loads(Path(args.summary).read_text(encoding="utf-8"))
        rope_radius_mm = summary.get("rope_radius_mm")
        if rope_radius_mm is not None and float(rope_radius_mm) > 0.0:
            return float(rope_radius_mm), float(args.tube_alpha)
    return None, float(args.tube_alpha)


def plot_view(
    map_path,
    path,
    footpoints,
    output,
    *,
    submap_arcsec,
    clip_percent,
    title,
    path_color,
    tube_radius_mm,
    tube_alpha,
    image_cmap,
    tube_color,
):
    raw_map = Map(map_path)
    smap = Map(np.asarray(raw_map.data, dtype=float), raw_map.meta).rotate(missing=0)
    submap = get_submap(smap, tuple(submap_arcsec)) if submap_arcsec is not None else smap
    observer = submap.coordinate_frame.observer
    obs_lon = observer.lon.to_value(u.degree)
    obs_lat = observer.lat.to_value(u.degree)

    projected_path = project_lonlatr_to_observer_hpc(path, obs_lon, obs_lat, submap.rsun_obs.to_value(u.arcsec))
    footpoint_path = ensure_radius_column(footpoints)
    projected_footpoints = project_lonlatr_to_observer_hpc(
        footpoint_path,
        obs_lon,
        obs_lat,
        submap.rsun_obs.to_value(u.arcsec),
    )

    rope_coords = SkyCoord(
        projected_path[:, 0] * u.arcsec,
        projected_path[:, 1] * u.arcsec,
        frame=submap.coordinate_frame,
    )
    footpoint_coords = SkyCoord(
        projected_footpoints[:, 0] * u.arcsec,
        projected_footpoints[:, 1] * u.arcsec,
        frame=submap.coordinate_frame,
    )

    fig = plt.figure(figsize=(6, 6), dpi=180)
    ax = fig.add_subplot(111, projection=submap)
    plot_kwargs = {"axes": ax, "clip_interval": clip_percent * u.percent}
    if image_cmap is not None:
        plot_kwargs["cmap"] = image_cmap
    submap.plot(**plot_kwargs)
    ax.grid(alpha=0)
    set_solar_xy_labels(ax)
    footpoint_radius_arcsec = None
    if tube_radius_mm is not None and tube_radius_mm > 0.0:
        footpoint_radius_arcsec = tube_radius_mm / 696.0 * submap.rsun_obs.to_value(u.arcsec)
        tube_polygon_arcsec = build_constant_width_tube_polygon(projected_path[:, :2], footpoint_radius_arcsec)
        tube_vertices_pixel = hpc_arcsec_polygon_to_pixel_vertices(submap, tube_polygon_arcsec)
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
    ax.plot_coord(rope_coords, color=path_color, linestyle="--", linewidth=1.6, marker=None)
    for coord, color, symbol, symbol_color in zip(
        footpoint_coords,
        ["white", "black"],
        ["+", "-"],
        ["black", "white"],
        strict=True,
    ):
        radius_arcsec = 10.0 if footpoint_radius_arcsec is None else footpoint_radius_arcsec
        circle_vertices_pixel = hpc_circle_arcsec_to_pixel_vertices(
            submap,
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
                linewidth=1.0,
                alpha=0.35,
                transform=ax.get_transform("pixel"),
                zorder=4,
            )
        )
        center_pixel = hpc_arcsec_polygon_to_pixel_vertices(
            submap,
            np.asarray([[coord.Tx.to_value(u.arcsec), coord.Ty.to_value(u.arcsec)]], dtype=float),
        )[0]
        ax.text(
            center_pixel[0],
            center_pixel[1],
            symbol,
            color=symbol_color,
            fontsize=13,
            fontweight="bold",
            ha="center",
            va="center",
            transform=ax.get_transform("pixel"),
            zorder=5,
        )
    wave = f"{submap.wavelength.to_value(u.AA):.0f} A" if submap.wavelength is not None else ""
    ax.set_title(f"{title} {wave} {submap.date.isot[:19]}")
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def project_lonlatr_to_observer_hpc(lon_lat_r, obs_lon, obs_lat, rsun_obs_arcsec):
    rotated = sun_rot_batch(np.asarray(lon_lat_r, dtype=float), obs_lon, obs_lat)
    return rotated[:, :2] * rsun_obs_arcsec


def ensure_radius_column(points):
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] not in {2, 3}:
        raise ValueError("points must have shape (N, 2) or (N, 3)")
    if points.shape[1] == 3:
        return points
    return np.column_stack([points, np.ones(points.shape[0])])


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


def sun_rot_batch(points, obs_lon, obs_lat):
    phi = np.radians(points[:, 0])
    theta = np.radians(90.0 - points[:, 1])
    radius = points[:, 2]
    xyz = np.column_stack(
        [
            radius * np.sin(theta) * np.cos(phi),
            radius * np.sin(theta) * np.sin(phi),
            radius * np.cos(theta),
        ]
    )
    rotation = rz(np.radians(obs_lon)) @ ry(np.radians(90.0 - obs_lat)) @ rz(np.radians(90.0))
    return xyz @ rotation


def ry(angle):
    ct = np.cos(angle)
    st = np.sin(angle)
    return np.array([[ct, 0.0, st], [0.0, 1.0, 0.0], [-st, 0.0, ct]])


def rz(angle):
    ct = np.cos(angle)
    st = np.sin(angle)
    return np.array([[ct, -st, 0.0], [st, ct, 0.0], [0.0, 0.0, 1.0]])


def set_solar_xy_labels(ax):
    ax.set_xlabel("Solar-X")
    ax.set_ylabel("Solar-Y")


if __name__ == "__main__":
    raise SystemExit(main())

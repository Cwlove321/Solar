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
from sunpy.map import Map

from solarutils.obs.maps import get_submap
from solarutils.sim.spherical import read_path_dat


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    path = read_path_dat(args.path_dat, subtract_lon_180=True)
    footpoints = path[[0, -1]]

    if args.fsi:
        plot_view(
            args.fsi,
            path,
            footpoints,
            output_dir / "solar_orbiter_fsi_path.png",
            submap_arcsec=args.fsi_submap_arcsec,
            clip_percent=tuple(args.fsi_clip_percent),
            title=args.fsi_title,
            path_color="lime",
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
        )

    print(f"wrote: {output_dir}")
    return 0


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--path-dat", required=True, help="Binary path.dat with lon+180, lat, r/Rsun.")
    parser.add_argument("--points", help="Optional pts.npy; first and last rows are used as footpoints.")
    parser.add_argument("--summary", help="Optional path_summary.json containing footpoints_stonyhurst.")
    parser.add_argument("--footpoints-stonyhurst", nargs=4, type=float, metavar=("LON0", "LAT0", "LON1", "LAT1"))
    parser.add_argument("--fsi", help="Solar Orbiter EUI/FSI FITS.")
    parser.add_argument("--euvi", help="STEREO EUVI FITS/FTS.")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--fsi-submap-arcsec", nargs=4, type=float, default=[2000.0, 3500.0, 250.0, 1750.0])
    parser.add_argument("--euvi-submap-arcsec", nargs=4, type=float, default=[-100.0, 410.0, 200.0, 500.0])
    parser.add_argument("--fsi-clip-percent", nargs=2, type=float, default=[10.0, 99.9])
    parser.add_argument("--euvi-clip-percent", nargs=2, type=float, default=[1.0, 99.9])
    parser.add_argument("--path-step", type=int, default=500)
    parser.add_argument("--fsi-title", default="Solar Orbiter FSI")
    parser.add_argument("--euvi-title", default="STEREO/EUVI")
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
        return np.vstack([points[0], points[-1]])
    return path[[0, -1], :2]


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
):
    raw_map = Map(map_path)
    smap = Map(np.asarray(raw_map.data, dtype=float), raw_map.meta).rotate(missing=0)
    submap = get_submap(smap, tuple(submap_arcsec))
    observer = submap.coordinate_frame.observer
    obs_lon = observer.lon.to_value(u.degree)
    obs_lat = observer.lat.to_value(u.degree)

    projected_path = project_lonlatr_to_observer_hpc(path, obs_lon, obs_lat, submap.rsun_obs.to_value(u.arcsec))
    footpoint_path = footpoints[:, :3]
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
    submap.plot(axes=ax, clip_interval=clip_percent * u.percent)
    ax.grid(alpha=0)
    set_solar_xy_labels(ax)
    ax.plot_coord(rope_coords, color=path_color, linestyle="--", linewidth=1.6, marker=None)
    for coord, color in zip(footpoint_coords, ["white", "black"], strict=True):
        ax.plot_coord(
            coord,
            "o",
            markersize=10,
            markerfacecolor=color,
            markeredgecolor=color,
            markeredgewidth=1.5,
            alpha=0.5,
        )
    wave = f"{submap.wavelength.to_value(u.AA):.0f} A" if submap.wavelength is not None else ""
    ax.set_title(f"{title} {wave} {submap.date.isot[:19]}")
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def project_lonlatr_to_observer_hpc(lon_lat_r, obs_lon, obs_lat, rsun_obs_arcsec):
    rotated = sun_rot_batch(np.asarray(lon_lat_r, dtype=float), obs_lon, obs_lat)
    return rotated[:, :2] * rsun_obs_arcsec


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

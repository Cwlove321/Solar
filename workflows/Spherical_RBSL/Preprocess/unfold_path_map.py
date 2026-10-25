#!/usr/bin/env python
"""Unfold an image along a projected 3D RBSL axis.

The horizontal coordinate is the physical arc length of the 3D axis.  The
vertical coordinate is the image-plane left normal to the projected axis when
walking from the first footpoint to the second footpoint.
"""

from __future__ import annotations

import argparse
import copy
import json
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.time import Time
from astropy.utils import iers
from astropy.visualization import ImageNormalize
from matplotlib.colors import LogNorm
from scipy.ndimage import map_coordinates
from sunpy.coordinates import frames
from sunpy.map import Map

iers.conf.auto_download = False


def cumulative_length_mm(path: np.ndarray, solar_radius_mm: float) -> np.ndarray:
    lon = np.radians(path[:, 0])
    lat = np.radians(path[:, 1])
    radius = path[:, 2] * solar_radius_mm
    xyz = np.column_stack(
        (
            radius * np.cos(lat) * np.cos(lon),
            radius * np.cos(lat) * np.sin(lon),
            radius * np.sin(lat),
        )
    )
    return np.concatenate(([0.0], np.cumsum(np.linalg.norm(np.diff(xyz, axis=0), axis=1))))


def resample_path(path: np.ndarray, distance: np.ndarray, step_mm: float):
    sample_distance = np.arange(0.0, distance[-1], step_mm)
    if sample_distance.size == 0 or sample_distance[-1] < distance[-1]:
        sample_distance = np.append(sample_distance, distance[-1])
    sampled = np.column_stack(
        [np.interp(sample_distance, distance, path[:, column]) for column in range(3)]
    )
    return sample_distance, sampled


def project_path(path: np.ndarray, smap, obstime: Time, solar_radius_mm: float):
    hgs = SkyCoord(
        path[:, 0] * u.deg,
        path[:, 1] * u.deg,
        radius=path[:, 2] * solar_radius_mm * u.Mm,
        frame=frames.HeliographicStonyhurst,
        obstime=obstime,
    )
    hpc = hgs.transform_to(smap.coordinate_frame)
    return np.column_stack((hpc.Tx.to_value(u.arcsec), hpc.Ty.to_value(u.arcsec)))


def unfold_map(
    map_path: Path,
    path: np.ndarray,
    *,
    path_obstime: Time,
    solar_radius_mm: float,
    half_width_mm: float,
    along_step_mm: float,
    normal_step_mm: float,
):
    smap = Map(map_path)
    distance = cumulative_length_mm(path, solar_radius_mm)
    x_mm, sampled_path = resample_path(path, distance, along_step_mm)
    projected = project_path(sampled_path, smap, path_obstime, solar_radius_mm)

    tangent = np.gradient(projected, x_mm, axis=0, edge_order=2)
    tangent /= np.linalg.norm(tangent, axis=1, keepdims=True)
    left_normal = np.column_stack((-tangent[:, 1], tangent[:, 0]))

    y_mm = np.arange(-half_width_mm, half_width_mm + 0.5 * normal_step_mm, normal_step_mm)
    arcsec_per_mm = smap.rsun_obs.to_value(u.arcsec) / solar_radius_mm
    sample_hpc = projected[None, :, :] + (
        y_mm[:, None, None] * arcsec_per_mm * left_normal[None, :, :]
    )
    world = SkyCoord(
        Tx=sample_hpc[..., 0] * u.arcsec,
        Ty=sample_hpc[..., 1] * u.arcsec,
        frame=smap.coordinate_frame,
    )
    px, py = smap.world_to_pixel(world)
    image = map_coordinates(
        np.asarray(smap.data, dtype=float),
        [py.to_value(u.pix), px.to_value(u.pix)],
        order=1,
        mode="constant",
        cval=np.nan,
        prefilter=False,
    )
    return smap, x_mm, y_mm, image, projected, left_normal, sample_hpc


def display_limits(image: np.ndarray, clip_percent: tuple[float, float]):
    finite = image[np.isfinite(image) & (image > 0)]
    if finite.size == 0:
        raise ValueError("The unfolded image has no finite positive samples")
    vmin, vmax = np.percentile(finite, clip_percent)
    if vmin <= 0 or vmax <= vmin:
        vmin = max(float(np.nanmin(finite)), np.finfo(float).tiny)
        vmax = float(np.nanmax(finite))
    return float(vmin), float(vmax)


def native_map_norm(smap, vmin: float, vmax: float):
    """Use the instrument's native SunPy stretch with explicit limits."""

    native = smap.plot_settings.get("norm")
    stretch = getattr(native, "stretch", None)
    if stretch is not None:
        return ImageNormalize(vmin=vmin, vmax=vmax, stretch=copy.deepcopy(stretch), clip=False)
    return LogNorm(vmin=vmin, vmax=vmax)


def plot_unfolded(smap, x_mm, y_mm, image, output, title, clip_percent):
    vmin, vmax = display_limits(image, clip_percent)
    cmap = smap.plot_settings.get("cmap", "gray")
    fig, ax = plt.subplots(figsize=(11.0, 4.2), constrained_layout=True)
    artist = ax.imshow(
        image,
        origin="lower",
        aspect="auto",
        extent=[x_mm[0], x_mm[-1], y_mm[0], y_mm[-1]],
        cmap=cmap,
        norm=native_map_norm(smap, vmin, vmax),
        interpolation="nearest",
    )
    ax.axhline(0.0, color="white", linestyle="--", linewidth=1.0, alpha=0.8)
    ax.set_xlabel("Distance along flux-rope axis [Mm]")
    ax.set_ylabel("Normal distance [Mm]")
    ax.set_title(title)
    colorbar = fig.colorbar(artist, ax=ax, pad=0.015)
    colorbar.set_label("Intensity [DN]")
    fig.savefig(output, dpi=200)
    plt.close(fig)
    return vmin, vmax


def plot_geometry(smap, sample_hpc, output, title, clip_percent):
    px, py = smap.world_to_pixel(
        SkyCoord(
            Tx=sample_hpc[..., 0] * u.arcsec,
            Ty=sample_hpc[..., 1] * u.arcsec,
            frame=smap.coordinate_frame,
        )
    )
    px = px.to_value(u.pix)
    py = py.to_value(u.pix)
    xmin, xmax = np.nanmin(px), np.nanmax(px)
    ymin, ymax = np.nanmin(py), np.nanmax(py)
    margin = 40
    x0 = max(0, int(np.floor(xmin - margin)))
    x1 = min(smap.data.shape[1], int(np.ceil(xmax + margin)))
    y0 = max(0, int(np.floor(ymin - margin)))
    y1 = min(smap.data.shape[0], int(np.ceil(ymax + margin)))
    data = np.asarray(smap.data[y0:y1, x0:x1], dtype=float)
    finite = data[np.isfinite(data) & (data > 0)]
    vmin, vmax = np.percentile(finite, clip_percent)
    fig, ax = plt.subplots(figsize=(9.0, 5.5), constrained_layout=True)
    ax.imshow(
        data,
        origin="lower",
        cmap=smap.plot_settings.get("cmap", "gray"),
        norm=native_map_norm(smap, float(vmin), float(vmax)),
    )
    for row, color, label in ((0, "cyan", "−20 Mm"), (sample_hpc.shape[0] // 2, "lime", "axis"), (-1, "cyan", "+20 Mm")):
        ax.plot(px[row] - x0, py[row] - y0, color=color, linewidth=1.5, label=label)
    ax.set_title(title)
    ax.set_xlabel("Image X [pixel]")
    ax.set_ylabel("Image Y [pixel]")
    ax.legend(loc="best", framealpha=0.7)
    fig.savefig(output, dpi=180)
    plt.close(fig)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--map", required=True)
    parser.add_argument("--path", required=True, help="Text file with lon_deg, lat_deg, radius/Rsun.")
    parser.add_argument("--path-obstime", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--name", required=True)
    parser.add_argument("--title", required=True)
    parser.add_argument("--solar-radius-mm", type=float, default=696.0)
    parser.add_argument("--half-width-mm", type=float, default=20.0)
    parser.add_argument("--along-step-mm", type=float, default=0.5)
    parser.add_argument("--normal-step-mm", type=float, default=0.5)
    parser.add_argument("--clip-percent", nargs=2, type=float, default=[1.0, 99.5])
    args = parser.parse_args(argv)

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = np.loadtxt(args.path)
    result = unfold_map(
        Path(args.map),
        path,
        path_obstime=Time(args.path_obstime),
        solar_radius_mm=args.solar_radius_mm,
        half_width_mm=args.half_width_mm,
        along_step_mm=args.along_step_mm,
        normal_step_mm=args.normal_step_mm,
    )
    smap, x_mm, y_mm, image, projected, left_normal, sample_hpc = result
    png = output_dir / f"{args.name}_axis_unfolded.png"
    geometry_png = output_dir / f"{args.name}_sampling_geometry.png"
    vmin, vmax = plot_unfolded(smap, x_mm, y_mm, image, png, args.title, tuple(args.clip_percent))
    plot_geometry(smap, sample_hpc, geometry_png, f"{args.title}: sampling band", tuple(args.clip_percent))
    np.savez_compressed(
        output_dir / f"{args.name}_axis_unfolded.npz",
        intensity=image,
        axis_distance_mm=x_mm,
        normal_distance_mm=y_mm,
        projected_axis_hpc_arcsec=projected,
        left_normal_hpc=left_normal,
        sample_hpc_arcsec=sample_hpc,
    )
    metadata = {
        "map": str(Path(args.map)),
        "path": str(Path(args.path)),
        "path_obstime": Time(args.path_obstime).isot,
        "solar_radius_mm": args.solar_radius_mm,
        "axis_length_mm": float(x_mm[-1]),
        "half_width_mm": args.half_width_mm,
        "along_step_mm": args.along_step_mm,
        "normal_step_mm": args.normal_step_mm,
        "shape_y_x": list(image.shape),
        "normal_definition": "+y is the image-plane left normal when walking from the first footpoint to the second.",
        "display_vmin_vmax": [vmin, vmax],
        "display_cmap": str(smap.plot_settings.get("cmap", "gray")),
        "display_stretch": type(getattr(smap.plot_settings.get("norm"), "stretch", None)).__name__,
        "display_clip_percent": list(args.clip_percent),
        "outputs": {"image": str(png), "geometry": str(geometry_png)},
    }
    (output_dir / f"{args.name}_axis_unfolded.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

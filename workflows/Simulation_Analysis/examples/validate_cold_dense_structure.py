#!/usr/bin/env python3
"""Validate whether a cold/dense structure follows magnetic field or AMR interfaces."""

from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, LogNorm, TwoSlopeNorm
from scipy.ndimage import distance_transform_edt

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import data_paths, load_config
from simulation_analysis.plotting import add_fieldlines, draw_scalar
from simulation_analysis.snapshot import extract_axis_slice, to_plot_plane


def radial_anomaly(values: np.ndarray) -> np.ndarray:
    reference = np.nanmedian(np.log10(np.maximum(values, 1.0e-300)), axis=1, keepdims=True)
    return np.log10(np.maximum(values, 1.0e-300)) - reference


def finite_summary(values: np.ndarray) -> dict[str, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return {"median": float("nan"), "p10": float("nan"), "p90": float("nan")}
    return {
        "median": float(np.median(finite)),
        "p10": float(np.percentile(finite, 10)),
        "p90": float(np.percentile(finite, 90)),
    }


def process(path: Path, config: dict, args) -> dict:
    section = extract_axis_slice(path, config, axis=args.axis, coordinate=args.coordinate)
    plane = to_plot_plane(section, config)
    temperature = plane.values["temperature_K"]
    density = plane.values["number_density_cm3"]
    pressure = plane.values["pressure_code"]
    entropy = plane.values["entropy_proxy"]
    dlogt = radial_anomaly(temperature)
    dlogn = radial_anomaly(density)
    dlogp = radial_anomaly(pressure)
    dlogs = radial_anomaly(entropy)
    yy = np.broadcast_to(plane.y[:, None], temperature.shape)
    domain = (yy >= args.min_height) & (yy <= args.max_height)
    structure = domain & (dlogt <= args.temperature_anomaly) & (dlogn >= args.density_anomaly)

    v = tuple(plane.values[name] for name in ("v1_code", "v2_code", "v3_code"))
    b = tuple(plane.values[name] for name in ("b1", "b2", "b3"))
    vmag = np.sqrt(sum(component**2 for component in v))
    bmag = np.sqrt(sum(component**2 for component in b))
    vdotb = sum(vc * bc for vc, bc in zip(v, b))
    angle = np.degrees(
        np.arccos(np.clip(np.abs(vdotb) / np.maximum(vmag * bmag, 1.0e-30), 0.0, 1.0))
    )
    angle[vmag < 1.0e-12] = np.nan

    level = plane.level.astype(float)
    interface = np.zeros_like(level, dtype=bool)
    interface[1:, :] |= level[1:, :] != level[:-1, :]
    interface[:-1, :] |= level[:-1, :] != level[1:, :]
    interface[:, 1:] |= level[:, 1:] != level[:, :-1]
    interface[:, :-1] |= level[:, :-1] != level[:, 1:]
    dy = float(np.median(np.diff(plane.y)))
    dx = float(np.median(np.diff(plane.x)))
    interface_distance = distance_transform_edt(~interface, sampling=(dy, dx))

    figure, axes = plt.subplots(2, 3, figsize=(17, 10), constrained_layout=True, sharex=True, sharey=True)
    panels = (
        (axes[0, 0], temperature, "Temperature [K]", "inferno", LogNorm(8.0e3, 4.0e6)),
        (axes[0, 1], density, r"Number density [cm$^{-3}$]", "viridis", LogNorm(2.0e9, 6.0e13)),
        (axes[0, 2], dlogs, r"Entropy anomaly [dex]", "coolwarm", TwoSlopeNorm(0.0, -1.5, 1.5)),
        (axes[1, 0], plane.values["vr_km_s"], r"$v_r$ [km s$^{-1}$]", "coolwarm", TwoSlopeNorm(0.0, -250.0, 250.0)),
        (axes[1, 1], angle, r"Acute angle$(v,B)$ [deg]", "viridis", plt.Normalize(0.0, 90.0)),
    )
    for axis, values, label, cmap, norm in panels:
        image = draw_scalar(axis, plane, values, cmap=cmap, norm=norm)
        axis.set_aspect("equal", adjustable="box")
        if np.any(structure):
            axis.contour(plane.x, plane.y, structure.astype(float), levels=[0.5], colors="cyan", linewidths=1.2)
        add_fieldlines(axis, plane)
        figure.colorbar(image, ax=axis, label=label, shrink=0.9)

    level_norm = BoundaryNorm(
        np.arange(int(np.nanmin(level)) - 0.5, int(np.nanmax(level)) + 1.5),
        int(np.nanmax(level) - np.nanmin(level) + 1),
    )
    image = draw_scalar(axes[1, 2], plane, level, cmap="viridis", norm=level_norm)
    axes[1, 2].set_aspect("equal", adjustable="box")
    if np.any(structure):
        axes[1, 2].contour(plane.x, plane.y, structure.astype(float), levels=[0.5], colors="red", linewidths=1.2)
    if np.any(interface):
        axes[1, 2].contour(plane.x, plane.y, interface.astype(float), levels=[0.5], colors="white", linewidths=0.45)
    figure.colorbar(image, ax=axes[1, 2], label="AMR level", shrink=0.9)

    titles = ("Temperature + B", "Density + B", "Entropy anomaly + B", "Radial velocity + B", "Velocity-field alignment + B", "AMR levels")
    for axis, title in zip(axes.flat, titles):
        axis.set_title(title)
        axis.set_xlabel(plane.x_name)
        axis.set_ylabel(plane.y_name)
        axis.set_ylim(args.min_plot_height, args.max_plot_height)
    figure.suptitle(f"{plane.fixed_label}, t={plane.time:.6g}; cyan/red=cold+dense mask")
    output_dir = Path(config["output_dir"]) / "cold_dense_validation_frames"
    output_dir.mkdir(parents=True, exist_ok=True)
    output = output_dir / f"cold_dense_validation_{path.stem.rsplit('_', 1)[-1]}.png"
    figure.savefig(output, dpi=int(config["plot"]["dpi"]), bbox_inches="tight")
    plt.close(figure)

    mask = structure & np.isfinite(angle)
    metrics = {
        "file": path.name,
        "time": plane.time,
        "pixels": int(np.count_nonzero(structure)),
        "temperature_K": finite_summary(temperature[structure]),
        "number_density_cm3": finite_summary(density[structure]),
        "pressure_anomaly_dex": finite_summary(dlogp[structure]),
        "entropy_anomaly_dex": finite_summary(dlogs[structure]),
        "radial_velocity_km_s": finite_summary(plane.values["vr_km_s"][structure]),
        "velocity_B_angle_degree": finite_summary(angle[mask]),
        "level4_fraction": float(np.mean(level[structure] == 4.0)) if np.any(structure) else float("nan"),
        "fraction_within_interface_distance": float(
            np.mean(interface_distance[structure] <= args.interface_distance)
        ) if np.any(structure) else float("nan"),
        "output": str(output),
    }
    del section, plane
    gc.collect()
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--axis", choices=("theta", "phi"))
    parser.add_argument("--coordinate", type=float)
    parser.add_argument("--all", action="store_true")
    parser.add_argument("--temperature-anomaly", type=float, default=-0.10)
    parser.add_argument("--density-anomaly", type=float, default=0.10)
    parser.add_argument("--min-height", type=float, default=1.0)
    parser.add_argument("--max-height", type=float, default=40.0)
    parser.add_argument("--min-plot-height", type=float, default=0.0)
    parser.add_argument("--max-plot-height", type=float, default=60.0)
    parser.add_argument("--interface-distance", type=float, default=0.75)
    args = parser.parse_args()
    config = load_config(args.config)
    paths = data_paths(config)
    if not paths:
        raise SystemExit("No DAT files selected")
    if not args.all:
        paths = [paths[-1]]
    metrics = []
    for position, path in enumerate(paths, start=1):
        result = process(path, config, args)
        metrics.append(result)
        print(f"[{position}/{len(paths)}] {path.name}: pixels={result['pixels']}", flush=True)
    output = Path(config["output_dir"]) / "cold_dense_validation_metrics.json"
    output.write_text(json.dumps(metrics, indent=2, allow_nan=True) + "\n")
    print(output, flush=True)


if __name__ == "__main__":
    main()

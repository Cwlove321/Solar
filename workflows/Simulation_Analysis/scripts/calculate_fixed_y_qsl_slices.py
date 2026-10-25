#!/usr/bin/env python3
"""Calculate and plot FastQSL2 vertical x-z slices at an exact fixed y.

The default run computes the AMRVAC-normalized plane y=0 for snapshots 10 and
95.  It writes one independently reusable NPZ per snapshot, a metadata JSON,
and full-domain plus flux-rope-region figures.  The magnetic axis intersection
is marked with an X (not a magnetic-null star) when an axis-feature file is
available.

Example
-------
/opt/anaconda3/bin/python \
  /Users/xiaoxingyun/Desktop/solar_code/workflows/Simulation_Analysis/scripts/\
calculate_fixed_y_qsl_slices.py \
  --magnetic-dir /Volumes/Seagate/w2_jet_mhd/out_path1/nlfff \
  --output-dir /Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/y0_qsl \
  --fixed-y 0 --threads 8
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import tempfile
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


REPOSITORY_ROOT = Path(__file__).resolve().parents[3]
FASTQSL2_DIR = Path(
    os.environ.get("FASTQSL2_DIR", REPOSITORY_ROOT / "Packages" / "FastQSL2")
).expanduser()
MAGNETIC_DIR = Path("/Volumes/Seagate/w2_jet_mhd/out_path1/nlfff")
OUTPUT_DIR = Path("/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/y0_qsl")
AXIS_FEATURE_DIR = Path(
    "/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/bottom_fastqsl2/"
    "overlying_flux_candidate"
)
SNAPSHOTS = (10, 95)

# AMRVAC normalized coordinates.  The 128^3 magnetic field is cell centred,
# one normalized length unit is 10 Mm, and z=0.1 is the bottom-boundary datum.
DOMAIN = {
    "x": (-4.0, 8.0),
    "y": (-6.0, 6.0),
    "z": (0.1, 12.1),
}
LENGTH_UNIT_MM = 10.0
BOTTOM_Z_CODE = 0.1
TIME_MINUTES_PER_SNAPSHOT = 0.71558

LOGQ_LEVELS = (2.0, 3.0, 4.0, 5.0, 6.0)
Q_COLORS = ("#E69F00", "#009E73", "#CC79A7", "#D55E00", "#0072B2")
Q_LINESTYLES = ("solid", "dashdot", "dotted", (0, (5, 1)), (0, (3, 1, 1, 1)))
CLOSED_COLOR = "#6A3D9A"
AXIS_COLOR = "#00BFC4"
ZOOM_LIMITS_MM = (-15.0, 30.0, 0.0, 30.0)

if str(FASTQSL2_DIR) not in sys.path:
    sys.path.insert(0, str(FASTQSL2_DIR))

from fastqsl import fastqsl  # noqa: E402


def read_binary_data(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read a QSL magnetic binary and return Bx, By, Bz in (z,y,x) order."""

    with path.open("rb") as stream:
        header = np.fromfile(stream, dtype=np.int32, count=3)
        if header.size != 3:
            raise ValueError(f"Missing (nx,ny,nz) header in {path}")
        nx, ny, nz = (int(value) for value in header)
        count = nx * ny * nz
        fields = []
        for component in ("Bx", "By", "Bz"):
            values = np.fromfile(stream, dtype=np.float64, count=count)
            if values.size != count:
                raise ValueError(f"Incomplete {component} array in {path}")
            fields.append(values.reshape((nz, ny, nx)))
        if stream.read(1):
            raise ValueError(f"Unexpected trailing bytes in {path}")
    return fields[0], fields[1], fields[2]


def cell_centres(bounds: tuple[float, float], count: int) -> np.ndarray:
    spacing = (bounds[1] - bounds[0]) / count
    return bounds[0] + (np.arange(count, dtype=float) + 0.5) * spacing


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def finite_stats(array: np.ndarray) -> dict[str, object]:
    finite = np.isfinite(array)
    values = array[finite]
    return {
        "shape": list(array.shape),
        "finite_count": int(finite.sum()),
        "nan_count": int(np.isnan(array).sum()),
        "inf_count": int(np.isinf(array).sum()),
        "min": float(values.min()) if values.size else None,
        "median": float(np.median(values)) if values.size else None,
        "max": float(values.max()) if values.size else None,
    }


def magnetic_axis_intersections(feature_path: Path, fixed_y_code: float) -> np.ndarray:
    """Find intersections between the stored magnetic axis and fixed-y plane."""

    if not feature_path.is_file():
        return np.empty((0, 2), dtype=float)
    with np.load(feature_path) as archive:
        if "axis" not in archive.files:
            return np.empty((0, 2), dtype=float)
        axis = np.asarray(archive["axis"], dtype=float)
    if axis.ndim != 2 or axis.shape[1] != 3 or len(axis) < 2:
        raise ValueError(f"Unexpected magnetic-axis shape in {feature_path}: {axis.shape}")

    # The feature file uses the old 0..128 index coordinates.  Convert the
    # requested AMRVAC y coordinate into that index system before intersecting.
    native_y = 128.0 * (fixed_y_code - DOMAIN["y"][0]) / (DOMAIN["y"][1] - DOMAIN["y"][0])
    hits = []
    for p0, p1 in zip(axis[:-1], axis[1:]):
        y0 = p0[1] - native_y
        y1 = p1[1] - native_y
        if np.isclose(y0, 0.0) and np.isclose(y1, 0.0):
            hits.extend((p0, p1))
            continue
        if np.isclose(y0, 0.0):
            hits.append(p0)
            continue
        if np.isclose(y1, 0.0):
            hits.append(p1)
            continue
        if y0 * y1 < 0.0:
            fraction = -y0 / (y1 - y0)
            hits.append(p0 + fraction * (p1 - p0))
    if not hits:
        return np.empty((0, 2), dtype=float)

    native_spacing_mm = (DOMAIN["x"][1] - DOMAIN["x"][0]) * LENGTH_UNIT_MM / 128.0
    converted = np.asarray(
        [
            [DOMAIN["x"][0] * LENGTH_UNIT_MM + hit[0] * native_spacing_mm,
             hit[2] * native_spacing_mm]
            for hit in hits
        ],
        dtype=float,
    )
    # A crossing exactly at an axis vertex is discovered by its two adjacent
    # segments; retain it once.
    rounded = np.round(converted, decimals=8)
    _, unique_indices = np.unique(rounded, axis=0, return_index=True)
    return converted[np.sort(unique_indices)]


def calculate_slice(
    snapshot: int,
    magnetic_path: Path,
    feature_path: Path,
    output_dir: Path,
    fixed_y_code: float,
    delta_index: float,
    threads: int,
) -> tuple[dict, dict[str, np.ndarray]]:
    """Calculate one exact fixed-y FastQSL2 x-z slice."""

    bx, by, bz = read_binary_data(magnetic_path)
    nz, ny, nx = bx.shape
    x = cell_centres(DOMAIN["x"], nx)
    y = cell_centres(DOMAIN["y"], ny)
    z = cell_centres(DOMAIN["z"], nz)
    dx = float(x[1] - x[0])
    if not np.allclose(np.diff(y), dx) or not np.allclose(np.diff(z), dx):
        raise ValueError("The QSL cube is expected to have equal grid spacing")
    if not (DOMAIN["y"][0] <= fixed_y_code <= DOMAIN["y"][1]):
        raise ValueError(f"fixed y={fixed_y_code} lies outside {DOMAIN['y']}")

    delta_code = delta_index * dx
    with tempfile.TemporaryDirectory(prefix=f"fastqsl2-y{fixed_y_code:+.3f}-t{snapshot:04d}-") as tmp_dir:
        result = fastqsl(
            bx,
            by,
            bz,
            xa=x,
            ya=y,
            za=z,
            xreg=[x[0], x[-1]],
            yreg=[fixed_y_code, fixed_y_code],
            zreg=[z[0], z[-1]],
            delta=delta_code,
            seed=True,
            nthreads=threads,
            silent=True,
            tmp_dir=tmp_dir,
        )

    q = np.asarray(result["q"], dtype=float)
    rboundary = np.asarray(result["rboundary"], dtype=np.int16)
    seed = np.asarray(result["seed"], dtype=float)
    if q.shape != rboundary.shape or seed.shape != q.shape + (3,):
        raise ValueError(f"Unexpected FastQSL2 output shapes at t={snapshot}")
    logq = np.full(q.shape, np.nan)
    valid = np.isfinite(q) & (q > 0.0)
    logq[valid] = np.log10(q[valid])
    x_mm = seed[..., 0] * LENGTH_UNIT_MM
    height_mm = (seed[..., 2] - BOTTOM_Z_CODE) * LENGTH_UNIT_MM
    y_actual = seed[..., 1]
    if not np.allclose(y_actual, fixed_y_code, atol=1.0e-12):
        raise ValueError("FastQSL2 seed plane differs from the requested fixed y")
    axis_hits = magnetic_axis_intersections(feature_path, fixed_y_code)

    result["logq"] = logq
    result["x_mm"] = x_mm
    result["height_mm"] = height_mm
    result["fixed_y_code_actual"] = y_actual
    result["fixed_y_Mm"] = np.asarray(fixed_y_code * LENGTH_UNIT_MM)
    result["magnetic_axis_intersections_Mm"] = axis_hits

    y_tag = f"{fixed_y_code:+.6f}".replace("+", "p").replace("-", "m").replace(".", "p")
    output_file = output_dir / f"fixed_y_qsl_t{snapshot:04d}_y{y_tag}_delta{delta_index:.2f}.npz"
    np.savez_compressed(output_file, **result)

    summary = {
        "snapshot": snapshot,
        "time_minutes": snapshot * TIME_MINUTES_PER_SNAPSHOT,
        "magnetic_binary": str(magnetic_path.resolve()),
        "magnetic_binary_sha256": sha256(magnetic_path),
        "magnetic_axis_source": str(feature_path.resolve()) if feature_path.is_file() else None,
        "magnetic_axis_intersections_Mm": axis_hits.tolist(),
        "slice_plane": "x-z at exact fixed y",
        "fixed_y_amrvac_normalized_requested": fixed_y_code,
        "fixed_y_amrvac_normalized_actual_minmax": [
            float(np.nanmin(y_actual)),
            float(np.nanmax(y_actual)),
        ],
        "fixed_y_Mm": fixed_y_code * LENGTH_UNIT_MM,
        "input_shape_z_y_x": [nz, ny, nx],
        "cell_center_ranges_amrvac_normalized": {
            "x": [float(x[0]), float(x[-1])],
            "y": [float(y[0]), float(y[-1])],
            "z_absolute": [float(z[0]), float(z[-1])],
        },
        "plotted_ranges_Mm": {
            "x": [float(np.nanmin(x_mm)), float(np.nanmax(x_mm))],
            "height": [float(np.nanmin(height_mm)), float(np.nanmax(height_mm))],
        },
        "grid_spacing_amrvac_normalized": dx,
        "fastqsl2_delta_original_cells": delta_index,
        "fastqsl2_delta_amrvac_normalized": delta_code,
        "fastqsl2_delta_Mm": delta_code * LENGTH_UNIT_MM,
        "q": finite_stats(q),
        "log10_q": finite_stats(logq),
        "closed_bottom_to_bottom_count": int(np.count_nonzero(rboundary == 11)),
        "closed_bottom_to_bottom_fraction": float(np.mean(rboundary == 11)),
        "output_npz": str(output_file.resolve()),
    }
    return summary, result


def draw_panel(axis, summary: dict, result: dict, *, limits, logq_min: float, logq_max: float):
    """Draw one Q slice with Q contours, closed-field boundary, and axis hit."""

    x_mm = result["x_mm"]
    height_mm = result["height_mm"]
    logq = result["logq"]
    closed = np.asarray(result["rboundary"]) == 11
    image = axis.pcolormesh(
        x_mm,
        height_mm,
        logq,
        shading="nearest",
        cmap="gray_r",
        vmin=logq_min,
        vmax=logq_max,
        rasterized=True,
        zorder=1,
    )
    closed_contour = axis.contour(
        x_mm,
        height_mm,
        closed.astype(float),
        levels=[0.5],
        colors=[CLOSED_COLOR],
        linewidths=1.2,
        zorder=2,
    )
    contour_counts = {}
    for level, color, style in zip(LOGQ_LEVELS, Q_COLORS, Q_LINESTYLES):
        if level == 2.0:
            axis.contour(
                x_mm,
                height_mm,
                logq,
                levels=[level],
                colors=["white"],
                linewidths=3.4,
                zorder=2.8,
            )
        contour = axis.contour(
            x_mm,
            height_mm,
            logq,
            levels=[level],
            colors=[color],
            linestyles=[style],
            linewidths=1.8 if level == 2.0 else 1.3,
            zorder=3,
        )
        contour_counts[str(int(level))] = len(contour.allsegs[0])

    axis_hits = np.asarray(result["magnetic_axis_intersections_Mm"], dtype=float)
    if axis_hits.size:
        axis.scatter(
            axis_hits[:, 0],
            axis_hits[:, 1],
            marker="X",
            s=105,
            facecolor=AXIS_COLOR,
            edgecolor="#101010",
            linewidth=0.9,
            zorder=5,
        )
    axis.set_title(f"t = {summary['time_minutes']:.1f} min", fontsize=17, pad=8)
    axis.set_xlim(limits[:2])
    axis.set_ylim(limits[2:])
    axis.set_aspect("equal")
    axis.tick_params(direction="in", top=True, right=True, labelsize=12)
    return image, len(closed_contour.allsegs[0]), contour_counts


def plot_slices(
    datasets: list[tuple[dict, dict[str, np.ndarray]]],
    output_base: Path,
    limits: tuple[float, float, float, float],
    logq_min: float,
    logq_max: float,
    dpi: int,
) -> list[dict]:
    """Create a shared-axis two-snapshot figure."""

    figure, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(17.6, 8.2),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    panel_metadata = []
    image = None
    for index, (axis, (summary, result)) in enumerate(zip(axes, datasets)):
        image, closed_count, q_counts = draw_panel(
            axis,
            summary,
            result,
            limits=limits,
            logq_min=logq_min,
            logq_max=logq_max,
        )
        axis.set_xlabel("Solar-X [Mm]", fontsize=14)
        if index == 0:
            axis.set_ylabel("Height above bottom boundary [Mm]", fontsize=14)
        panel_metadata.append(
            {"closed_contour_path_count": closed_count, "q_contour_path_counts": q_counts}
        )

    handles = [
        Line2D([0], [0], color=CLOSED_COLOR, linewidth=2.0, label="Both terminals at bottom"),
        Line2D(
            [0], [0], marker="X", linestyle="none", markersize=9,
            markerfacecolor=AXIS_COLOR, markeredgecolor="#101010",
            label="Magnetic-axis intersection",
        ),
    ]
    handles.extend(
        Line2D(
            [0], [0], color=color, linestyle=style, linewidth=2.1,
            label=rf"$\log_{{10}}Q={int(level)}$",
        )
        for level, color, style in zip(LOGQ_LEVELS, Q_COLORS, Q_LINESTYLES)
    )
    figure.legend(
        handles=handles,
        loc="outside lower center",
        ncol=4,
        frameon=False,
        handlelength=3.0,
        columnspacing=1.25,
        fontsize=11,
    )
    colorbar = figure.colorbar(image, ax=axes, shrink=0.91, pad=0.02)
    colorbar.set_label(r"$\log_{10}Q$", fontsize=14)
    colorbar.ax.tick_params(labelsize=11)
    figure.suptitle(
        rf"FastQSL2 vertical slice at $y={datasets[0][0]['fixed_y_amrvac_normalized_requested']:g}$",
        fontsize=18,
    )
    figure.savefig(output_base.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    figure.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return panel_metadata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--magnetic-dir", type=Path, default=MAGNETIC_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--axis-feature-dir", type=Path, default=AXIS_FEATURE_DIR)
    parser.add_argument("--snapshots", type=int, nargs="+", default=list(SNAPSHOTS))
    parser.add_argument("--fixed-y", type=float, default=0.0,
                        help="Fixed y in AMRVAC normalized coordinates (1 unit = 10 Mm).")
    parser.add_argument("--delta", type=float, default=0.25,
                        help="FastQSL2 sampling interval in original magnetic-grid cells.")
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--logq-min", type=float, default=float(np.log10(2.0)))
    parser.add_argument("--logq-max", type=float, default=8.0)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if not np.isfinite(args.fixed_y):
        parser.error("--fixed-y must be finite")
    if not np.isfinite(args.delta) or args.delta <= 0.0:
        parser.error("--delta must be positive and finite")
    if args.logq_min >= args.logq_max:
        parser.error("--logq-min must be smaller than --logq-max")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    datasets = []
    for snapshot in args.snapshots:
        magnetic_path = args.magnetic_dir / f"nlfff_{snapshot:04d}.binary"
        feature_path = args.axis_feature_dir / f"axis_features_factor2_t{snapshot:04d}.npz"
        if not magnetic_path.is_file():
            raise FileNotFoundError(magnetic_path)
        print(f"Calculating t={snapshot} at fixed y={args.fixed_y:g}...", flush=True)
        datasets.append(
            calculate_slice(
                snapshot,
                magnetic_path,
                feature_path,
                args.output_dir,
                args.fixed_y,
                args.delta,
                args.threads,
            )
        )

    snapshot_tag = "_".join(f"t{snapshot:04d}" for snapshot in args.snapshots)
    y_tag = f"{args.fixed_y:+.3f}".replace("+", "p").replace("-", "m").replace(".", "p")
    full_base = args.output_dir / f"fixed_y_qsl_y{y_tag}_{snapshot_tag}_full"
    zoom_base = args.output_dir / f"fixed_y_qsl_y{y_tag}_{snapshot_tag}_zoom"
    first_result = datasets[0][1]
    full_limits = (
        float(np.nanmin(first_result["x_mm"])),
        float(np.nanmax(first_result["x_mm"])),
        0.0,
        float(np.nanmax(first_result["height_mm"])),
    )
    full_panels = plot_slices(
        datasets, full_base, full_limits, args.logq_min, args.logq_max, args.dpi
    )
    zoom_panels = plot_slices(
        datasets, zoom_base, ZOOM_LIMITS_MM, args.logq_min, args.logq_max, args.dpi
    )

    metadata = {
        "method": "FastQSL2 CPU/OpenMP at an exact fixed-y x-z plane",
        "fastqsl2_directory": str(FASTQSL2_DIR.resolve()),
        "coordinate_definition": {
            "length_unit_Mm": LENGTH_UNIT_MM,
            "bottom_z_amrvac_normalized": BOTTOM_Z_CODE,
            "magnetic_grid": "128^3 cell centered",
            "fixed_y_amrvac_normalized": args.fixed_y,
            "fixed_y_Mm": args.fixed_y * LENGTH_UNIT_MM,
            "height": "(z_absolute - 0.1) * 10 Mm",
            "lowest_Q_seed_height_Mm": float(np.nanmin(first_result["height_mm"])),
        },
        "plot": {
            "quantity": "log10(Q)",
            "color_limits": [args.logq_min, args.logq_max],
            "contour_levels_log10Q": list(LOGQ_LEVELS),
            "closed_contour": "rboundary==11; both field-line terminals at bottom",
            "axis_marker": "X; magnetic-axis intersection, not a magnetic null",
            "full_limits_Mm": list(full_limits),
            "zoom_limits_Mm": list(ZOOM_LIMITS_MM),
            "full_png": str(full_base.with_suffix(".png").resolve()),
            "full_pdf": str(full_base.with_suffix(".pdf").resolve()),
            "zoom_png": str(zoom_base.with_suffix(".png").resolve()),
            "zoom_pdf": str(zoom_base.with_suffix(".pdf").resolve()),
        },
        "snapshots": [],
    }
    for (summary, _), full_panel, zoom_panel in zip(datasets, full_panels, zoom_panels):
        summary["plot_contours"] = {"full": full_panel, "zoom": zoom_panel}
        metadata["snapshots"].append(summary)

    metadata_path = args.output_dir / f"fixed_y_qsl_y{y_tag}_{snapshot_tag}_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(full_base.with_suffix(".png").resolve())
    print(zoom_base.with_suffix(".png").resolve())
    print(metadata_path.resolve())


if __name__ == "__main__":
    main()

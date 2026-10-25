#!/usr/bin/env python3
"""Calculate FastQSL2 x-z slices at the y coordinate of coronal nulls."""

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
JET_MHD_ROOT = Path(
    os.environ.get("JET_MHD_ROOT", REPOSITORY_ROOT.parent / "jet_mhd")
).expanduser()
MAGNETIC_DIR = Path(
    os.environ.get("FASTQSL_MAGNETIC_DIR", JET_MHD_ROOT / "out_path1" / "nlfff")
).expanduser()
NULL_DIR = Path(
    os.environ.get("FASTQSL_NULL_DIR", JET_MHD_ROOT / "qfactor" / "magnetic_nulls")
).expanduser()
OUTPUT_DIR = Path(
    os.environ.get("FASTQSL_OUTPUT_DIR", JET_MHD_ROOT / "qfactor" / "null_y_slices")
).expanduser()
FASTQSL2_DIR = Path(
    os.environ.get("FASTQSL2_DIR", REPOSITORY_ROOT / "Packages" / "FastQSL2")
).expanduser()
SNAPSHOTS = (10, 95)

# AMRVAC normalized coordinates of the 128^3 QSL cube.  The magnetic field is
# cell centred.  One normalized length unit is 10 Mm.
DOMAIN = {
    "x": (-4.0, 8.0),
    "y": (-6.0, 6.0),
    "z": (0.1, 12.1),
}
LENGTH_UNIT_MM = 10.0
BOTTOM_Z_CODE = 0.1
TIME_MINUTES_PER_SNAPSHOT = 0.71558

if str(FASTQSL2_DIR) not in sys.path:
    sys.path.insert(0, str(FASTQSL2_DIR))

from fastqsl import fastqsl  # noqa: E402


def read_binary_data(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read the QSL magnetic binary and return arrays in (z, y, x) order."""
    with path.open("rb") as stream:
        nx, ny, nz = (int(value) for value in np.fromfile(stream, np.int32, 3))
        count = nx * ny * nz
        fields = []
        for _ in range(3):
            values = np.fromfile(stream, np.float64, count)
            if values.size != count:
                raise ValueError(f"Incomplete magnetic field in {path}")
            fields.append(values.reshape((nz, ny, nx)))
    return fields[0], fields[1], fields[2]


def cell_centres(bounds: tuple[float, float], count: int) -> np.ndarray:
    spacing = (bounds[1] - bounds[0]) / count
    return bounds[0] + (np.arange(count, dtype=float) + 0.5) * spacing


def read_coronal_null(path: Path) -> tuple[dict, dict]:
    """Return the highest validated null and the complete null-search record."""
    record = json.loads(path.read_text())
    if not record["nulls"]:
        raise ValueError(f"No magnetic nulls found in {path}")
    return max(record["nulls"], key=lambda item: item["position"]["z"]), record


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


def calculate_slice(
    snapshot: int,
    magnetic_path: Path,
    null_path: Path,
    output_dir: Path,
    delta_index: float,
    threads: int,
) -> tuple[dict, dict[str, np.ndarray]]:
    bx, by, bz = read_binary_data(magnetic_path)
    nz, ny, nx = bx.shape
    x = cell_centres(DOMAIN["x"], nx)
    y = cell_centres(DOMAIN["y"], ny)
    z = cell_centres(DOMAIN["z"], nz)
    dx = float(x[1] - x[0])
    if not np.allclose([np.diff(y), np.diff(z)], dx):
        raise ValueError("The QSL cube is expected to have equal grid spacing")

    null, null_record = read_coronal_null(null_path)
    null_mm = null["position"]
    null_code = {
        "x": null_mm["x"] / LENGTH_UNIT_MM,
        "y": null_mm["y"] / LENGTH_UNIT_MM,
        "z": BOTTOM_Z_CODE + null_mm["z"] / LENGTH_UNIT_MM,
    }
    if not (y[0] <= null_code["y"] <= y[-1]):
        raise ValueError(f"Null y={null_code['y']} is outside the magnetic cube")

    delta_code = delta_index * dx
    with tempfile.TemporaryDirectory(prefix=f"fastqsl2-null-y-t{snapshot:04d}-") as tmp_dir:
        result = fastqsl(
            bx,
            by,
            bz,
            xa=x,
            ya=y,
            za=z,
            xreg=[x[0], x[-1]],
            yreg=[null_code["y"], null_code["y"]],
            zreg=[z[0], z[-1]],
            delta=delta_code,
            seed=True,
            nthreads=threads,
            silent=True,
            tmp_dir=tmp_dir,
        )

    q = np.asarray(result["q"], dtype=float)
    logq = np.full(q.shape, np.nan)
    valid = np.isfinite(q) & (q > 0.0)
    logq[valid] = np.log10(q[valid])
    result["logq"] = logq

    seed = np.asarray(result["seed"], dtype=float)
    x_mm = seed[..., 0] * LENGTH_UNIT_MM
    y_code_actual = seed[..., 1]
    height_mm = (seed[..., 2] - BOTTOM_Z_CODE) * LENGTH_UNIT_MM
    result["x_mm"] = x_mm
    result["height_mm"] = height_mm
    result["y_code_actual"] = y_code_actual
    result["null_position_code"] = np.array(
        [null_code["x"], null_code["y"], null_code["z"]], dtype=float
    )
    result["null_position_mm"] = np.array(
        [null_mm["x"], null_mm["y"], null_mm["z"]], dtype=float
    )

    output_file = output_dir / (
        f"fixed_y_qsl_t{snapshot:04d}_y{null_code['y']:.6f}_delta{delta_index:.2f}.npz"
    )
    np.savez_compressed(output_file, **result)

    summary = {
        "snapshot": snapshot,
        "time_minutes": snapshot * TIME_MINUTES_PER_SNAPSHOT,
        "magnetic_binary": str(magnetic_path.resolve()),
        "magnetic_binary_sha256": sha256(magnetic_path),
        "null_json": str(null_path.resolve()),
        "null_json_sha256": sha256(null_path),
        "null_selection": "highest validated null (the coronal-null star)",
        "null_count_in_source": int(null_record["null_count"]),
        "null_index_position": null["index_position"],
        "null_position_amrvac_normalized": null_code,
        "null_position_Mm": {
            "x": null_mm["x"],
            "y": null_mm["y"],
            "height": null_mm["z"],
        },
        "slice_plane": "x-z at exact fixed y",
        "fixed_y_amrvac_normalized_requested": null_code["y"],
        "fixed_y_amrvac_normalized_actual_minmax": [
            float(np.nanmin(y_code_actual)),
            float(np.nanmax(y_code_actual)),
        ],
        "fixed_y_Mm": null_mm["y"],
        "input_shape_z_y_x": [nz, ny, nx],
        "cell_center_ranges_amrvac_normalized": {
            "x": [float(x[0]), float(x[-1])],
            "y": [float(y[0]), float(y[-1])],
            "z_absolute": [float(z[0]), float(z[-1])],
        },
        "grid_spacing_amrvac_normalized": dx,
        "fastqsl2_delta_original_cells": delta_index,
        "fastqsl2_delta_amrvac_normalized": delta_code,
        "fastqsl2_delta_Mm": delta_code * LENGTH_UNIT_MM,
        "q": finite_stats(q),
        "log10_q": finite_stats(logq),
        "output_npz": str(output_file.resolve()),
    }
    return summary, result


def plot_slices(
    datasets: list[tuple[dict, dict[str, np.ndarray]]],
    output_path: Path,
    logq_min: float,
    logq_max: float,
) -> None:
    figure, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(13.2, 6.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    image = None
    for panel, (axis, (summary, result)) in enumerate(zip(axes, datasets)):
        image = axis.pcolormesh(
            result["x_mm"],
            result["height_mm"],
            result["logq"],
            shading="nearest",
            cmap="gray_r",
            vmin=logq_min,
            vmax=logq_max,
            rasterized=True,
        )
        null = summary["null_position_Mm"]
        axis.scatter(
            null["x"],
            null["height"],
            marker="*",
            s=250,
            facecolor="#3ee7e7",
            edgecolor="#101010",
            linewidth=1.1,
            zorder=5,
        )
        axis.set_title(
            f"({chr(97 + panel)})  t = {summary['time_minutes']:.1f} min\n"
            f"fixed y = {summary['fixed_y_amrvac_normalized_requested']:.5f} "
            f"({summary['fixed_y_Mm']:.2f} Mm)"
        )
        axis.set_xlabel("Solar-X [Mm]")
        axis.tick_params(direction="in", top=True, right=True)
        axis.set_aspect("equal")
    axes[0].set_ylabel("Height above bottom boundary [Mm]")

    legend = Line2D(
        [0],
        [0],
        marker="*",
        linestyle="none",
        markersize=13,
        markerfacecolor="#3ee7e7",
        markeredgecolor="#101010",
        label="Magnetic null",
    )
    axes[0].legend(handles=[legend], loc="upper right", fontsize=9, framealpha=0.92)
    colorbar = figure.colorbar(image, ax=axes, shrink=0.91, pad=0.02)
    colorbar.set_label(r"$\log_{10} Q$")
    figure.suptitle("FastQSL2 vertical slices through the coronal magnetic nulls")
    figure.savefig(output_path, dpi=240, bbox_inches="tight")
    figure.savefig(output_path.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--magnetic-dir", type=Path, default=MAGNETIC_DIR)
    parser.add_argument("--null-dir", type=Path, default=NULL_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--snapshots", type=int, nargs="+", default=list(SNAPSHOTS))
    parser.add_argument("--delta", type=float, default=0.25)
    parser.add_argument("--threads", type=int, default=0)
    parser.add_argument("--logq-min", type=float, default=0.3)
    parser.add_argument("--logq-max", type=float, default=8.0)
    args = parser.parse_args()
    if args.delta <= 0.0:
        raise ValueError("--delta must be positive")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    datasets = []
    for snapshot in args.snapshots:
        print(f"Calculating t={snapshot} fixed-y Q slice...", flush=True)
        datasets.append(
            calculate_slice(
                snapshot,
                args.magnetic_dir / f"nlfff_{snapshot:04d}.binary",
                args.null_dir / f"magnetic_nulls_t{snapshot:04d}.json",
                args.output_dir,
                args.delta,
                args.threads,
            )
        )

    plot_path = args.output_dir / "fixed_y_qsl_at_coronal_nulls_t0010_t0095.png"
    plot_slices(datasets, plot_path, args.logq_min, args.logq_max)
    metadata = {
        "method": "FastQSL2 CPU/OpenMP on interpolated AMRVAC-normalized coordinates",
        "fastqsl2_directory": str(FASTQSL2_DIR.resolve()),
        "coordinate_definition": {
            "length_unit_Mm": LENGTH_UNIT_MM,
            "bottom_z_amrvac_normalized": BOTTOM_Z_CODE,
            "magnetic_grid": "cell centered",
            "slice": "vertical x-z plane at each coronal null's exact y",
        },
        "plot": {
            "png": str(plot_path.resolve()),
            "pdf": str(plot_path.with_suffix(".pdf").resolve()),
            "quantity": "log10(Q)",
            "color_limits": [args.logq_min, args.logq_max],
            "null_markers": "only the selected coronal null is shown, as a star",
        },
        "snapshots": [item[0] for item in datasets],
    }
    metadata_path = args.output_dir / "fixed_y_qsl_at_coronal_nulls_metadata.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(f"Wrote {plot_path}")
    print(f"Wrote {metadata_path}")


if __name__ == "__main__":
    main()

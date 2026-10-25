#!/usr/bin/env python3
"""Plot two exact-y=0 |curl(B)|/|B| slices with one shared color scale."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import numpy as np


DEFAULT_PRODUCT_DIR = Path(
    "/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/y0_qsl/j_over_b_evolution"
)
DEFAULT_STACK = DEFAULT_PRODUCT_DIR / "y0_j_over_b_t0010_t0095_stack.npz"
DEFAULT_METADATA = DEFAULT_PRODUCT_DIR / "y0_j_over_b_t0010_t0095_metadata.json"
DEFAULT_QSL_DIR = DEFAULT_PRODUCT_DIR.parent
DEFAULT_OUTPUT_STEM = (
    DEFAULT_PRODUCT_DIR
    / "y0_j_over_b_t0010_t0095_side_by_side_logq124_current_sheet_peaks"
)
DEFAULT_CONTOUR_LEVELS = (1.0, 2.0, 4.0)
CONTOUR_COLORS = ("#0072B2", "#E69F00", "#CC79A7")
CONTOUR_LINESTYLES = ("solid", "dashed", "dashdot")
DEFAULT_PEAK_WINDOW_MM = (12.0, 18.0, 8.0, 60.0)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack", type=Path, default=DEFAULT_STACK)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--qsl-dir", type=Path, default=DEFAULT_QSL_DIR)
    parser.add_argument("--snapshots", type=int, nargs=2, default=(10, 95))
    parser.add_argument(
        "--contour-levels",
        type=float,
        nargs="+",
        default=DEFAULT_CONTOUR_LEVELS,
        help="log10(Q) contour levels (default: 1 2 4)",
    )
    parser.add_argument(
        "--peak-window",
        type=float,
        nargs=4,
        metavar=("XMIN", "XMAX", "HMIN", "HMAX"),
        default=DEFAULT_PEAK_WINDOW_MM,
        help="window in Mm used to locate the current-sheet J/B maximum",
    )
    parser.add_argument("--output-stem", type=Path, default=DEFAULT_OUTPUT_STEM)
    parser.add_argument("--dpi", type=int, default=220)
    return parser.parse_args()


def snapshot_index(all_snapshots: np.ndarray, requested: int) -> int:
    found = np.flatnonzero(all_snapshots == requested)
    if found.size != 1:
        raise ValueError(
            f"Snapshot {requested} is not uniquely present in {all_snapshots.tolist()}"
        )
    return int(found[0])


def load_qsl_slice(qsl_dir: Path, snapshot: int) -> dict[str, object]:
    candidates = sorted(
        qsl_dir.glob(f"fixed_y_qsl_t{snapshot:04d}_yp0p000000_delta*.npz")
    )
    if len(candidates) != 1:
        raise FileNotFoundError(
            f"Expected one exact-y=0 QSL file for snapshot {snapshot}, found "
            f"{len(candidates)} in {qsl_dir}"
        )
    path = candidates[0]
    with np.load(path) as archive:
        logq = np.asarray(archive["logq"], dtype=float)
        x_mm = np.asarray(archive["x_mm"], dtype=float)
        height_mm = np.asarray(archive["height_mm"], dtype=float)
        fixed_y_mm = float(np.asarray(archive["fixed_y_Mm"]))

    if logq.shape != x_mm.shape or logq.shape != height_mm.shape:
        raise ValueError(
            f"QSL coordinate/field shape mismatch in {path}: "
            f"logq={logq.shape}, x={x_mm.shape}, height={height_mm.shape}"
        )
    if not np.isclose(fixed_y_mm, 0.0, atol=1.0e-12):
        raise ValueError(f"Expected exact y=0 QSL data, got y={fixed_y_mm} Mm in {path}")

    return {
        "logq": logq,
        "x_mm": x_mm,
        "height_mm": height_mm,
        "fixed_y_mm": fixed_y_mm,
        "source": str(path.resolve()),
    }


def find_windowed_peak(
    values: np.ndarray,
    x_mm: np.ndarray,
    height_mm: np.ndarray,
    window: tuple[float, float, float, float],
) -> dict[str, float | int]:
    xmin, xmax, hmin, hmax = window
    inside = (
        (x_mm[None, :] >= xmin)
        & (x_mm[None, :] <= xmax)
        & (height_mm[:, None] >= hmin)
        & (height_mm[:, None] <= hmax)
        & np.isfinite(values)
    )
    if not np.any(inside):
        raise ValueError(f"No finite J/B values in peak-search window {window}")
    selected = np.where(inside, values, -np.inf)
    iz, ix = np.unravel_index(int(np.argmax(selected)), selected.shape)
    return {
        "x_Mm": float(x_mm[ix]),
        "height_Mm": float(height_mm[iz]),
        "j_over_b_per_Mm": float(values[iz, ix]),
        "x_index": int(ix),
        "height_index": int(iz),
    }


def main() -> None:
    args = parse_args()
    args.output_stem.parent.mkdir(parents=True, exist_ok=True)

    mpl_cache = Path("/private/tmp/matplotlib-j-over-b-comparison")
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(mpl_cache))

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.lines import Line2D

    with np.load(args.stack) as stack:
        all_snapshots = np.asarray(stack["snapshots"], dtype=int)
        all_times = np.asarray(stack["time_minutes"], dtype=float)
        x_mm = np.asarray(stack["x_Mm"], dtype=float)
        height_mm = np.asarray(stack["height_Mm"], dtype=float)
        ratio = np.asarray(stack["j_over_b_per_Mm"], dtype=float)

    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))
    visual = metadata["visualization"]
    limits = tuple(float(value) for value in visual["view_limits_Mm"])
    vmin = float(visual["global_vmin_per_Mm"])
    vmax = float(visual["global_vmax_per_Mm"])
    colormap_name = str(visual["colormap"])

    if visual["scale"] != "logarithmic":
        raise ValueError("Expected a logarithmic color scale in the source metadata")
    if colormap_name != "rainbow":
        raise ValueError(f"Expected the continuous rainbow colormap, got {colormap_name!r}")

    selected_indices = [snapshot_index(all_snapshots, value) for value in args.snapshots]
    selected_times = [float(all_times[index]) for index in selected_indices]
    selected_ratio = [ratio[index] for index in selected_indices]
    selected_qsl = [load_qsl_slice(args.qsl_dir, value) for value in args.snapshots]

    contour_levels = tuple(float(level) for level in args.contour_levels)
    if len(contour_levels) > len(CONTOUR_COLORS):
        raise ValueError(
            f"At most {len(CONTOUR_COLORS)} contour levels are supported; "
            f"got {len(contour_levels)}"
        )
    contour_colors = CONTOUR_COLORS[: len(contour_levels)]
    contour_linestyles = CONTOUR_LINESTYLES[: len(contour_levels)]
    peak_window = tuple(float(value) for value in args.peak_window)
    if not (peak_window[0] < peak_window[1] and peak_window[2] < peak_window[3]):
        raise ValueError(f"Invalid peak-search window: {peak_window}")
    selected_peaks = [
        find_windowed_peak(values, x_mm, height_mm, peak_window)
        for values in selected_ratio
    ]

    cmap = plt.get_cmap(colormap_name).copy()
    cmap.set_bad("white")
    norm = LogNorm(vmin=vmin, vmax=vmax, clip=True)

    figure, axes = plt.subplots(
        1,
        2,
        figsize=(13.2, 7.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    image = None
    contour_path_counts: dict[str, dict[str, int]] = {}
    for panel, (axis, values, qsl, peak, snapshot, time_minutes) in enumerate(
        zip(
            axes,
            selected_ratio,
            selected_qsl,
            selected_peaks,
            args.snapshots,
            selected_times,
            strict=True,
        )
    ):
        image = axis.pcolormesh(
            x_mm,
            height_mm,
            values,
            shading="nearest",
            cmap=cmap,
            norm=norm,
            rasterized=True,
        )
        panel_counts: dict[str, int] = {}
        for level, color, linestyle in zip(
            contour_levels, contour_colors, contour_linestyles, strict=True
        ):
            # A broader white duplicate acts as a halo so each contour remains
            # visible over every hue in the rainbow background.
            axis.contour(
                qsl["x_mm"],
                qsl["height_mm"],
                qsl["logq"],
                levels=[level],
                colors=["white"],
                linestyles=[linestyle],
                linewidths=[3.4],
                zorder=4,
            )
            contours = axis.contour(
                qsl["x_mm"],
                qsl["height_mm"],
                qsl["logq"],
                levels=[level],
                colors=[color],
                linestyles=[linestyle],
                linewidths=[1.8],
                zorder=5,
            )
            panel_counts[f"{level:g}"] = len(contours.allsegs[0])
        contour_path_counts[str(snapshot)] = panel_counts
        peak_x = float(peak["x_Mm"])
        peak_height = float(peak["height_Mm"])
        axis.scatter(
            [peak_x],
            [peak_height],
            marker="*",
            s=260,
            facecolor="#F0E442",
            edgecolor="#101010",
            linewidth=1.3,
            zorder=9,
        )
        axis.annotate(
            f"J/B peak\n({peak_x:.1f}, {peak_height:.1f}) Mm",
            xy=(peak_x, peak_height),
            xytext=(28, 30),
            textcoords="offset points",
            fontsize=10.5,
            ha="left",
            va="bottom",
            color="#101010",
            bbox={
                "boxstyle": "round,pad=0.25",
                "facecolor": "white",
                "edgecolor": "#303030",
                "alpha": 0.9,
            },
            arrowprops={
                "arrowstyle": "->",
                "color": "#101010",
                "linewidth": 1.2,
            },
            zorder=10,
        )
        axis.set_xlim(limits[0], limits[1])
        axis.set_ylim(limits[2], limits[3])
        axis.set_aspect("equal", adjustable="box")
        axis.set_xlabel("x [Mm]", fontsize=15)
        axis.set_title(
            f"({chr(97 + panel)})  t = {time_minutes:.1f} min",
            fontsize=16,
            pad=8,
        )
        axis.tick_params(axis="both", which="major", labelsize=13)

    axes[0].set_ylabel("Height [Mm]", fontsize=15)
    figure.suptitle(
        r"$y=0$: $|\nabla\times\mathbf{B}|/|\mathbf{B}|$",
        fontsize=17,
    )
    legend_handles = [
        Line2D(
            [0],
            [0],
            color=color,
            linestyle=linestyle,
            linewidth=2.4,
            label=rf"$\log_{{10}}Q={level:g}$",
        )
        for level, color, linestyle in zip(
            contour_levels, contour_colors, contour_linestyles, strict=True
        )
    ]
    figure.legend(
        handles=legend_handles,
        loc="outside lower center",
        ncol=len(legend_handles),
        frameon=False,
        fontsize=12,
        handlelength=3.2,
        columnspacing=2.0,
    )
    assert image is not None
    colorbar = figure.colorbar(image, ax=axes, fraction=0.035, pad=0.025)
    colorbar.set_label(
        r"$|\nabla\times\mathbf{B}|/|\mathbf{B}|$ [Mm$^{-1}$]",
        fontsize=15,
    )
    colorbar.ax.tick_params(labelsize=12)

    png_path = args.output_stem.with_suffix(".png")
    pdf_path = args.output_stem.with_suffix(".pdf")
    figure.savefig(png_path, dpi=args.dpi, facecolor="white")
    figure.savefig(pdf_path, dpi=args.dpi, facecolor="white")
    plt.close(figure)

    sidecar = {
        "quantity": "|curl(B)|/|B|",
        "plane": "exact y=0",
        "snapshots": list(args.snapshots),
        "time_minutes": selected_times,
        "view_limits_Mm": list(limits),
        "scale": "logarithmic",
        "colormap": colormap_name,
        "color_limits_per_Mm": [vmin, vmax],
        "source_stack": str(args.stack),
        "source_metadata": str(args.metadata),
        "qsl_contours": {
            "levels_log10Q": list(contour_levels),
            "equivalent_Q": [10.0**level for level in contour_levels],
            "colors": dict(
                zip((f"{level:g}" for level in contour_levels), contour_colors)
            ),
            "linestyles": dict(
                zip((f"{level:g}" for level in contour_levels), contour_linestyles)
            ),
            "white_halo_linewidth": 3.4,
            "colored_linewidth": 1.8,
            "path_counts_by_snapshot": contour_path_counts,
            "source_files": [qsl["source"] for qsl in selected_qsl],
        },
        "current_sheet_peak_proxy": {
            "definition": (
                "maximum |curl(B)|/|B| on the native J/B grid inside the "
                "specified x-height window"
            ),
            "search_window_Mm": {
                "x": [peak_window[0], peak_window[1]],
                "height": [peak_window[2], peak_window[3]],
            },
            "grid_half_cell_uncertainty_Mm": {
                "x": float(np.median(np.diff(x_mm)) / 2.0),
                "height": float(np.median(np.diff(height_mm)) / 2.0),
            },
            "by_snapshot": {
                str(snapshot): peak
                for snapshot, peak in zip(args.snapshots, selected_peaks, strict=True)
            },
        },
        "outputs": {"png": str(png_path), "pdf": str(pdf_path)},
    }
    json_path = args.output_stem.with_suffix(".json")
    json_path.write_text(json.dumps(sidecar, indent=2) + "\n", encoding="utf-8")

    print(f"PNG: {png_path}")
    print(f"PDF: {pdf_path}")
    print(f"Metadata: {json_path}")


if __name__ == "__main__":
    main()

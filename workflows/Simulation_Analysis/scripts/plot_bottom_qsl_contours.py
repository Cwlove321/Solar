#!/usr/bin/env python3
"""Plot bottom-boundary QSL maps with colored log10(Q) contours."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


PROJECT_CANDIDATES = (
    Path("/Users/xiaoxingyun/Desktop/jet_mhd"),
    Path("/Users/xiaoxingyun/Desktop/w2_jet_mhd"),
)
PROJECT_DIR = next((path for path in PROJECT_CANDIDATES if path.is_dir()), PROJECT_CANDIDATES[0])
QSL_DIR = PROJECT_DIR / "qfactor/bottom_fastqsl2"
OUTPUT_DIR = QSL_DIR / "colored_contours"
SNAPSHOTS = (10, 95)
TIME_MINUTES_PER_SNAPSHOT = 0.71558
EXTENT_MM = (-40.0, 80.0, -60.0, 60.0)
LOGQ_LEVELS = np.arange(1.0, 7.0)

# Bright, colorblind-friendly colors selected to remain visible on a grayscale
# QSL background.  Line styles provide a second, non-color distinction.
COLORS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")
LINESTYLES = ("solid", "dashed", "dashdot", "dotted", (0, (5, 1)), (0, (3, 1, 1, 1)))


def read_logq(path: Path) -> tuple[np.ndarray, dict[str, object]]:
    with np.load(path) as archive:
        q = np.asarray(archive["q"], dtype=float)
        xreg = np.asarray(archive["xreg"], dtype=float)
        yreg = np.asarray(archive["yreg"], dtype=float)
        delta = float(archive["delta"])
    logq = np.full(q.shape, np.nan)
    valid = np.isfinite(q) & (q > 0.0)
    logq[valid] = np.log10(q[valid])
    metadata = {
        "shape": list(q.shape),
        "finite_count": int(valid.sum()),
        "nan_count": int(np.isnan(q).sum()),
        "q_min": float(np.nanmin(q)),
        "q_max": float(np.nanmax(q)),
        "log10_q_min": float(np.nanmin(logq)),
        "log10_q_max": float(np.nanmax(logq)),
        "xreg": xreg.tolist(),
        "yreg": yreg.tolist(),
        "delta": delta,
    }
    return logq, metadata


def contour_path_counts(contour_set) -> list[int]:
    """Count disconnected contour paths at each requested level."""
    return [len(paths) for paths in contour_set.allsegs]


def draw_panel(
    axis,
    logq: np.ndarray,
    snapshot: int,
    *,
    show_ylabel: bool,
    background_max: float,
) -> tuple[object, list[int]]:
    ny, nx = logq.shape
    x = np.linspace(EXTENT_MM[0], EXTENT_MM[1], nx)
    y = np.linspace(EXTENT_MM[2], EXTENT_MM[3], ny)
    image = axis.imshow(
        logq,
        origin="lower",
        extent=EXTENT_MM,
        cmap="gray_r",
        vmin=0.3,
        vmax=background_max,
        interpolation="nearest",
        rasterized=True,
    )
    contours = axis.contour(
        x,
        y,
        logq,
        levels=LOGQ_LEVELS,
        colors=COLORS,
        linestyles=LINESTYLES,
        linewidths=1.45,
        antialiased=True,
    )
    axis.set_title(
        f"t = {snapshot * TIME_MINUTES_PER_SNAPSHOT:.1f} min",
        fontsize=15,
        pad=9,
    )
    axis.set_xlabel("Solar-X [Mm]")
    if show_ylabel:
        axis.set_ylabel("Solar-Y [Mm]")
    axis.set_xlim(EXTENT_MM[:2])
    axis.set_ylim(EXTENT_MM[2:])
    axis.set_aspect("equal")
    axis.tick_params(direction="in", top=True, right=True)
    return image, contour_path_counts(contours)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qsl-dir", type=Path, default=QSL_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--snapshots", type=int, nargs="+", default=list(SNAPSHOTS))
    parser.add_argument("--background-max", type=float, default=6.5)
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    datasets = []
    for snapshot in args.snapshots:
        source = args.qsl_dir / f"bottom_fastqsl2_t{snapshot:04d}_delta0.25.npz"
        if not source.is_file():
            raise FileNotFoundError(source)
        logq, stats = read_logq(source)
        datasets.append({"snapshot": snapshot, "source": source, "logq": logq, "stats": stats})

    figure, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(18.0, 8.0),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    image = None
    for index, (axis, dataset) in enumerate(zip(axes, datasets)):
        image, counts = draw_panel(
            axis,
            dataset["logq"],
            dataset["snapshot"],
            show_ylabel=index == 0,
            background_max=args.background_max,
        )
        dataset["contour_path_counts"] = counts

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=color,
            linestyle=style,
            linewidth=2.2,
            label=rf"$\log_{{10}}Q={level:.0f}$  ($Q=10^{{{level:.0f}}}$)",
        )
        for level, color, style in zip(LOGQ_LEVELS, COLORS, LINESTYLES)
    ]
    figure.legend(
        handles=legend_handles,
        loc="outside lower center",
        ncol=6,
        frameon=False,
        handlelength=3.0,
        columnspacing=1.35,
        fontsize=10,
    )
    colorbar = figure.colorbar(image, ax=axes, shrink=0.91, pad=0.02)
    colorbar.set_label(r"Bottom-boundary $\log_{10}Q$")
    figure.suptitle(
        r"Bottom-boundary QSL with colored $\log_{10}Q$ contours",
        fontsize=18,
    )

    output_base = args.output_dir / "bottom_qsl_logq_colored_contours_t0010_t0095"
    figure.savefig(output_base.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight")
    figure.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)

    metadata = {
        "interpretation": "Requested Q=1..6 contours are plotted as log10(Q)=1..6",
        "contour_levels_log10_q": LOGQ_LEVELS.tolist(),
        "equivalent_q_values": [float(10**level) for level in LOGQ_LEVELS],
        "colors": dict(zip((str(int(v)) for v in LOGQ_LEVELS), COLORS)),
        "linestyles": [str(style) for style in LINESTYLES],
        "background": {
            "quantity": "log10(Q)",
            "colormap": "gray_r",
            "limits": [0.3, args.background_max],
        },
        "extent_Mm": list(EXTENT_MM),
        "outputs": {
            "png": str(output_base.with_suffix(".png").resolve()),
            "pdf": str(output_base.with_suffix(".pdf").resolve()),
        },
        "snapshots": [
            {
                "snapshot": dataset["snapshot"],
                "time_minutes": dataset["snapshot"] * TIME_MINUTES_PER_SNAPSHOT,
                "source": str(dataset["source"].resolve()),
                "statistics": dataset["stats"],
                "contour_path_counts_by_level": dict(
                    zip((str(int(v)) for v in LOGQ_LEVELS), dataset["contour_path_counts"])
                ),
            }
            for dataset in datasets
        ],
    }
    metadata_path = output_base.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(output_base.with_suffix(".png").resolve())
    print(output_base.with_suffix(".pdf").resolve())
    print(metadata_path.resolve())


if __name__ == "__main__":
    main()

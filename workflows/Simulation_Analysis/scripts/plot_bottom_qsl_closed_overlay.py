#!/usr/bin/env python3
"""Overlay FastQSL2 bottom-to-bottom closed-field regions on bottom QSL maps."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch


QSL_DIR = Path("/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/bottom_fastqsl2")
OUTPUT_DIR = Path(
    "/Users/xiaoxingyun/.codex/visualizations/2026/08/21/"
    "01a02330-d5b3-7762-968c-7c51d9bfcca7/bottom_qsl_closed_overlay"
)
SNAPSHOTS = (10, 95)
TIME_MINUTES_PER_SNAPSHOT = 0.71558
EXTENT_MM = (-40.0, 80.0, -60.0, 60.0)
LOGQ_LEVELS = np.arange(1.0, 7.0)
Q_COLORS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")
Q_LINESTYLES = ("solid", "dashed", "dashdot", "dotted", (0, (5, 1)), (0, (3, 1, 1, 1)))
CLOSED_FILL = "#8F77B5"
CLOSED_EDGE = "#49306B"


def load_bottom_result(path: Path) -> dict[str, np.ndarray | float | list[int]]:
    with np.load(path) as archive:
        required = {"q", "rboundary", "rFs", "rFe", "xreg", "yreg", "delta"}
        missing = required.difference(archive.files)
        if missing:
            raise KeyError(f"{path} is missing FastQSL2 outputs: {sorted(missing)}")
        q = np.asarray(archive["q"], dtype=float)
        rboundary = np.asarray(archive["rboundary"], dtype=np.int16)
        rfs = np.asarray(archive["rFs"], dtype=float)
        rfe = np.asarray(archive["rFe"], dtype=float)
        xreg = np.asarray(archive["xreg"], dtype=float)
        yreg = np.asarray(archive["yreg"], dtype=float)
        delta = float(archive["delta"])

    if q.shape != rboundary.shape or rfs.shape != q.shape + (3,) or rfe.shape != q.shape + (3,):
        raise ValueError(f"Inconsistent FastQSL2 array shapes in {path}")

    logq = np.full(q.shape, np.nan)
    valid_q = np.isfinite(q) & (q > 0.0)
    logq[valid_q] = np.log10(q[valid_q])

    # FastQSL2 encodes the two terminal boundaries as 10*rbs+rbe.
    # Boundary code 1 is the photospheric/bottom boundary, so 11 means
    # that both integrated directions terminate at the bottom boundary.
    closed = rboundary == 11
    start_bottom = np.isclose(rfs[..., 2], 0.0, atol=1.0e-3)
    end_bottom = np.isclose(rfe[..., 2], 0.0, atol=1.0e-3)
    endpoint_bottom_bottom = start_bottom & end_bottom
    # Every code-11 line must have z=0 at both terminals. The converse is not
    # sufficient: a side/bottom edge can also have z=0 while rboundary
    # correctly records the side or corner rather than the bottom.
    invalid_closed = closed & ~endpoint_bottom_bottom
    if np.any(invalid_closed):
        count = int(np.count_nonzero(invalid_closed))
        raise ValueError(f"{count} code-11 lines do not have two z=0 terminals in {path}")

    return {
        "q": q,
        "logq": logq,
        "rboundary": rboundary,
        "rFs": rfs,
        "rFe": rfe,
        "closed": closed,
        "valid_q": valid_q,
        "xreg": xreg,
        "yreg": yreg,
        "delta": delta,
    }


def draw_panel(axis, data: dict, snapshot: int, *, show_ylabel: bool, background_max: float):
    logq = data["logq"]
    closed = data["closed"]
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
        zorder=1,
    )
    closed_layer = np.ma.masked_where(~closed, np.ones_like(closed, dtype=float))
    axis.imshow(
        closed_layer,
        origin="lower",
        extent=EXTENT_MM,
        cmap=ListedColormap([CLOSED_FILL]),
        vmin=0.0,
        vmax=1.0,
        interpolation="nearest",
        alpha=0.34,
        rasterized=True,
        zorder=2,
    )
    axis.contour(
        x,
        y,
        closed.astype(float),
        levels=[0.5],
        colors=[CLOSED_EDGE],
        linewidths=0.95,
        zorder=3,
    )
    q_contours = axis.contour(
        x,
        y,
        logq,
        levels=LOGQ_LEVELS,
        colors=Q_COLORS,
        linestyles=Q_LINESTYLES,
        linewidths=1.45,
        antialiased=True,
        zorder=4,
    )

    axis.set_title(f"t = {snapshot * TIME_MINUTES_PER_SNAPSHOT:.1f} min", fontsize=15, pad=9)
    axis.set_xlabel("Solar-X [Mm]")
    if show_ylabel:
        axis.set_ylabel("Solar-Y [Mm]")
    axis.set_xlim(EXTENT_MM[:2])
    axis.set_ylim(EXTENT_MM[2:])
    axis.set_aspect("equal")
    axis.tick_params(direction="in", top=True, right=True)
    return image, [len(paths) for paths in q_contours.allsegs]


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
        data = load_bottom_result(source)
        mask_path = args.output_dir / f"closed_bottom_to_bottom_t{snapshot:04d}.npz"
        np.savez_compressed(
            mask_path,
            closed=data["closed"].astype(np.uint8),
            rboundary=data["rboundary"].astype(np.int8),
            rFs=data["rFs"].astype(np.float32),
            rFe=data["rFe"].astype(np.float32),
            extent_Mm=np.asarray(EXTENT_MM),
            source=np.asarray(str(source.resolve())),
            definition=np.asarray("closed iff FastQSL2 rboundary == 11 (both terminals at bottom)"),
        )
        datasets.append({"snapshot": snapshot, "source": source, "mask": mask_path, "data": data})

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
            dataset["data"],
            dataset["snapshot"],
            show_ylabel=index == 0,
            background_max=args.background_max,
        )
        dataset["contour_counts"] = counts

    legend_handles = [
        Patch(
            facecolor=CLOSED_FILL,
            edgecolor=CLOSED_EDGE,
            alpha=0.55,
            label="Closed field: both terminals at bottom",
        )
    ]
    legend_handles.extend(
        Line2D(
            [0],
            [0],
            color=color,
            linestyle=style,
            linewidth=2.2,
            label=rf"$\log_{{10}}Q={level:.0f}$",
        )
        for level, color, style in zip(LOGQ_LEVELS, Q_COLORS, Q_LINESTYLES)
    )
    figure.legend(
        handles=legend_handles,
        loc="outside lower center",
        ncol=7,
        frameon=False,
        handlelength=3.0,
        columnspacing=1.25,
        fontsize=10,
    )
    colorbar = figure.colorbar(image, ax=axes, shrink=0.91, pad=0.02)
    colorbar.set_label(r"Bottom-boundary $\log_{10}Q$")
    figure.suptitle(
        r"Bottom-boundary QSL, colored contours, and bottom-to-bottom closed field",
        fontsize=18,
    )

    output_base = args.output_dir / "bottom_qsl_contours_closed_field_t0010_t0095"
    figure.savefig(output_base.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight")
    figure.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)

    total_pixels = int(np.prod(datasets[0]["data"]["closed"].shape))
    metadata = {
        "closed_field_definition": {
            "source": "FastQSL2 field-line integration performed during bottom-QSL calculation",
            "criterion": "rboundary == 11",
            "meaning": "both terminal points are on the photospheric/bottom boundary",
            "independent_check": "all rboundary=11 pixels have rFs_z=rFe_z=0 within 1e-3",
            "side_or_top_boundary_and_unfinished_lines": "excluded from closed mask",
        },
        "q_contours": {
            "levels_log10_q": LOGQ_LEVELS.tolist(),
            "equivalent_q": [float(10**level) for level in LOGQ_LEVELS],
            "colors": dict(zip((str(int(v)) for v in LOGQ_LEVELS), Q_COLORS)),
        },
        "closed_overlay": {"fill": CLOSED_FILL, "edge": CLOSED_EDGE, "alpha": 0.34},
        "extent_Mm": list(EXTENT_MM),
        "outputs": {
            "png": str(output_base.with_suffix(".png").resolve()),
            "pdf": str(output_base.with_suffix(".pdf").resolve()),
        },
        "snapshots": [],
    }
    for dataset in datasets:
        data = dataset["data"]
        rb_values, rb_counts = np.unique(data["rboundary"], return_counts=True)
        closed_count = int(np.count_nonzero(data["closed"]))
        metadata["snapshots"].append(
            {
                "snapshot": dataset["snapshot"],
                "time_minutes": dataset["snapshot"] * TIME_MINUTES_PER_SNAPSHOT,
                "source": str(dataset["source"].resolve()),
                "mask_npz": str(dataset["mask"].resolve()),
                "shape": list(data["closed"].shape),
                "closed_pixel_count": closed_count,
                "closed_fraction": closed_count / total_pixels,
                "rboundary_counts": {
                    str(int(value)): int(count) for value, count in zip(rb_values, rb_counts)
                },
                "q_contour_path_counts": dict(
                    zip((str(int(v)) for v in LOGQ_LEVELS), dataset["contour_counts"])
                ),
            }
        )

    metadata_path = output_base.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(output_base.with_suffix(".png").resolve())
    print(output_base.with_suffix(".pdf").resolve())
    print(metadata_path.resolve())


if __name__ == "__main__":
    main()

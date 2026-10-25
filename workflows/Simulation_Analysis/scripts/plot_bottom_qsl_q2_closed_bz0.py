#!/usr/bin/env python3
"""Plot log10(Q)=2, Bz=0, signed closed regions, and open-field annotations."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import distance_transform_edt


QSL_DIR = Path("/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/bottom_fastqsl2")
MASK_DIR = QSL_DIR / "colored_contours"
MAGNETIC_DIR = Path("/Volumes/Seagate/jet_mhd/out_path1/nlfff")
OUTPUT_DIR = Path(
    "/Users/xiaoxingyun/.codex/visualizations/2026/08/21/"
    "01a02330-d5b3-7762-968c-7c51d9bfcca7/bottom_qsl_q2_closed_bz0"
)
SNAPSHOTS = (10, 95)
TIME_MINUTES_PER_SNAPSHOT = 0.71558
EXTENT_MM = (-40.0, 80.0, -60.0, 60.0)
Q2_LEVEL = 2.0
Q2_COLOR = "#E69F00"
BZ0_COLOR = "#D62728"
CLOSED_POSITIVE = "#CC79A7"
CLOSED_NEGATIVE = "#0072B2"
CLOSED_ALPHA = 0.30
CLOSED_EDGE = "#49306B"
MAGNETIC_UNIT_G = 2.0
PANEL_TITLE_FONTSIZE = 20
AXIS_LABEL_FONTSIZE = 18
TICK_FONTSIZE = 16


def split_closed_polarity(bz, x_b, y_b, closed, x_q, y_q):
    """Sample native cell-center Bz on the QSL grid; never change closed mask."""
    x_mesh, y_mesh = np.meshgrid(x_q, y_q)
    # Clamp only the outer half-cell; the displayed cropped domain is interior.
    points = np.column_stack((np.clip(y_mesh.ravel(), y_b[0], y_b[-1]),
                              np.clip(x_mesh.ravel(), x_b[0], x_b[-1])))
    bz_on_q = RegularGridInterpolator((y_b, x_b), bz)(points).reshape(closed.shape)
    if not np.all(np.isfinite(bz_on_q[closed])):
        raise ValueError("Non-finite Bz in the closed-field mask")
    positive = closed & (bz_on_q > 0)
    negative = closed & (bz_on_q < 0)
    neutral = closed & (bz_on_q == 0)
    assert np.array_equal(positive | negative | neutral, closed)
    assert not np.any(positive & negative)
    return bz_on_q, positive, negative, neutral


def label_polarities(axis, positive, negative, bz_on_q, x_q, y_q, view_extent, q_segments):
    """Place one sign inside each polarity, away from boundaries and Q lines."""
    dx, dy = x_q[1] - x_q[0], y_q[1] - y_q[0]
    visible = ((y_q[:, None] > view_extent[2]) & (y_q[:, None] < view_extent[3])
               & (x_q[None, :] > view_extent[0]) & (x_q[None, :] < view_extent[1]))
    q_edge = np.zeros(positive.shape, dtype=bool)
    for segment in q_segments:
        points = segment[np.all(np.isfinite(segment), axis=1)]
        ix = np.clip(np.rint((points[:, 0] - x_q[0]) / dx).astype(int), 0, len(x_q) - 1)
        iy = np.clip(np.rint((points[:, 1] - y_q[0]) / dy).astype(int), 0, len(y_q) - 1)
        q_edge[iy, ix] = True
    q_clearance = (distance_transform_edt(~q_edge, sampling=(dy, dx))
                   if q_edge.any() else np.full(positive.shape, np.inf))
    positions = []
    for polarity, mask, sign in (("positive", positive, "+"), ("negative", negative, "-")):
        interior = mask & visible
        if not interior.any():
            continue
        clearance = np.minimum(distance_transform_edt(interior, sampling=(dy, dx)), q_clearance)
        clearance[~interior] = -1
        iy, ix = np.unravel_index(np.argmax(clearance), clearance.shape)
        axis.text(x_q[ix], y_q[iy], f"${sign}$", fontsize=34, fontweight="bold",
                  color="#202020", ha="center", va="center", zorder=6, clip_on=True)
        positions.append({"polarity": polarity, "x_Mm": float(x_q[ix]),
                          "y_Mm": float(y_q[iy]), "bz_G": float(bz_on_q[iy, ix]),
                          "clearance_Mm": float(clearance[iy, ix])})
    return positions


def label_open_region(axis, rboundary, bz_on_q, x_q, y_q, view_extent):
    """Annotate the requested open-background locations after checking the data."""
    annotations = []
    for label, fx, fy in (("Open field", 0.95, 0.95), ("-", 0.87, 0.12)):
        x = view_extent[0] + fx * (view_extent[1] - view_extent[0])
        y = view_extent[2] + fy * (view_extent[3] - view_extent[2])
        ix, iy = int(np.argmin(abs(x_q - x))), int(np.argmin(abs(y_q - y)))
        code = int(rboundary[iy, ix])
        # Boundary 1 is bottom; 2..6 are side/top boundaries of the box.
        terminals = divmod(code, 10)
        if sorted(terminals)[0] != 1 or not 2 <= sorted(terminals)[1] <= 6:
            raise ValueError(f"Open-field label at {(x, y)} has unexpected boundary code {code}")
        if label == "-" and not bz_on_q[iy, ix] < 0:
            raise ValueError(f"Open negative label at {(x, y)} does not have Bz<0")
        axis.text(x, y, label if label == "Open field" else "$-$",
                  fontsize=20 if label == "Open field" else 34,
                  ha="right" if label == "Open field" else "center",
                  va="top" if label == "Open field" else "center",
                  color="#202020", clip_on=True, zorder=6)
        annotations.append({"label": label, "x_Mm": x, "y_Mm": y,
                            "rboundary": code, "bz_G": float(bz_on_q[iy, ix])})
    return annotations


def read_bottom_bz(path: Path) -> np.ndarray:
    """Read bottom Bz from a QSL binary stored as three (z,y,x) float64 arrays."""
    with path.open("rb") as stream:
        header = np.fromfile(stream, dtype=np.int32, count=3)
        if header.size != 3:
            raise ValueError(f"Missing magnetic-field header in {path}")
        nx, ny, nz = (int(value) for value in header)
        count = nx * ny * nz
        stream.seek(2 * count * np.dtype(np.float64).itemsize, 1)
        bz = np.fromfile(stream, dtype=np.float64, count=count)
    if bz.size != count:
        raise ValueError(f"Incomplete Bz array in {path}")
    return bz.reshape((nz, ny, nx))[0] * MAGNETIC_UNIT_G


def read_qsl_and_closed(qsl_path: Path, mask_path: Path) -> dict[str, np.ndarray]:
    with np.load(qsl_path) as archive:
        q = np.asarray(archive["q"], dtype=float)
        rboundary = np.asarray(archive["rboundary"], dtype=np.int16)
    with np.load(mask_path) as archive:
        closed = np.asarray(archive["closed"], dtype=bool)
    if q.shape != closed.shape or q.shape != rboundary.shape:
        raise ValueError(f"Shape mismatch among Q and closed-field arrays for {qsl_path}")
    if not np.array_equal(closed, rboundary == 11):
        raise ValueError(f"Closed mask is not identical to rboundary==11 for {qsl_path}")
    logq = np.full(q.shape, np.nan)
    valid = np.isfinite(q) & (q > 0.0)
    logq[valid] = np.log10(q[valid])
    return {"q": q, "logq": logq, "closed": closed, "rboundary": rboundary}


def draw_panel(
    axis, dataset: dict, *, show_ylabel: bool, background_max: float,
    background: str = "qsl", bz_vmax: float = 50.0,
    view_extent: tuple[float, float, float, float] = EXTENT_MM,
):
    logq = dataset["logq"]
    closed = dataset["closed"]
    bz = dataset["bz"]
    ny_q, nx_q = logq.shape
    ny_b, nx_b = bz.shape
    x_q = np.linspace(EXTENT_MM[0], EXTENT_MM[1], nx_q)
    y_q = np.linspace(EXTENT_MM[2], EXTENT_MM[3], ny_q)
    dx_b = (EXTENT_MM[1] - EXTENT_MM[0]) / nx_b
    dy_b = (EXTENT_MM[3] - EXTENT_MM[2]) / ny_b
    x_b = EXTENT_MM[0] + (np.arange(nx_b) + 0.5) * dx_b
    y_b = EXTENT_MM[2] + (np.arange(ny_b) + 0.5) * dy_b

    image = None
    axis.set_facecolor("white")
    if background != "white":
        image = axis.imshow(
            bz if background == "bz" else logq,
            origin="lower",
            extent=EXTENT_MM,
            cmap="gray" if background == "bz" else "gray_r",
            vmin=-bz_vmax if background == "bz" else 0.3,
            vmax=bz_vmax if background == "bz" else background_max,
            interpolation="nearest",
            rasterized=True,
            zorder=1,
        )
    bz_on_q, positive, negative, neutral = split_closed_polarity(
        bz, x_b, y_b, closed, x_q, y_q,
    )
    categories = np.zeros_like(bz_on_q)
    categories[positive] = 1
    categories[neutral] = 2
    closed_layer = np.ma.masked_where(~closed, categories)
    axis.pcolormesh(
        x_q, y_q, closed_layer,
        cmap=ListedColormap([CLOSED_NEGATIVE, CLOSED_POSITIVE, "#BBBBBB"]),
        vmin=-0.5,
        vmax=2.5,
        shading="nearest",
        alpha=CLOSED_ALPHA,
        rasterized=True,
        zorder=2,
    )
    closed_contour = axis.contour(
        x_q,
        y_q,
        closed.astype(float),
        levels=[0.5],
        colors=[CLOSED_EDGE],
        linewidths=0.95,
        zorder=3,
    )
    q2_contour = axis.contour(
        x_q,
        y_q,
        logq,
        levels=[Q2_LEVEL],
        colors=[Q2_COLOR],
        linestyles="solid",
        linewidths=1.8,
        antialiased=True,
        zorder=4,
    )
    bz0_contour = axis.contour(
        x_b, y_b, bz,
        levels=[0.0], colors=[BZ0_COLOR], linestyles="dashed",
        linewidths=1.35, antialiased=True, zorder=5,
    )
    polarity_labels = label_polarities(
        axis, positive, negative, bz_on_q, x_q, y_q, view_extent, q2_contour.allsegs[0],
    )
    open_annotations = label_open_region(
        axis, dataset["rboundary"], bz_on_q, x_q, y_q, view_extent,
    )

    snapshot = dataset["snapshot"]
    axis.set_title(f"t = {snapshot * TIME_MINUTES_PER_SNAPSHOT:.1f} min", fontsize=PANEL_TITLE_FONTSIZE, pad=9)
    axis.set_xlabel("Solar-X [Mm]", fontsize=AXIS_LABEL_FONTSIZE)
    if show_ylabel:
        axis.set_ylabel("Solar-Y [Mm]", fontsize=AXIS_LABEL_FONTSIZE)
    # Crop the view only: retain the original coordinates for every data layer.
    axis.set_xlim(view_extent[:2])
    axis.set_ylim(view_extent[2:])
    axis.set_aspect("equal")
    axis.tick_params(direction="in", top=True, right=True, labelsize=TICK_FONTSIZE)
    return {
        "image": image,
        "closed_path_count": len(closed_contour.allsegs[0]),
        "q2_path_count": len(q2_contour.allsegs[0]),
        "bz0_path_count": len(bz0_contour.allsegs[0]),
        "closed_positive_pixels": int(positive.sum()),
        "closed_negative_pixels": int(negative.sum()),
        "closed_zero_bz_pixels": int(neutral.sum()),
        "polarity_labels": polarity_labels,
        "open_annotations": open_annotations,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qsl-dir", type=Path, default=QSL_DIR)
    parser.add_argument("--mask-dir", type=Path, default=MASK_DIR)
    parser.add_argument("--magnetic-dir", type=Path, default=MAGNETIC_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--snapshots", type=int, nargs="+", default=list(SNAPSHOTS))
    parser.add_argument("--background-max", type=float, default=6.5)
    parser.add_argument("--background", choices=["qsl", "bz", "white"], default="qsl")
    parser.add_argument(
        "--bz-vmax", type=float, default=50.0,
        help="Symmetric Bz grayscale limit in G; identical for all panels.",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--crop-mm", type=float, default=0.0,
                        help="Inset each edge of the displayed spatial domain by this many Mm.")
    parser.add_argument("--legend-fontsize", type=float, default=20.0)
    args = parser.parse_args()
    if not np.isfinite(args.bz_vmax) or args.bz_vmax <= 0:
        parser.error("--bz-vmax must be positive and finite")
    max_crop = min(EXTENT_MM[1] - EXTENT_MM[0], EXTENT_MM[3] - EXTENT_MM[2]) / 2
    if not np.isfinite(args.crop_mm) or not 0 <= args.crop_mm < max_crop:
        parser.error(f"--crop-mm must be finite and in [0, {max_crop})")
    if not np.isfinite(args.legend_fontsize) or args.legend_fontsize <= 0:
        parser.error("--legend-fontsize must be positive and finite")
    view_extent = (
        EXTENT_MM[0] + args.crop_mm, EXTENT_MM[1] - args.crop_mm,
        EXTENT_MM[2] + args.crop_mm, EXTENT_MM[3] - args.crop_mm,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)

    datasets = []
    for snapshot in args.snapshots:
        qsl_path = args.qsl_dir / f"bottom_fastqsl2_t{snapshot:04d}_delta0.25.npz"
        mask_path = args.mask_dir / f"closed_bottom_to_bottom_t{snapshot:04d}.npz"
        magnetic_path = args.magnetic_dir / f"nlfff_{snapshot:04d}.binary"
        for path in (qsl_path, mask_path, magnetic_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        data = read_qsl_and_closed(qsl_path, mask_path)
        data.update(
            {
                "snapshot": snapshot,
                "bz": read_bottom_bz(magnetic_path),
                "qsl_path": qsl_path,
                "mask_path": mask_path,
                "magnetic_path": magnetic_path,
            }
        )
        datasets.append(data)

    figure, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(14.5, 9.0) if args.background == "white" else (18.0, 8.0),
        sharex=True,
        sharey=True,
        layout="none" if args.background == "white" else "constrained",
    )
    axes = np.atleast_1d(axes)
    if args.background == "white":
        # Explicit bands keep enlarged title/axis labels/legend from colliding.
        figure.subplots_adjust(left=0.08, right=0.985, bottom=0.15, top=0.885, wspace=0.07)
    figure.set_facecolor("white")
    panel_results = []
    for index, (axis, dataset) in enumerate(zip(axes, datasets)):
        panel_results.append(
            draw_panel(
                axis, dataset, show_ylabel=index == 0,
                background_max=args.background_max,
                background=args.background, bz_vmax=args.bz_vmax,
                view_extent=view_extent,
            )
        )
        axis.tick_params(axis="y", labelleft=index == 0, labelright=False)

    legend_handles = [
        Patch(
            facecolor=CLOSED_POSITIVE,
            edgecolor=CLOSED_EDGE,
            alpha=CLOSED_ALPHA,
            label=r"Closed field ($+$)",
        ),
        Patch(
            facecolor=CLOSED_NEGATIVE,
            edgecolor=CLOSED_EDGE,
            alpha=CLOSED_ALPHA,
            label=r"Closed field ($-$)",
        ),
        Line2D([0], [0], color=Q2_COLOR, linewidth=2.3, label=r"$\log_{10}Q=2$"),
        Line2D([0], [0], color=BZ0_COLOR, linestyle="dashed", linewidth=2.0,
               label=r"$B_z=0$"),
    ]
    legend = figure.legend(
        handles=legend_handles,
        loc="lower center" if args.background == "white" else "outside lower center",
        bbox_to_anchor=(0.5, 0.005) if args.background == "white" else None,
        ncol=4,
        frameon=False,
        handlelength=2.1,
        columnspacing=1.2,
        fontsize=args.legend_fontsize,
    )
    if args.background != "white":
        colorbar = figure.colorbar(panel_results[0]["image"], ax=axes, shrink=0.91, pad=0.02)
        colorbar.set_label(
            r"Bottom-boundary $B_z$ [G]" if args.background == "bz"
            else r"Bottom-boundary $\log_{10}Q$"
        )
    if args.background == "white":
        figure.canvas.draw()
        renderer = figure.canvas.get_renderer()
        for axis in axes:
            bounds = axis.get_tightbbox(renderer)
            assert not bounds.overlaps(legend.get_window_extent(renderer)), "Axis/legend overlap"

    snapshot_tag = "_".join(f"t{snapshot:04d}" for snapshot in args.snapshots)
    output_base = args.output_dir / f"bottom_{args.background}_q2_closed_field_bz0_{snapshot_tag}"
    figure.savefig(output_base.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight", facecolor="white", transparent=False)
    figure.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight", facecolor="white", transparent=False)
    plt.close(figure)

    metadata = {
        "layout": {
            "share_y": True,
            "y_tick_labels": "leftmost panel only",
            "x_range_Mm": list(view_extent[:2]),
            "y_range_Mm": list(view_extent[2:]),
            "crop_each_edge_Mm": args.crop_mm,
            "legend_fontsize_pt": args.legend_fontsize,
            "title": None,
            "panel_title_fontsize_pt": PANEL_TITLE_FONTSIZE,
            "axis_label_fontsize_pt": AXIS_LABEL_FONTSIZE,
            "tick_fontsize_pt": TICK_FONTSIZE,
            "closed_legend_labels": ["Closed field (+)", "Closed field (-)"],
            "q_legend_label": "log10(Q)=2",
            "engine": "manual shared-y with reserved text bands" if args.background == "white" else "constrained",
        },
        "layers_bottom_to_top": [
            "pure white (no scalar background)" if args.background == "white" else
            "grayscale Bz [G]" if args.background == "bz" else "grayscale log10(Q)",
            "closed-field fill: pink Bz>0, blue Bz<0; unchanged purple boundary",
            "solid orange log10(Q)=2 contour",
            "dashed red Bz=0 contour",
            "positive and negative polarity signs inside closed regions",
            "Open field text upper right and negative sign lower right",
        ],
        "background": {
            "quantity": "Bz" if args.background == "bz" else "log10(Q)",
            "colormap": "gray" if args.background == "bz" else "gray_r",
            "limits": [-args.bz_vmax, args.bz_vmax] if args.background == "bz"
            else [0.3, args.background_max],
            "unit": "G" if args.background == "bz" else "dimensionless",
        },
        "q_contour": {"log10_q": Q2_LEVEL, "q": 10**Q2_LEVEL, "color": Q2_COLOR},
        "bz_contour": {"shown": True, "value_G": 0.0, "color": BZ0_COLOR,
                       "linestyle": "dashed"},
        "closed_definition": "FastQSL2 rboundary==11, both field-line terminals at bottom",
        "polarity_definition": {
            "quantity": "bottom-boundary Bz",
            "sampling": "bilinear from native Bz cell centers to QSL grid; clamp outer half-cell",
            "positive": {"condition": "Bz>0", "color": CLOSED_POSITIVE},
            "negative": {"condition": "Bz<0", "color": CLOSED_NEGATIVE},
            "exact_zero": "neutral gray if present",
            "alpha": CLOSED_ALPHA,
            "label_placement": "maximize clearance from polarity/closed edges and Q contour within viewport",
        },
        "extent_Mm": list(EXTENT_MM),
        "view_extent_Mm": list(view_extent),
        "statistics_scope": "full source domain, not the cropped viewport",
        "outputs": {
            "png": str(output_base.with_suffix(".png").resolve()),
            "pdf": str(output_base.with_suffix(".pdf").resolve()),
        },
        "snapshots": [],
    }
    if args.background == "white":
        metadata["background"] = {
            "quantity": None, "color": "#FFFFFF", "colorbar": False,
            "colormap": None, "limits": None, "unit": None,
        }
    for dataset, panel in zip(datasets, panel_results):
        metadata["snapshots"].append(
            {
                "snapshot": dataset["snapshot"],
                "time_minutes": dataset["snapshot"] * TIME_MINUTES_PER_SNAPSHOT,
                "qsl_source": str(dataset["qsl_path"].resolve()),
                "closed_mask_source": str(dataset["mask_path"].resolve()),
                "magnetic_source": str(dataset["magnetic_path"].resolve()),
                "bz_shape": list(dataset["bz"].shape),
                "bz_range_G": [float(np.min(dataset["bz"])), float(np.max(dataset["bz"]))],
                "closed_fraction": float(np.mean(dataset["closed"])),
                "closed_boundary_path_count": panel["closed_path_count"],
                "log10_q_2_path_count": panel["q2_path_count"],
                "bz_0_path_count": panel["bz0_path_count"],
                "closed_positive_pixels": panel["closed_positive_pixels"],
                "closed_negative_pixels": panel["closed_negative_pixels"],
                "closed_zero_bz_pixels": panel["closed_zero_bz_pixels"],
                "polarity_labels": panel["polarity_labels"],
                "open_annotations": panel["open_annotations"],
            }
        )
    metadata_path = output_base.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(output_base.with_suffix(".png").resolve())
    print(output_base.with_suffix(".pdf").resolve())
    print(metadata_path.resolve())


if __name__ == "__main__":
    main()

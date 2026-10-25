#!/usr/bin/env python3
"""Plot a candidate overlying-rope flux-integration footprint on closed/Q maps.

The input ``axis_features_factor2_tXXXX.npz`` files contain field-line features
measured in a local coordinate system tied to the three-dimensional magnetic
axis/current channel.  This script selects positive-polarity footpoints whose
field lines

* cross both lateral sides of the current channel,
* pass above the central part of the magnetic axis,
* cover at least 25 per cent of that central axis, and
* do not wind once around the axis (``abs(winding_number) < 1``).

The resulting footprint is a physically motivated *candidate* for integrating
the flux that runs across and over the erupting current channel.  It is shown
separately from the full bottom-to-bottom closed-field region and is not defined
as the interior of a chosen Q contour.

Example
-------
/opt/anaconda3/bin/python \
  /Users/xiaoxingyun/Desktop/solar_code/workflows/Simulation_Analysis/scripts/\
plot_candidate_overlying_flux_on_closed_q.py
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from scipy.ndimage import distance_transform_edt


QSL_DIR = Path("/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/bottom_fastqsl2")
FEATURE_DIR = QSL_DIR / "overlying_flux_candidate"
OUTPUT_DIR = FEATURE_DIR
SNAPSHOTS = (10, 95)

TIME_MINUTES_PER_SNAPSHOT = 0.71558
PHYSICAL_EXTENT_MM = (-40.0, 80.0, -60.0, 60.0)
MAIN_LIMITS_MM = (-20.0, 60.0, -40.0, 40.0)
ZOOM_LIMITS_MM = (-10.0, 15.0, -12.0, 15.0)

LOGQ_LEVELS = np.arange(1.0, 7.0)
Q_COLORS = ("#0072B2", "#E69F00", "#009E73", "#CC79A7", "#D55E00", "#56B4E9")
Q_LINESTYLES = ("solid", "dashed", "dashdot", "dotted", (0, (5, 1)), (0, (3, 1, 1, 1)))
CLOSED_FILL = "#8F77B5"
CLOSED_EDGE = "#49306B"
CANDIDATE_FILL = "#F0E442"
CANDIDATE_EDGE = "#615600"
PAPER_CANDIDATE_FILL = "#00BFC4"
PAPER_CANDIDATE_EDGE = "#006B70"
PAPER_POSITIVE_FILL = "#CC79A7"
PAPER_NEGATIVE_FILL = "#0072B2"
PAPER_BZ0_COLOR = "#D62728"

# The magnetic field stored in the simulation is normalized by 2 G.  The
# factor-2 seed grid has a spacing of 0.5 native cells, while the 128-cell-wide
# bottom plane spans 120 Mm.
B_UNIT_G = 2.0
NATIVE_CELL_MM = 120.0 / 128.0
SEED_STEP_NATIVE = 0.5

# Parameters determined from the initial current-channel geometry at each
# snapshot.  Distances are native-grid units in the local axis coordinates.
# ``current_density_fraction`` records how the supplied feature files were made.
SELECTION_PARAMETERS = {
    10: {
        "a_minus": 5.0,
        "a_plus": 7.0,
        "height": 5.0,
        "current_density_fraction": 0.25,
        "axis_coverage_min": 0.25,
        "edge_tolerance": 0.5,
        "winding_abs_max": 1.0,
    },
    95: {
        "a_minus": 3.0,
        "a_plus": 7.0,
        "height": 5.0,
        "current_density_fraction": 0.25,
        "axis_coverage_min": 0.25,
        "edge_tolerance": 0.5,
        "winding_abs_max": 1.0,
    },
}

FEATURE_COLUMNS = (
    "eta_max",
    "eta_min",
    "xi_min",
    "xi_max",
    "eta_center_max",
    "min_axis_distance",
    "central_axis_coverage",
    "winding_number",
    "median_axis_distance",
)


def load_bottom_result(path: Path) -> dict[str, np.ndarray | float]:
    """Load and validate the FastQSL2 bottom-plane product."""

    with np.load(path) as archive:
        required = {"q", "rboundary", "rFs", "rFe", "xreg", "yreg", "delta", "sign2d"}
        missing = required.difference(archive.files)
        if missing:
            raise KeyError(f"{path} is missing FastQSL2 outputs: {sorted(missing)}")
        q = np.asarray(archive["q"], dtype=float)
        rboundary = np.asarray(archive["rboundary"], dtype=np.int16)
        rfs = np.asarray(archive["rFs"], dtype=float)
        rfe = np.asarray(archive["rFe"], dtype=float)
        sign2d = np.asarray(archive["sign2d"], dtype=np.int8)
        xreg = np.asarray(archive["xreg"], dtype=float)
        yreg = np.asarray(archive["yreg"], dtype=float)
        delta = float(archive["delta"])

    if (
        q.shape != rboundary.shape
        or sign2d.shape != q.shape
        or rfs.shape != q.shape + (3,)
        or rfe.shape != q.shape + (3,)
    ):
        raise ValueError(f"Inconsistent FastQSL2 array shapes in {path}")

    logq = np.full(q.shape, np.nan)
    valid = np.isfinite(q) & (q > 0.0)
    logq[valid] = np.log10(q[valid])
    closed = rboundary == 11
    endpoint_bottom_bottom = np.isclose(rfs[..., 2], 0.0, atol=1.0e-3) & np.isclose(
        rfe[..., 2], 0.0, atol=1.0e-3
    )
    if np.any(closed & ~endpoint_bottom_bottom):
        count = int(np.count_nonzero(closed & ~endpoint_bottom_bottom))
        raise ValueError(f"{count} code-11 lines lack two bottom terminals in {path}")

    ny, nx = q.shape
    x_native = np.linspace(xreg[0], xreg[1], nx)
    y_native = np.linspace(yreg[0], yreg[1], ny)
    x_mm = native_to_mm(x_native, xreg, PHYSICAL_EXTENT_MM[:2])
    y_mm = native_to_mm(y_native, yreg, PHYSICAL_EXTENT_MM[2:])
    return {
        "q": q,
        "logq": logq,
        "rboundary": rboundary,
        "sign2d": sign2d,
        "closed": closed,
        "xreg": xreg,
        "yreg": yreg,
        "x_native": x_native,
        "y_native": y_native,
        "x_mm": x_mm,
        "y_mm": y_mm,
        "delta": delta,
    }


def native_to_mm(values: np.ndarray, native_range: np.ndarray, physical_range: tuple[float, float]):
    """Map native grid coordinates linearly to the plotted AMRVAC coordinates."""

    scale = (physical_range[1] - physical_range[0]) / (native_range[1] - native_range[0])
    return physical_range[0] + (values - native_range[0]) * scale


def select_candidate(feature_path: Path, snapshot: int, qsl: dict) -> dict:
    """Apply the documented 3-D axis/current-channel selection."""

    if snapshot not in SELECTION_PARAMETERS:
        raise KeyError(f"No current-channel selection parameters for snapshot {snapshot}")
    with np.load(feature_path) as archive:
        required = {"seeds", "bzseed", "status", "features", "ftw", "rb", "axis"}
        missing = required.difference(archive.files)
        if missing:
            raise KeyError(f"{feature_path} is missing feature arrays: {sorted(missing)}")
        seeds = np.asarray(archive["seeds"], dtype=float)
        bzseed = np.asarray(archive["bzseed"], dtype=float)
        status = np.asarray(archive["status"], dtype=np.int8)
        features = np.asarray(archive["features"], dtype=float)
        ftw = np.asarray(archive["ftw"], dtype=float)
        rb = np.asarray(archive["rb"], dtype=np.int16)
        axis = np.asarray(archive["axis"], dtype=float)

    if features.shape != (len(seeds), len(FEATURE_COLUMNS)):
        raise ValueError(f"Unexpected feature shape {features.shape} in {feature_path}")
    for name, array in (("bzseed", bzseed), ("status", status), ("ftw", ftw), ("rb", rb)):
        if len(array) != len(seeds):
            raise ValueError(f"{name} length does not match seeds in {feature_path}")

    p = SELECTION_PARAMETERS[snapshot]
    tol = p["edge_tolerance"]
    selected = (
        (status == 1)
        & np.isfinite(features[:, 4])
        & (features[:, 2] <= -(p["a_minus"] - tol))
        & (features[:, 3] >= p["a_plus"] - tol)
        & (features[:, 4] >= p["height"] + tol)
        & (np.abs(features[:, 7]) < p["winding_abs_max"])
        & (features[:, 6] >= p["axis_coverage_min"])
    )
    if not np.any(selected):
        raise RuntimeError(f"The 3-D selection returned no seeds for snapshot {snapshot}")
    if np.any(rb[selected] != 11):
        raise ValueError("Candidate includes a field line without two bottom terminals")
    if np.any(bzseed[selected] <= 0.0):
        raise ValueError("Candidate must be defined on the positive-polarity footpoints")

    x_native = np.arange(
        qsl["xreg"][0], qsl["xreg"][1] + 0.5 * SEED_STEP_NATIVE, SEED_STEP_NATIVE
    )
    y_native = np.arange(
        qsl["yreg"][0], qsl["yreg"][1] + 0.5 * SEED_STEP_NATIVE, SEED_STEP_NATIVE
    )
    ix = np.rint((seeds[selected, 0] - x_native[0]) / SEED_STEP_NATIVE).astype(int)
    iy = np.rint((seeds[selected, 1] - y_native[0]) / SEED_STEP_NATIVE).astype(int)
    if np.any(ix < 0) or np.any(ix >= len(x_native)) or np.any(iy < 0) or np.any(iy >= len(y_native)):
        raise IndexError("Selected seed lies outside the bottom plane")
    if not np.allclose(x_native[ix], seeds[selected, 0], atol=1.0e-8) or not np.allclose(
        y_native[iy], seeds[selected, 1], atol=1.0e-8
    ):
        raise ValueError("Selected seeds are not aligned with the factor-2 grid")

    mask = np.zeros((len(y_native), len(x_native)), dtype=bool)
    mask[iy, ix] = True
    x_mm = native_to_mm(x_native, qsl["xreg"], PHYSICAL_EXTENT_MM[:2])
    y_mm = native_to_mm(y_native, qsl["yreg"], PHYSICAL_EXTENT_MM[2:])
    seed_x_mm = native_to_mm(seeds[selected, 0], qsl["xreg"], PHYSICAL_EXTENT_MM[:2])
    seed_y_mm = native_to_mm(seeds[selected, 1], qsl["yreg"], PHYSICAL_EXTENT_MM[2:])

    # The Q product has twice the sampling frequency of the factor-2 feature grid.
    qix = np.rint((seeds[selected, 0] - qsl["xreg"][0]) / qsl["delta"]).astype(int)
    qiy = np.rint((seeds[selected, 1] - qsl["yreg"][0]) / qsl["delta"]).astype(int)
    rb_qsl = qsl["rboundary"][qiy, qix]
    if np.any(rb_qsl != 11):
        raise ValueError("Candidate disagrees with the QSL product's closed-field mask")

    seed_area_cm2 = (SEED_STEP_NATIVE * NATIVE_CELL_MM * 1.0e8) ** 2
    positive_flux_mx = float(np.sum(bzseed[selected]) * B_UNIT_G * seed_area_cm2)
    return {
        "selected": selected,
        "seeds": seeds,
        "bzseed": bzseed,
        "features": features,
        "ftw": ftw,
        "rb": rb,
        "axis": axis,
        "mask": mask,
        "x_native": x_native,
        "y_native": y_native,
        "x_mm": x_mm,
        "y_mm": y_mm,
        "seed_x_mm": seed_x_mm,
        "seed_y_mm": seed_y_mm,
        "positive_flux_mx": positive_flux_mx,
        "seed_area_cm2": seed_area_cm2,
        "q_values": qsl["q"][qiy, qix],
    }


def save_candidate(path: Path, snapshot: int, candidate: dict, source: Path, qsl_source: Path) -> None:
    """Save the selected footprint and enough metadata to reuse it."""

    selected = candidate["selected"]
    np.savez_compressed(
        path,
        mask=candidate["mask"].astype(np.uint8),
        x_native=candidate["x_native"].astype(np.float32),
        y_native=candidate["y_native"].astype(np.float32),
        x_Mm=candidate["x_mm"].astype(np.float32),
        y_Mm=candidate["y_mm"].astype(np.float32),
        selected_seeds_native=candidate["seeds"][selected].astype(np.float32),
        selected_seeds_Mm=np.column_stack(
            (candidate["seed_x_mm"], candidate["seed_y_mm"])
        ).astype(np.float32),
        selected_bz_code=candidate["bzseed"][selected].astype(np.float32),
        selected_features=candidate["features"][selected].astype(np.float32),
        feature_column_names=np.asarray(FEATURE_COLUMNS),
        selected_ftw=candidate["ftw"][selected].astype(np.float32),
        magnetic_axis_native=candidate["axis"].astype(np.float32),
        positive_flux_Mx=np.asarray(candidate["positive_flux_mx"]),
        seed_area_cm2=np.asarray(candidate["seed_area_cm2"]),
        B_unit_G=np.asarray(B_UNIT_G),
        selection_parameters_json=np.asarray(json.dumps(SELECTION_PARAMETERS[snapshot])),
        feature_source=np.asarray(str(source.resolve())),
        qsl_source=np.asarray(str(qsl_source.resolve())),
        definition=np.asarray(
            "candidate positive footpoints of field lines that cross and overlie the 3-D current channel"
        ),
    )


def draw_layers(axis, dataset: dict, *, limits: tuple[float, float, float, float]):
    """Draw Q, all closed field, and the selected integration footprint."""

    qsl = dataset["qsl"]
    candidate = dataset["candidate"]
    image = axis.imshow(
        qsl["logq"],
        origin="lower",
        extent=PHYSICAL_EXTENT_MM,
        cmap="gray_r",
        vmin=0.3,
        vmax=6.5,
        interpolation="nearest",
        rasterized=True,
        zorder=1,
    )
    closed_layer = np.ma.masked_where(~qsl["closed"], np.ones_like(qsl["closed"], dtype=float))
    axis.imshow(
        closed_layer,
        origin="lower",
        extent=PHYSICAL_EXTENT_MM,
        cmap=ListedColormap([CLOSED_FILL]),
        vmin=0.0,
        vmax=1.0,
        interpolation="nearest",
        alpha=0.32,
        rasterized=True,
        zorder=2,
    )
    axis.contour(
        qsl["x_mm"],
        qsl["y_mm"],
        qsl["closed"].astype(float),
        levels=[0.5],
        colors=[CLOSED_EDGE],
        linewidths=0.9,
        zorder=3.0,
    )

    masked_candidate = np.ma.masked_where(~candidate["mask"], np.ones_like(candidate["mask"], dtype=float))
    axis.pcolormesh(
        candidate["x_mm"],
        candidate["y_mm"],
        masked_candidate,
        cmap=ListedColormap([CANDIDATE_FILL]),
        vmin=0.0,
        vmax=1.0,
        shading="nearest",
        alpha=0.88,
        rasterized=True,
        zorder=3.3,
    )
    axis.contour(
        candidate["x_mm"],
        candidate["y_mm"],
        candidate["mask"].astype(float),
        levels=[0.5],
        colors=[CANDIDATE_EDGE],
        linewidths=0.9,
        zorder=3.5,
    )

    contour_paths = {}
    for level, color, style in zip(LOGQ_LEVELS, Q_COLORS, Q_LINESTYLES):
        if level == 2.0:
            axis.contour(
                qsl["x_mm"],
                qsl["y_mm"],
                qsl["logq"],
                levels=[level],
                colors=["white"],
                linestyles=[style],
                linewidths=3.4,
                antialiased=True,
                zorder=3.9,
            )
        contour = axis.contour(
            qsl["x_mm"],
            qsl["y_mm"],
            qsl["logq"],
            levels=[level],
            colors=[color],
            linestyles=[style],
            linewidths=1.7 if level == 2.0 else 1.35,
            antialiased=True,
            zorder=4.0,
        )
        contour_paths[int(level)] = [segment for segment in contour.allsegs[0] if len(segment) > 1]

    snapshot = dataset["snapshot"]
    axis.set_title(f"t = {snapshot * TIME_MINUTES_PER_SNAPSHOT:.1f} min", fontsize=16, pad=8)
    axis.set_xlim(limits[0], limits[1])
    axis.set_ylim(limits[2], limits[3])
    axis.set_aspect("equal")
    axis.tick_params(direction="in", top=True, right=True, labelsize=12)
    return image, contour_paths


def distances_to_contours(points: np.ndarray, paths: list[np.ndarray]) -> np.ndarray:
    """Return exact minimum point-to-polyline distances in plotted Mm."""

    segment_starts = []
    segment_ends = []
    for path in paths:
        if len(path) > 1:
            segment_starts.append(path[:-1])
            segment_ends.append(path[1:])
    if not segment_starts:
        return np.full(len(points), np.nan)
    starts = np.concatenate(segment_starts, axis=0)
    vectors = np.concatenate(segment_ends, axis=0) - starts
    lengths2 = np.einsum("ij,ij->i", vectors, vectors)
    valid = lengths2 > 0.0
    starts = starts[valid]
    vectors = vectors[valid]
    lengths2 = lengths2[valid]
    result = []
    for point in points:
        offsets = point - starts
        t = np.clip(np.einsum("ij,ij->i", offsets, vectors) / lengths2, 0.0, 1.0)
        nearest = starts + t[:, None] * vectors
        result.append(float(np.min(np.linalg.norm(point - nearest, axis=1))))
    return np.asarray(result)


def make_figure(datasets: list[dict], output_base: Path, limits, dpi: int, *, legend: bool):
    """Create one shared-axis, two-time figure and return contour geometry."""

    figure, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(17.5, 8.0 if legend else 7.2),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    all_paths = []
    image = None
    for index, (axis, dataset) in enumerate(zip(axes, datasets)):
        image, paths = draw_layers(axis, dataset, limits=limits)
        all_paths.append(paths)
        axis.set_xlabel("Solar-X [Mm]", fontsize=14)
        if index == 0:
            axis.set_ylabel("Solar-Y [Mm]", fontsize=14)

    handles = [
        Patch(facecolor=CLOSED_FILL, edgecolor=CLOSED_EDGE, alpha=0.55, label="All closed field"),
        Patch(
            facecolor=CANDIDATE_FILL,
            edgecolor=CANDIDATE_EDGE,
            alpha=0.9,
            label=r"Candidate $\Phi_{\rm over}$ footpoints",
        ),
    ]
    handles.extend(
        Line2D(
            [0],
            [0],
            color=color,
            linestyle=style,
            linewidth=2.2,
            label=rf"$\log_{{10}}Q={int(level)}$",
        )
        for level, color, style in zip(LOGQ_LEVELS, Q_COLORS, Q_LINESTYLES)
    )
    if legend:
        figure.legend(
            handles=handles,
            loc="outside lower center",
            ncol=4,
            frameon=False,
            handlelength=3.0,
            columnspacing=1.4,
            fontsize=11,
        )
    colorbar = figure.colorbar(image, ax=axes, shrink=0.91, pad=0.02)
    colorbar.set_label(r"Bottom-boundary $\log_{10}Q$", fontsize=13)
    colorbar.ax.tick_params(labelsize=11)
    figure.savefig(output_base.with_suffix(".png"), dpi=dpi, bbox_inches="tight")
    figure.savefig(output_base.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(figure)
    return all_paths


def place_polarity_labels(axis, qsl: dict, limits: tuple[float, float, float, float]) -> None:
    """Place one large sign well inside each visible closed-polarity region."""

    visible = (
        (qsl["y_mm"][:, None] > limits[2])
        & (qsl["y_mm"][:, None] < limits[3])
        & (qsl["x_mm"][None, :] > limits[0])
        & (qsl["x_mm"][None, :] < limits[1])
    )
    dx = float(qsl["x_mm"][1] - qsl["x_mm"][0])
    dy = float(qsl["y_mm"][1] - qsl["y_mm"][0])
    for mask, sign in (
        (qsl["closed"] & (qsl["sign2d"] > 0), "+"),
        (qsl["closed"] & (qsl["sign2d"] < 0), "-"),
    ):
        interior = mask & visible
        if not np.any(interior):
            continue
        clearance = distance_transform_edt(interior, sampling=(dy, dx))
        iy, ix = np.unravel_index(np.argmax(clearance), clearance.shape)
        axis.text(
            qsl["x_mm"][ix],
            qsl["y_mm"][iy],
            rf"${sign}$",
            fontsize=32,
            fontweight="bold",
            color="#202020",
            ha="center",
            va="center",
            zorder=6,
            clip_on=True,
        )


def draw_white_q2_layers(
    axis,
    dataset: dict,
    *,
    limits: tuple[float, float, float, float],
    annotate_regions: bool,
):
    """Match the final white closed/Q paper style and add the candidate mask."""

    qsl = dataset["qsl"]
    candidate = dataset["candidate"]
    positive = qsl["closed"] & (qsl["sign2d"] > 0)
    negative = qsl["closed"] & (qsl["sign2d"] < 0)
    if not np.array_equal(positive | negative, qsl["closed"]):
        raise ValueError("sign2d does not split every closed-field pixel into +/- polarity")

    axis.set_facecolor("white")
    categories = np.zeros_like(qsl["sign2d"], dtype=float)
    categories[positive] = 1.0
    categories[negative] = 0.0
    signed_closed = np.ma.masked_where(~qsl["closed"], categories)
    axis.pcolormesh(
        qsl["x_mm"],
        qsl["y_mm"],
        signed_closed,
        cmap=ListedColormap([PAPER_NEGATIVE_FILL, PAPER_POSITIVE_FILL]),
        vmin=-0.5,
        vmax=1.5,
        shading="nearest",
        alpha=0.30,
        rasterized=True,
        zorder=2,
    )
    axis.contour(
        qsl["x_mm"],
        qsl["y_mm"],
        qsl["closed"].astype(float),
        levels=[0.5],
        colors=[CLOSED_EDGE],
        linewidths=0.95,
        zorder=3.0,
    )

    masked_candidate = np.ma.masked_where(~candidate["mask"], np.ones_like(candidate["mask"], dtype=float))
    axis.pcolormesh(
        candidate["x_mm"],
        candidate["y_mm"],
        masked_candidate,
        cmap=ListedColormap([PAPER_CANDIDATE_FILL]),
        vmin=0.0,
        vmax=1.0,
        shading="nearest",
        alpha=0.82,
        rasterized=True,
        zorder=3.3,
    )
    axis.contour(
        candidate["x_mm"],
        candidate["y_mm"],
        candidate["mask"].astype(float),
        levels=[0.5],
        colors=[PAPER_CANDIDATE_EDGE],
        linewidths=1.05,
        zorder=3.5,
    )

    axis.contour(
        qsl["x_mm"],
        qsl["y_mm"],
        qsl["logq"],
        levels=[2.0],
        colors=["white"],
        linewidths=3.3,
        zorder=3.9,
    )
    q2 = axis.contour(
        qsl["x_mm"],
        qsl["y_mm"],
        qsl["logq"],
        levels=[2.0],
        colors=[Q_COLORS[1]],
        linewidths=1.9,
        zorder=4.0,
    )
    # FastQSL2's sign2d is the sign of the bottom normal field used by the Q
    # calculation.  Its zero interface reproduces Bz=0 to Q-grid accuracy even
    # when the native magnetic volume is offline.
    axis.contour(
        qsl["x_mm"],
        qsl["y_mm"],
        qsl["sign2d"].astype(float),
        levels=[0.0],
        colors=[PAPER_BZ0_COLOR],
        linestyles="dashed",
        linewidths=1.35,
        zorder=5.0,
    )

    if annotate_regions:
        place_polarity_labels(axis, qsl, limits)
        open_x = limits[0] + 0.95 * (limits[1] - limits[0])
        open_y = limits[2] + 0.95 * (limits[3] - limits[2])
        minus_x = limits[0] + 0.87 * (limits[1] - limits[0])
        minus_y = limits[2] + 0.12 * (limits[3] - limits[2])
        axis.text(
            open_x,
            open_y,
            "Open field",
            fontsize=20,
            color="#202020",
            ha="right",
            va="top",
            zorder=6,
        )
        axis.text(
            minus_x,
            minus_y,
            "$-$",
            fontsize=34,
            color="#202020",
            ha="center",
            va="center",
            zorder=6,
        )

    snapshot = dataset["snapshot"]
    axis.set_title(f"t = {snapshot * TIME_MINUTES_PER_SNAPSHOT:.1f} min", fontsize=20, pad=9)
    axis.set_xlim(limits[:2])
    axis.set_ylim(limits[2:])
    axis.set_aspect("equal")
    axis.tick_params(direction="in", top=True, right=True, labelsize=16)
    return [segment for segment in q2.allsegs[0] if len(segment) > 1]


def make_white_q2_figure(
    datasets: list[dict],
    output_base: Path,
    limits: tuple[float, float, float, float],
    dpi: int,
    *,
    annotate_regions: bool,
) -> None:
    """Create the paper-style white-background comparison."""

    figure, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(16.8, 9.0),
        sharex=True,
        sharey=True,
        layout="none",
    )
    figure.set_facecolor("white")
    figure.subplots_adjust(left=0.075, right=0.985, bottom=0.15, top=0.885, wspace=0.07)
    axes = np.atleast_1d(axes)
    for index, (axis, dataset) in enumerate(zip(axes, datasets)):
        draw_white_q2_layers(
            axis,
            dataset,
            limits=limits,
            annotate_regions=annotate_regions,
        )
        axis.set_xlabel("Solar-X [Mm]", fontsize=18)
        if index == 0:
            axis.set_ylabel("Solar-Y [Mm]", fontsize=18)
        axis.tick_params(axis="y", labelleft=index == 0, labelright=False)

    handles = [
        Patch(
            facecolor=PAPER_POSITIVE_FILL,
            edgecolor=CLOSED_EDGE,
            alpha=0.30,
            label=r"Closed field ($+$)",
        ),
        Patch(
            facecolor=PAPER_NEGATIVE_FILL,
            edgecolor=CLOSED_EDGE,
            alpha=0.30,
            label=r"Closed field ($-$)",
        ),
        Patch(
            facecolor=PAPER_CANDIDATE_FILL,
            edgecolor=PAPER_CANDIDATE_EDGE,
            alpha=0.82,
            label=r"Candidate $\Phi_{\rm over}$",
        ),
        Line2D([0], [0], color=Q_COLORS[1], linewidth=2.3, label=r"$\log_{10}Q=2$"),
        Line2D(
            [0],
            [0],
            color=PAPER_BZ0_COLOR,
            linestyle="dashed",
            linewidth=2.0,
            label=r"$B_z=0$",
        ),
    ]
    figure.legend(
        handles=handles,
        loc="lower center",
        bbox_to_anchor=(0.5, 0.005),
        ncol=5,
        frameon=False,
        handlelength=2.1,
        columnspacing=1.05,
        fontsize=16,
    )
    figure.savefig(
        output_base.with_suffix(".png"),
        dpi=dpi,
        bbox_inches="tight",
        facecolor="white",
        transparent=False,
    )
    figure.savefig(
        output_base.with_suffix(".pdf"),
        bbox_inches="tight",
        facecolor="white",
        transparent=False,
    )
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--qsl-dir", type=Path, default=QSL_DIR)
    parser.add_argument("--feature-dir", type=Path, default=FEATURE_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--snapshots", type=int, nargs="+", default=list(SNAPSHOTS))
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    datasets = []
    for snapshot in args.snapshots:
        qsl_source = args.qsl_dir / f"bottom_fastqsl2_t{snapshot:04d}_delta0.25.npz"
        feature_source = args.feature_dir / f"axis_features_factor2_t{snapshot:04d}.npz"
        if not qsl_source.is_file():
            raise FileNotFoundError(qsl_source)
        if not feature_source.is_file():
            raise FileNotFoundError(feature_source)
        qsl = load_bottom_result(qsl_source)
        candidate = select_candidate(feature_source, snapshot, qsl)
        selected_output = args.output_dir / f"candidate_overlying_positive_t{snapshot:04d}.npz"
        save_candidate(selected_output, snapshot, candidate, feature_source, qsl_source)
        datasets.append(
            {
                "snapshot": snapshot,
                "qsl": qsl,
                "candidate": candidate,
                "qsl_source": qsl_source,
                "feature_source": feature_source,
                "candidate_output": selected_output,
            }
        )

    suffix = "_".join(f"t{snapshot:04d}" for snapshot in args.snapshots)
    main_base = args.output_dir / f"bottom_closed_q_overlying_candidate_{suffix}"
    zoom_base = args.output_dir / f"bottom_closed_q_overlying_candidate_zoom_{suffix}"
    paper_base = args.output_dir / f"bottom_white_q2_closed_overlying_candidate_{suffix}"
    paper_zoom_base = args.output_dir / f"bottom_white_q2_closed_overlying_candidate_zoom_{suffix}"
    contour_paths = make_figure(datasets, main_base, MAIN_LIMITS_MM, args.dpi, legend=True)
    make_figure(datasets, zoom_base, ZOOM_LIMITS_MM, args.dpi, legend=True)
    make_white_q2_figure(
        datasets,
        paper_base,
        MAIN_LIMITS_MM,
        args.dpi,
        annotate_regions=True,
    )
    make_white_q2_figure(
        datasets,
        paper_zoom_base,
        ZOOM_LIMITS_MM,
        args.dpi,
        annotate_regions=False,
    )

    metadata = {
        "method": {
            "name": "3-D magnetic-axis/current-channel candidate overlying-flux footprint",
            "interpretation": (
                "positive-polarity footpoints of closed field lines that cross both current-channel "
                "flanks, pass above the central magnetic axis, cover >=25% of the central axis, "
                "and have abs(winding number)<1"
            ),
            "status": "candidate integration region; not a universal topological boundary",
            "feature_columns": list(FEATURE_COLUMNS),
            "B_unit_G": B_UNIT_G,
            "native_cell_Mm": NATIVE_CELL_MM,
            "factor2_seed_step_native": SEED_STEP_NATIVE,
            "negative_conjugate_flux_check": (
                "pending because the raw magnetic-field volume was unavailable when this figure was made"
            ),
            "paper_style_Bz0_source": (
                "zero interface of FastQSL2 sign2d; the exact native-Bz contour should be regenerated "
                "after the magnetic volume is remounted"
            ),
        },
        "layers": {
            "background": "bottom-boundary log10(Q)",
            "purple": "all FastQSL2 rboundary=11 bottom-to-bottom closed field",
            "yellow": "candidate positive-polarity Phi_over integration footprint",
            "contours_log10_Q": LOGQ_LEVELS.tolist(),
        },
        "limits_Mm": {"main": list(MAIN_LIMITS_MM), "zoom": list(ZOOM_LIMITS_MM)},
        "outputs": {
            "main_png": str(main_base.with_suffix(".png").resolve()),
            "main_pdf": str(main_base.with_suffix(".pdf").resolve()),
            "zoom_png": str(zoom_base.with_suffix(".png").resolve()),
            "zoom_pdf": str(zoom_base.with_suffix(".pdf").resolve()),
            "paper_white_q2_png": str(paper_base.with_suffix(".png").resolve()),
            "paper_white_q2_pdf": str(paper_base.with_suffix(".pdf").resolve()),
            "paper_white_q2_zoom_png": str(paper_zoom_base.with_suffix(".png").resolve()),
            "paper_white_q2_zoom_pdf": str(paper_zoom_base.with_suffix(".pdf").resolve()),
        },
        "snapshots": [],
    }

    for dataset, paths in zip(datasets, contour_paths):
        candidate = dataset["candidate"]
        selected = candidate["selected"]
        points = np.column_stack((candidate["seed_x_mm"], candidate["seed_y_mm"]))
        q2_dist = distances_to_contours(points, paths[2])
        entry = {
            "snapshot": dataset["snapshot"],
            "time_minutes": dataset["snapshot"] * TIME_MINUTES_PER_SNAPSHOT,
            "parameters": SELECTION_PARAMETERS[dataset["snapshot"]],
            "selected_seed_count": int(np.count_nonzero(selected)),
            "positive_flux_Mx": candidate["positive_flux_mx"],
            "positive_flux_1e18_Mx": candidate["positive_flux_mx"] / 1.0e18,
            "footpoint_bbox_Mm": {
                "x": [float(np.min(points[:, 0])), float(np.max(points[:, 0]))],
                "y": [float(np.min(points[:, 1])), float(np.max(points[:, 1]))],
            },
            "distance_to_log10Q2_Mm": {
                "median": float(np.nanmedian(q2_dist)),
                "p90": float(np.nanpercentile(q2_dist, 90.0)),
                "maximum": float(np.nanmax(q2_dist)),
                "fraction_le_0p5": float(np.mean(q2_dist <= 0.5)),
                "fraction_le_1": float(np.mean(q2_dist <= 1.0)),
                "fraction_le_2": float(np.mean(q2_dist <= 2.0)),
            },
            "selected_q_range": [
                float(np.nanmin(candidate["q_values"])),
                float(np.nanmax(candidate["q_values"])),
            ],
            "all_selected_lines_closed": bool(np.all(candidate["rb"][selected] == 11)),
            "all_selected_footpoints_positive": bool(np.all(candidate["bzseed"][selected] > 0.0)),
            "qsl_source": str(dataset["qsl_source"].resolve()),
            "feature_source": str(dataset["feature_source"].resolve()),
            "candidate_npz": str(dataset["candidate_output"].resolve()),
        }
        metadata["snapshots"].append(entry)

    metadata_path = args.output_dir / f"bottom_closed_q_overlying_candidate_{suffix}.json"
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(paper_base.with_suffix(".png").resolve())
    print(paper_zoom_base.with_suffix(".png").resolve())
    print(main_base.with_suffix(".png").resolve())
    print(zoom_base.with_suffix(".png").resolve())
    print(metadata_path.resolve())


if __name__ == "__main__":
    main()

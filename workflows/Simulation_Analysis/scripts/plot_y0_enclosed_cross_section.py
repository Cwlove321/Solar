#!/usr/bin/env python3
"""Fill a y=0 magnetic cross-section that contains the magnetic axis.

Two reproducible selections are available.  ``qsl-envelope`` uses the
four-connected component satisfying all of:

* both field-line terminals reach the bottom boundary (rboundary == 11),
* log10(Q) is below the selected ridge threshold (default 2), and
* the component contains the magnetic-axis intersection with the y=0 plane.

When the high-Q ridge is open at its lower end, the lowest Q-sampling row is
treated as connected to the physical bottom boundary.  On every intersected
height row, the hybrid envelope runs from the leftmost point of the closed
component to the rightmost point of the axis-containing low-Q component.  Thus
the left/top edge is always the ``rboundary == 11`` connectivity boundary and
only the central high-Q ridge acts as the right edge.  The envelope includes
internal or left-side secondary high-Q layers.  Finally, holes fully enclosed
by this hybrid boundary are filled geometrically; the saved NPZ also retains a
separate bottom-to-bottom-closed mask so those pixels remain auditable.

``closed-component`` instead fills the complete four-connected
``rboundary == 11`` component containing the magnetic axis.  This second
selection is the whole bottom-to-bottom closed dome, not only flux directly
overlying the current channel.
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
from scipy.ndimage import binary_fill_holes, label


INPUT_DIR = Path("/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/y0_qsl")
OUTPUT_DIR = INPUT_DIR / "enclosed_region"
SNAPSHOTS = (10, 95)
TIME_MINUTES_PER_SNAPSHOT = 0.71558
VIEW_LIMITS_MM = (-15.0, 30.0, 0.0, 30.0)
CLOSED_VIEW_LIMITS_MM = (-10.0, 47.5, 0.0, 25.0)

LOGQ_MIN = float(np.log10(2.0))
LOGQ_MAX = 8.0
Q2_COLOR = "#E69F00"
CLOSED_COLOR = "#6A3D9A"
REGION_FILL = "#43C7CA"
REGION_EDGE = "#006B70"
AXIS_COLOR = "#F0E442"
FOUR_CONNECTED = np.asarray([[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.uint8)


def grid_vectors(x_mm: np.ndarray, height_mm: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Return 1-D x and height vectors from FastQSL2 seed-coordinate arrays."""

    x = np.asarray(x_mm, dtype=float)
    height = np.asarray(height_mm, dtype=float)
    if x.ndim == 2:
        x = x[0, :]
    if height.ndim == 2:
        height = height[:, 0]
    if x.ndim != 1 or height.ndim != 1:
        raise ValueError("Expected separable fixed-y x and height coordinates")
    return x, height


def component_containing_axis(
    logq: np.ndarray,
    rboundary: np.ndarray,
    x: np.ndarray,
    height: np.ndarray,
    axis_point: np.ndarray,
    logq_cut: float,
) -> tuple[np.ndarray, dict]:
    """Extract the conservative four-connected low-Q closed component."""

    ix = int(np.argmin(np.abs(x - axis_point[0])))
    iz = int(np.argmin(np.abs(height - axis_point[1])))
    eligible = (rboundary == 11) & np.isfinite(logq) & (logq < logq_cut)
    if not eligible[iz, ix]:
        raise ValueError(
            f"Magnetic-axis grid point is outside the eligible set: "
            f"log10Q={logq[iz, ix]}, rboundary={rboundary[iz, ix]}"
        )
    components, component_count = label(eligible, structure=FOUR_CONNECTED)
    component_id = int(components[iz, ix])
    if component_id == 0:
        raise RuntimeError("The magnetic-axis component was not labeled")
    raw = components == component_id

    # The purple connectivity contour, not a secondary high-Q ridge, defines
    # the left/top boundary.  Select the complete closed component containing
    # the magnetic axis, then extend each row of the raw low-Q basin leftward
    # to that component's outer edge.  The raw basin's rightmost point retains
    # the central high-Q ridge as the right boundary.
    closed_components, closed_component_count = label(
        rboundary == 11, structure=FOUR_CONNECTED
    )
    closed_component_id = int(closed_components[iz, ix])
    closed_component = closed_components == closed_component_id
    geometric_span = np.zeros_like(raw)
    for row in range(raw.shape[0]):
        raw_columns = np.flatnonzero(raw[row])
        closed_columns = np.flatnonzero(closed_component[row])
        if raw_columns.size and closed_columns.size:
            geometric_span[
                row, closed_columns.min() : raw_columns.max() + 1
            ] = True
    closed_only = geometric_span & closed_component
    # The requested region is a solid geometric cross-section.  Fill only
    # holes completely enclosed by the hybrid boundary; intersecting with the
    # geometric span guarantees that neither outer boundary can move.
    filled = binary_fill_holes(closed_only) & geometric_span
    internal_hole_fill = filled & ~closed_only
    added = closed_only & ~raw
    excluded_open = geometric_span & ~filled

    selected_rows = np.flatnonzero(np.any(filled, axis=1))
    left_boundary_matches = 0
    right_boundary_matches = 0
    for row in selected_rows:
        columns = np.flatnonzero(filled[row])
        left = int(columns.min())
        right = int(columns.max())
        if left == 0 or not closed_component[row, left - 1]:
            left_boundary_matches += 1
        if right + 1 < filled.shape[1]:
            next_logq = logq[row, right + 1]
            if (
                closed_component[row, right + 1]
                and np.isfinite(next_logq)
                and next_logq >= logq_cut
            ):
                right_boundary_matches += 1

    izs, ixs = np.where(filled)
    dx = float(x[1] - x[0])
    dz = float(height[1] - height[0])
    statistics = {
        "log10Q_cut": logq_cut,
        "connectivity": 4,
        "eligible_component_count": int(component_count),
        "component_id": component_id,
        "closed_component_count": int(closed_component_count),
        "closed_component_id": closed_component_id,
        "raw_pixel_count": int(np.count_nonzero(raw)),
        "geometric_span_pixel_count": int(np.count_nonzero(geometric_span)),
        "closed_only_pixel_count": int(np.count_nonzero(closed_only)),
        "filled_pixel_count": int(np.count_nonzero(filled)),
        "added_closed_internal_pixel_count": int(np.count_nonzero(added)),
        "filled_internal_hole_pixel_count": int(np.count_nonzero(internal_hole_fill)),
        "remaining_excluded_span_pixel_count": int(np.count_nonzero(excluded_open)),
        "selected_height_row_count": int(selected_rows.size),
        "left_rows_matching_closed_boundary": int(left_boundary_matches),
        "right_rows_matching_high_Q_ridge": int(right_boundary_matches),
        "area_Mm2": float(np.count_nonzero(filled) * dx * dz),
        "bbox_Mm": {
            "x": [float(x[ixs].min()), float(x[ixs].max())],
            "height": [float(height[izs].min()), float(height[izs].max())],
        },
        "axis_grid_point_Mm": [float(x[ix]), float(height[iz])],
        "axis_grid_log10Q": float(logq[iz, ix]),
        "touches_lowest_Q_row": bool(np.any(filled[0, :])),
        "touches_top_Q_row": bool(np.any(filled[-1, :])),
        "touches_left_Q_edge": bool(np.any(filled[:, 0])),
        "touches_right_Q_edge": bool(np.any(filled[:, -1])),
    }
    return filled, statistics


def closed_component_containing_axis(
    rboundary: np.ndarray,
    x: np.ndarray,
    height: np.ndarray,
    axis_point: np.ndarray,
) -> tuple[np.ndarray, dict]:
    """Extract the bottom-to-bottom closed component containing the axis."""

    ix = int(np.argmin(np.abs(x - axis_point[0])))
    iz = int(np.argmin(np.abs(height - axis_point[1])))
    closed = rboundary == 11
    if not closed[iz, ix]:
        raise ValueError(
            f"Magnetic-axis grid point is not bottom-to-bottom closed: "
            f"rboundary={rboundary[iz, ix]}"
        )
    components, component_count = label(closed, structure=FOUR_CONNECTED)
    component_id = int(components[iz, ix])
    if component_id == 0:
        raise RuntimeError("The magnetic-axis closed-field component was not labeled")
    mask = components == component_id
    izs, ixs = np.where(mask)
    dx = float(x[1] - x[0])
    dz = float(height[1] - height[0])
    statistics = {
        "connectivity": 4,
        "closed_component_count": int(component_count),
        "component_id": component_id,
        "filled_pixel_count": int(np.count_nonzero(mask)),
        "area_Mm2": float(np.count_nonzero(mask) * dx * dz),
        "bbox_Mm": {
            "x": [float(x[ixs].min()), float(x[ixs].max())],
            "height": [float(height[izs].min()), float(height[izs].max())],
        },
        "axis_grid_point_Mm": [float(x[ix]), float(height[iz])],
        "touches_lowest_Q_row": bool(np.any(mask[0, :])),
        "touches_top_Q_row": bool(np.any(mask[-1, :])),
        "touches_left_Q_edge": bool(np.any(mask[:, 0])),
        "touches_right_Q_edge": bool(np.any(mask[:, -1])),
    }
    return mask, statistics


def load_snapshot(
    path: Path,
    snapshot: int,
    logq_cut: float,
    selection_mode: str,
) -> dict:
    with np.load(path) as archive:
        required = {
            "logq",
            "rboundary",
            "x_mm",
            "height_mm",
            "magnetic_axis_intersections_Mm",
        }
        missing = required.difference(archive.files)
        if missing:
            raise KeyError(f"{path} is missing arrays: {sorted(missing)}")
        logq = np.asarray(archive["logq"], dtype=float)
        rboundary = np.asarray(archive["rboundary"], dtype=np.int16)
        x_mm = np.asarray(archive["x_mm"], dtype=float)
        height_mm = np.asarray(archive["height_mm"], dtype=float)
        axis_hits = np.asarray(archive["magnetic_axis_intersections_Mm"], dtype=float)
    if logq.shape != rboundary.shape or x_mm.shape != logq.shape or height_mm.shape != logq.shape:
        raise ValueError(f"Inconsistent fixed-y slice shapes in {path}")
    if axis_hits.shape != (1, 2):
        raise ValueError(f"Expected one magnetic-axis intersection in {path}; got {axis_hits.shape}")
    x, height = grid_vectors(x_mm, height_mm)
    if selection_mode == "qsl-envelope":
        mask, statistics = component_containing_axis(
            logq, rboundary, x, height, axis_hits[0], logq_cut
        )
    elif selection_mode == "closed-component":
        mask, statistics = closed_component_containing_axis(
            rboundary, x, height, axis_hits[0]
        )
    else:
        raise ValueError(f"Unknown selection mode: {selection_mode}")
    return {
        "snapshot": snapshot,
        "source": path,
        "logq": logq,
        "rboundary": rboundary,
        "x": x,
        "height": height,
        "axis_point": axis_hits[0],
        "mask": mask,
        "statistics": statistics,
        "selection_mode": selection_mode,
    }


def sensitivity_scan(dataset: dict, cuts: tuple[float, ...]) -> list[dict]:
    results = []
    for cut in cuts:
        mask, stats = component_containing_axis(
            dataset["logq"],
            dataset["rboundary"],
            dataset["x"],
            dataset["height"],
            dataset["axis_point"],
            cut,
        )
        results.append(
            {
                "log10Q_cut": cut,
                "pixel_count": int(np.count_nonzero(mask)),
                "area_Mm2": stats["area_Mm2"],
                "bbox_Mm": stats["bbox_Mm"],
                "touches_lowest_Q_row": stats["touches_lowest_Q_row"],
            }
        )
    return results


def draw_panel(axis, dataset: dict, logq_cut: float, view_limits: tuple[float, ...]):
    image = axis.pcolormesh(
        dataset["x"],
        dataset["height"],
        dataset["logq"],
        shading="nearest",
        cmap="gray_r",
        vmin=LOGQ_MIN,
        vmax=LOGQ_MAX,
        rasterized=True,
        zorder=1,
    )
    region_layer = np.ma.masked_where(
        ~dataset["mask"], np.ones_like(dataset["mask"], dtype=float)
    )
    axis.pcolormesh(
        dataset["x"],
        dataset["height"],
        region_layer,
        shading="nearest",
        cmap=ListedColormap([REGION_FILL]),
        vmin=0.0,
        vmax=1.0,
        alpha=0.48,
        rasterized=True,
        zorder=2,
    )
    axis.contour(
        dataset["x"],
        dataset["height"],
        dataset["mask"].astype(float),
        levels=[0.5],
        colors=[REGION_EDGE],
        linewidths=1.25,
        zorder=2.4,
    )
    closed_contour_mask = dataset["rboundary"] == 11
    if dataset["selection_mode"] == "qsl-envelope":
        # Internal rboundary holes are geometrically filled in the selected
        # region, so draw only the outer purple connectivity boundary.
        closed_contour_mask = binary_fill_holes(closed_contour_mask)
    axis.contour(
        dataset["x"],
        dataset["height"],
        closed_contour_mask.astype(float),
        levels=[0.5],
        colors=[CLOSED_COLOR],
        linewidths=1.45,
        zorder=3.2,
    )
    axis.contour(
        dataset["x"],
        dataset["height"],
        dataset["logq"],
        levels=[logq_cut],
        colors=["white"],
        linewidths=3.5,
        zorder=2.6,
    )
    axis.contour(
        dataset["x"],
        dataset["height"],
        dataset["logq"],
        levels=[logq_cut],
        colors=[Q2_COLOR],
        linewidths=1.9,
        zorder=2.8,
    )
    axis.scatter(
        dataset["axis_point"][0],
        dataset["axis_point"][1],
        marker="X",
        s=115,
        facecolor=AXIS_COLOR,
        edgecolor="#101010",
        linewidth=0.9,
        zorder=5,
    )
    axis.set_title(
        f"t = {dataset['snapshot'] * TIME_MINUTES_PER_SNAPSHOT:.1f} min",
        fontsize=18,
        pad=8,
    )
    axis.set_xlim(view_limits[:2])
    axis.set_ylim(view_limits[2:])
    axis.set_aspect("equal")
    axis.tick_params(direction="in", top=True, right=True, labelsize=13)
    return image


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--snapshots", type=int, nargs="+", default=list(SNAPSHOTS))
    parser.add_argument("--logq-cut", type=float, default=4.0)
    parser.add_argument(
        "--selection-mode",
        choices=("qsl-envelope", "closed-component"),
        default="qsl-envelope",
        help=(
            "qsl-envelope selects the axis-containing low-Q basin and fills its internal "
            "layers; closed-component fills the complete rboundary==11 component that "
            "contains the magnetic axis"
        ),
    )
    parser.add_argument("--dpi", type=int, default=300)
    args = parser.parse_args()
    if not np.isfinite(args.logq_cut) or args.logq_cut <= LOGQ_MIN:
        parser.error(f"--logq-cut must be finite and larger than log10(2)={LOGQ_MIN}")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    datasets = []
    for snapshot in args.snapshots:
        source = args.input_dir / f"fixed_y_qsl_t{snapshot:04d}_yp0p000000_delta0.25.npz"
        if not source.is_file():
            raise FileNotFoundError(source)
        dataset = load_snapshot(source, snapshot, args.logq_cut, args.selection_mode)
        if args.selection_mode == "qsl-envelope":
            mask_name = f"y0_axis_envelope_logq{args.logq_cut:g}_t{snapshot:04d}.npz"
            mask_definition = (
                "hybrid boundary around the magnetic-axis side: the left/top boundary is the "
                "rboundary==11 connectivity boundary, the central/right boundary is the "
                "log10Q=cut ridge, and the lower opening is closed by the physical bottom; "
                "internal high-Q layers and enclosed holes are filled without crossing either "
                "outer boundary"
            )
        else:
            mask_name = f"y0_axis_closed_component_t{snapshot:04d}.npz"
            mask_definition = (
                "four-connected rboundary==11 component containing the magnetic-axis "
                "intersection; both field-line terminals reach the bottom boundary"
            )
        mask_output = args.output_dir / mask_name
        np.savez_compressed(
            mask_output,
            mask=dataset["mask"].astype(np.uint8),
            bottom_to_bottom_closed_within_mask=(
                dataset["mask"] & (dataset["rboundary"] == 11)
            ).astype(np.uint8),
            filled_internal_holes=(
                dataset["mask"] & (dataset["rboundary"] != 11)
            ).astype(np.uint8),
            x_Mm=dataset["x"].astype(np.float32),
            height_Mm=dataset["height"].astype(np.float32),
            magnetic_axis_intersection_Mm=dataset["axis_point"].astype(np.float32),
            log10Q_cut=np.asarray(args.logq_cut),
            connectivity=np.asarray(4),
            source=np.asarray(str(source.resolve())),
            selection_mode=np.asarray(args.selection_mode),
            definition=np.asarray(mask_definition),
        )
        dataset["mask_output"] = mask_output
        dataset["sensitivity"] = (
            sensitivity_scan(dataset, (1.5, 2.0, 2.5, 3.0, 4.0))
            if args.selection_mode == "qsl-envelope"
            else None
        )
        datasets.append(dataset)

    view_limits = (
        VIEW_LIMITS_MM
        if args.selection_mode == "qsl-envelope"
        else CLOSED_VIEW_LIMITS_MM
    )

    figure, axes = plt.subplots(
        1,
        len(datasets),
        figsize=(17.0, 7.8),
        sharex=True,
        sharey=True,
        constrained_layout=True,
    )
    axes = np.atleast_1d(axes)
    image = None
    for index, (axis, dataset) in enumerate(zip(axes, datasets)):
        image = draw_panel(axis, dataset, args.logq_cut, view_limits)
        axis.set_xlabel("Solar-X [Mm]", fontsize=15)
        if index == 0:
            axis.set_ylabel("Height above bottom boundary [Mm]", fontsize=15)

    region_label = (
        "Hybrid-boundary region on the magnetic-axis side"
        if args.selection_mode == "qsl-envelope"
        else "Axis-connected bottom-to-bottom closed field"
    )
    handles = [
        Patch(
            facecolor=REGION_FILL,
            edgecolor=REGION_EDGE,
            alpha=0.62,
            label=region_label,
        ),
        Line2D(
            [0], [0], color=Q2_COLOR, linewidth=2.3,
            label=rf"$\log_{{10}}Q={args.logq_cut:g}$",
        ),
        Line2D(
            [0], [0], color=CLOSED_COLOR, linewidth=2.0,
            label="Both terminals at bottom",
        ),
        Line2D(
            [0], [0], marker="X", linestyle="none", markersize=10,
            markerfacecolor=AXIS_COLOR, markeredgecolor="#101010",
            label="Magnetic-axis intersection",
        ),
    ]
    figure.legend(
        handles=handles,
        loc="outside lower center",
        ncol=4,
        frameon=False,
        fontsize=12,
        columnspacing=1.35,
    )
    colorbar = figure.colorbar(image, ax=axes, shrink=0.91, pad=0.02)
    colorbar.set_label(r"$\log_{10}Q$", fontsize=14)
    colorbar.ax.tick_params(labelsize=12)
    figure_title = (
        r"Hybrid closed-field/QSL boundary at $y=0$"
        if args.selection_mode == "qsl-envelope"
        else r"Bottom-to-bottom closed-field cross-section at $y=0$"
    )
    figure.suptitle(figure_title, fontsize=19)

    output_stem = (
        f"y0_axis_enclosed_region_logq{args.logq_cut:g}_t0010_t0095"
        if args.selection_mode == "qsl-envelope"
        else "y0_axis_closed_component_t0010_t0095"
    )
    output_base = args.output_dir / output_stem
    figure.savefig(output_base.with_suffix(".png"), dpi=args.dpi, bbox_inches="tight")
    plt.close(figure)

    if args.selection_mode == "qsl-envelope":
        definition = {
            "slice": "exact AMRVAC-normalized y=0 x-z plane",
            "eligible": f"rboundary==11 and finite log10Q<{args.logq_cut:g}",
            "selection": (
                "on each height row, fill from the leftmost point of the axis-connected "
                "rboundary==11 component to the rightmost point of the axis-containing "
                "log10Q-below-threshold component, then clip to that closed component"
            ),
            "boundary_construction": {
                "left_and_top": "rboundary==11 connectivity boundary (purple)",
                "central_and_right": f"log10Q={args.logq_cut:g} ridge (orange)",
                "bottom": "physical bottom boundary h=0",
            },
            "lower_closure": (
                "the component reaches the lowest cell-centered Q row; the physical h=0 bottom "
                "boundary closes the high-Q ridge's lower opening"
            ),
            "internal_Q_layers": (
                "included inside the outer row-wise envelope so the mask represents the complete "
                "geometric rope-plus-overlying-flux cross-section"
            ),
            "internal_holes": (
                "filled to make the geometric cross-section solid; the NPZ separately records "
                "bottom_to_bottom_closed_within_mask and filled_internal_holes"
            ),
            "interpretation": "candidate cross-section containing the flux rope and its overlying closed flux",
        }
    else:
        definition = {
            "slice": "exact AMRVAC-normalized y=0 x-z plane",
            "eligible": "rboundary==11",
            "selection": (
                "four-connected bottom-to-bottom closed-field component containing the "
                "magnetic-axis intersection"
            ),
            "interpretation": (
                "the complete closed-dome cross-section; this is broader than the flux that "
                "directly overlies the erupting current channel"
            ),
        }
    metadata = {
        "definition": definition,
        "selection_mode": args.selection_mode,
        "view_limits_Mm": list(view_limits),
        "output_png": str(output_base.with_suffix(".png").resolve()),
        "snapshots": [],
    }
    for dataset in datasets:
        metadata["snapshots"].append(
            {
                "snapshot": dataset["snapshot"],
                "time_minutes": dataset["snapshot"] * TIME_MINUTES_PER_SNAPSHOT,
                "source": str(dataset["source"].resolve()),
                "mask_npz": str(dataset["mask_output"].resolve()),
                "magnetic_axis_intersection_Mm": dataset["axis_point"].tolist(),
                "selected_component": dataset["statistics"],
                "threshold_sensitivity": dataset["sensitivity"],
            }
        )
    metadata_path = output_base.with_suffix(".json")
    metadata_path.write_text(json.dumps(metadata, indent=2, ensure_ascii=False) + "\n")
    print(output_base.with_suffix(".png").resolve())
    print(metadata_path.resolve())


if __name__ == "__main__":
    main()

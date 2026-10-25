#!/usr/bin/env python3
"""Build a time-resolved y=0 hybrid closed-field/QSL mask and GIF.

The default run covers every snapshot from 10 through 95.  Existing exact-y
FastQSL2 slice NPZs are reused; only missing slices are calculated, one at a
time, with :mod:`calculate_fixed_y_qsl_slices`.

The low-Q component is anchored by the real magnetic-axis intersection at
snapshot 10 in a forward pass and at snapshot 95 in a reverse pass.  Their
shared many-to-one component identifies the merge window.  Inside that window
a dynamic-programming path follows the same bottom-rooted
``log10(Q) == 4`` ridge through contemporary low-to-high Q crossings, rather
than switching to the much larger right-hand component.

For each height row, the final hybrid mask extends from the leftmost point of
the selected ``rboundary == 11`` dome to the tracked central high-Q ridge.
The region is filled solid, including enclosed holes.  The saved mask NPZ
explicitly distinguishes originally bottom-to-bottom-closed pixels from
geometrically filled holes and includes the mask boundary plus seed
coordinates ready for subsequent field-line tracing.

Example
-------
/opt/anaconda3/bin/python \
  /Users/xiaoxingyun/Desktop/solar_code/workflows/Simulation_Analysis/scripts/\
make_y0_hybrid_mask_evolution.py
"""

from __future__ import annotations

import argparse
import csv
import json
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Any

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from scipy.ndimage import binary_fill_holes, distance_transform_edt, label
from skimage.segmentation import watershed

from q4_seam_tracking import (
    SeamState,
    seam_from_mask,
    solid_mask_from_seam,
    track_q4_seam,
)


MAGNETIC_DIR = Path("/Volumes/Seagate/w2_jet_mhd/out_path1/nlfff")
SLICE_DIR = Path("/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/y0_qsl")
AXIS_FEATURE_DIR = Path(
    "/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/bottom_fastqsl2/"
    "overlying_flux_candidate"
)
OUTPUT_DIR = SLICE_DIR / "enclosed_region" / "evolution_logq4"
ENDPOINT_REFERENCE_DIR = SLICE_DIR / "enclosed_region"

START_SNAPSHOT = 10
STOP_SNAPSHOT = 95
FIXED_Y_CODE = 0.0
DELTA_INDEX = 0.25
LOGQ_CUT = 4.0
THREADS = 8
TIME_MINUTES_PER_SNAPSHOT = 0.71558
LENGTH_UNIT_MM = 10.0
BOTTOM_Z_CODE = 0.1
VIEW_LIMITS_MM = (-15.0, 30.0, 0.0, 30.0)
FOUR_CONNECTED = np.asarray(
    [[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.uint8
)

MASK_COLOR = "#43C7CA"
MASK_EDGE_COLOR = "#006B70"
CLOSED_COLOR = "#6A3D9A"
Q_COLOR = "#E69F00"
LEFT_MARKER_X_HEIGHT_MM = (0.0, 3.0)
RIGHT_MARKER_X_HEIGHT_MM = (30.0, 3.0)


@dataclass(frozen=True)
class SliceData:
    """Arrays required from one fixed-y FastQSL2 product."""

    snapshot: int
    source: Path
    logq: np.ndarray
    rboundary: np.ndarray
    x_mm: np.ndarray
    height_mm: np.ndarray
    axis_hits_mm: np.ndarray


@dataclass(frozen=True)
class ComponentChoice:
    """One four-connected low-Q component (kept for tracking diagnostics)."""

    component_id: int
    component_count: int
    mask: np.ndarray
    pixel_count: int
    overlap_pixels: int | None
    overlap_fraction_previous: float | None
    overlap_fraction_current: float | None
    overlap_iou: float | None


@dataclass(frozen=True)
class WatershedChoice:
    """Marker-controlled partition and low-Q pixels on the selected side."""

    selected_partition: np.ndarray
    selected_lowq: np.ndarray
    native_closed_component: np.ndarray
    filled_closed_component: np.ndarray
    left_marker_index: tuple[int, int]
    right_marker_index: tuple[int, int]
    left_marker_grid_Mm: tuple[float, float]
    right_marker_grid_Mm: tuple[float, float]


@dataclass(frozen=True)
class FrameMask:
    """Final mask and reproducibility information for one snapshot."""

    mask: np.ndarray
    selected_lowq: np.ndarray
    bottom_to_bottom: np.ndarray
    filled_internal_holes: np.ndarray
    closed_component: np.ndarray
    selection_mode: str
    seam_state: SeamState | None
    statistics: dict[str, Any]


def fixed_y_tag(value: float) -> str:
    """Match the six-decimal y tag used by calculate_fixed_y_qsl_slices.py."""

    return f"{value:+.6f}".replace("+", "p").replace("-", "m").replace(".", "p")


def slice_path(slice_dir: Path, snapshot: int, fixed_y: float, delta: float) -> Path:
    return slice_dir / (
        f"fixed_y_qsl_t{snapshot:04d}_y{fixed_y_tag(fixed_y)}_delta{delta:.2f}.npz"
    )


def grid_vectors(
    x_mm: np.ndarray, height_mm: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Return separable 1-D coordinates from the two seed-coordinate grids."""

    x_array = np.asarray(x_mm, dtype=float)
    height_array = np.asarray(height_mm, dtype=float)
    x = x_array[0, :] if x_array.ndim == 2 else x_array
    height = height_array[:, 0] if height_array.ndim == 2 else height_array
    if x.ndim != 1 or height.ndim != 1:
        raise ValueError("Expected separable fixed-y x and height coordinates")
    if x.size < 2 or height.size < 2:
        raise ValueError("The fixed-y grid must contain at least two points per axis")
    if not np.all(np.diff(x) > 0.0) or not np.all(np.diff(height) > 0.0):
        raise ValueError("Fixed-y x and height coordinates must increase monotonically")
    return x, height


def load_slice(path: Path, snapshot: int) -> SliceData:
    """Load and validate one cached or newly calculated FastQSL2 slice."""

    with np.load(path, allow_pickle=False) as archive:
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
        x_mm, height_mm = grid_vectors(archive["x_mm"], archive["height_mm"])
        axis_hits_mm = np.asarray(
            archive["magnetic_axis_intersections_Mm"], dtype=float
        )

    expected_shape = (height_mm.size, x_mm.size)
    if logq.shape != expected_shape or rboundary.shape != expected_shape:
        raise ValueError(
            f"Inconsistent shapes in {path}: logq={logq.shape}, "
            f"rboundary={rboundary.shape}, coordinates={expected_shape}"
        )
    if axis_hits_mm.size == 0:
        axis_hits_mm = np.empty((0, 2), dtype=float)
    elif axis_hits_mm.ndim != 2 or axis_hits_mm.shape[1] != 2:
        raise ValueError(
            f"Unexpected magnetic-axis intersection shape in {path}: "
            f"{axis_hits_mm.shape}"
        )
    return SliceData(
        snapshot=snapshot,
        source=path.resolve(),
        logq=logq,
        rboundary=rboundary,
        x_mm=x_mm,
        height_mm=height_mm,
        axis_hits_mm=axis_hits_mm,
    )


def ensure_slice(
    snapshot: int,
    *,
    magnetic_dir: Path,
    slice_dir: Path,
    axis_feature_dir: Path,
    fixed_y: float,
    delta: float,
    threads: int,
) -> tuple[Path, str]:
    """Reuse a valid cached slice or calculate the missing snapshot once."""

    path = slice_path(slice_dir, snapshot, fixed_y, delta)
    if path.is_file():
        # Validate now so a corrupt cache is never silently accepted.
        load_slice(path, snapshot)
        return path, "cached"

    magnetic_path = magnetic_dir / f"nlfff_{snapshot:04d}.binary"
    feature_path = axis_feature_dir / f"axis_features_factor2_t{snapshot:04d}.npz"
    if not magnetic_path.is_file():
        raise FileNotFoundError(
            f"Missing both cached slice and magnetic binary for t={snapshot}: "
            f"{magnetic_path}"
        )

    # Import lazily so --help and cached-only validation do not initialize
    # FastQSL2 or require its executable.
    from calculate_fixed_y_qsl_slices import calculate_slice

    print(f"Calculating missing fixed-y QSL slice t={snapshot:04d} ...", flush=True)
    calculate_slice(
        snapshot,
        magnetic_path,
        feature_path,
        slice_dir,
        fixed_y,
        delta,
        threads,
    )
    if not path.is_file():
        raise RuntimeError(f"FastQSL2 calculation did not create the expected file: {path}")
    load_slice(path, snapshot)
    return path, "calculated"


def component_labels(data: SliceData, logq_cut: float) -> tuple[np.ndarray, int]:
    eligible = (
        (data.rboundary == 11)
        & np.isfinite(data.logq)
        & (data.logq < logq_cut)
    )
    return label(eligible, structure=FOUR_CONNECTED)


def axis_grid_index(data: SliceData) -> tuple[int, int, np.ndarray]:
    """Map the unique real magnetic-axis intersection to its nearest grid point."""

    if data.axis_hits_mm.shape != (1, 2):
        raise ValueError(
            f"Snapshot {data.snapshot} must have exactly one real magnetic-axis "
            f"intersection; found {data.axis_hits_mm.shape[0]} in {data.source}"
        )
    point = data.axis_hits_mm[0]
    ix = int(np.argmin(np.abs(data.x_mm - point[0])))
    iz = int(np.argmin(np.abs(data.height_mm - point[1])))
    return iz, ix, point


def choice_from_id(
    components: np.ndarray,
    component_count: int,
    component_id: int,
    *,
    previous_mask: np.ndarray | None,
) -> ComponentChoice:
    if not 1 <= component_id <= component_count:
        raise ValueError(f"Invalid component id {component_id} of {component_count}")
    mask = components == component_id
    pixel_count = int(np.count_nonzero(mask))
    if previous_mask is None:
        return ComponentChoice(
            component_id=component_id,
            component_count=component_count,
            mask=mask,
            pixel_count=pixel_count,
            overlap_pixels=None,
            overlap_fraction_previous=None,
            overlap_fraction_current=None,
            overlap_iou=None,
        )

    overlap = int(np.count_nonzero(mask & previous_mask))
    previous_count = int(np.count_nonzero(previous_mask))
    union = pixel_count + previous_count - overlap
    return ComponentChoice(
        component_id=component_id,
        component_count=component_count,
        mask=mask,
        pixel_count=pixel_count,
        overlap_pixels=overlap,
        overlap_fraction_previous=overlap / previous_count,
        overlap_fraction_current=overlap / pixel_count,
        overlap_iou=overlap / union,
    )


def choose_axis_component(data: SliceData, logq_cut: float) -> ComponentChoice:
    components, component_count = component_labels(data, logq_cut)
    iz, ix, _ = axis_grid_index(data)
    component_id = int(components[iz, ix])
    if component_id == 0:
        raise ValueError(
            f"The real magnetic-axis grid point at t={data.snapshot} is outside "
            f"rboundary==11 and finite log10Q<{logq_cut:g}; "
            f"rboundary={data.rboundary[iz, ix]}, log10Q={data.logq[iz, ix]}"
        )
    return choice_from_id(
        components, component_count, component_id, previous_mask=None
    )


def choose_maximum_overlap_component(
    data: SliceData,
    logq_cut: float,
    previous_mask: np.ndarray,
) -> ComponentChoice:
    """Select the component with maximum raw overlap, then IoU, then size."""

    components, component_count = component_labels(data, logq_cut)
    if component_count == 0:
        raise RuntimeError(f"No eligible low-Q closed component at t={data.snapshot}")
    overlap_counts = np.bincount(
        components[previous_mask].ravel(), minlength=component_count + 1
    )
    overlap_counts[0] = 0
    maximum_overlap = int(overlap_counts.max())
    if maximum_overlap == 0:
        raise RuntimeError(
            f"No four-connected low-Q component at t={data.snapshot} overlaps "
            "the previously tracked component"
        )

    candidate_ids = np.flatnonzero(overlap_counts == maximum_overlap)
    choices = [
        choice_from_id(
            components,
            component_count,
            int(component_id),
            previous_mask=previous_mask,
        )
        for component_id in candidate_ids
    ]
    # Ties in raw overlap are resolved deterministically without changing the
    # primary maximum-overlap definition.
    return max(
        choices,
        key=lambda item: (
            float(item.overlap_iou),
            item.pixel_count,
            -item.component_id,
        ),
    )


def forward_track(
    sources: dict[int, Path], snapshots: list[int], logq_cut: float
) -> dict[int, ComponentChoice]:
    choices: dict[int, ComponentChoice] = {}
    first = load_slice(sources[snapshots[0]], snapshots[0])
    choices[snapshots[0]] = choose_axis_component(first, logq_cut)
    previous_mask = choices[snapshots[0]].mask
    for snapshot in snapshots[1:]:
        data = load_slice(sources[snapshot], snapshot)
        choice = choose_maximum_overlap_component(data, logq_cut, previous_mask)
        choices[snapshot] = choice
        previous_mask = choice.mask
    return choices


def reverse_track(
    sources: dict[int, Path], snapshots: list[int], logq_cut: float
) -> dict[int, ComponentChoice]:
    """Independently track backward from the real final-axis component."""

    reverse: dict[int, ComponentChoice] = {}
    last_snapshot = snapshots[-1]
    last = load_slice(sources[last_snapshot], last_snapshot)
    reverse[last_snapshot] = choose_axis_component(last, logq_cut)
    next_mask = reverse[last_snapshot].mask
    for snapshot in reversed(snapshots[:-1]):
        data = load_slice(sources[snapshot], snapshot)
        choice = choose_maximum_overlap_component(data, logq_cut, next_mask)
        reverse[snapshot] = choice
        next_mask = choice.mask
    return reverse


def select_closed_component(
    data: SliceData, lowq_component: np.ndarray
) -> tuple[np.ndarray, int, int]:
    """Find the unique rboundary==11 component containing the tracked basin."""

    closed_labels, closed_count = label(
        data.rboundary == 11, structure=FOUR_CONNECTED
    )
    ids, counts = np.unique(closed_labels[lowq_component], return_counts=True)
    nonzero = [(int(i), int(n)) for i, n in zip(ids, counts) if i != 0]
    if len(nonzero) != 1:
        raise RuntimeError(
            f"Tracked low-Q component at t={data.snapshot} belongs to "
            f"{len(nonzero)} nonzero closed components, expected exactly one"
        )
    component_id = nonzero[0][0]
    closed_component = closed_labels == component_id
    if not np.all(closed_component[lowq_component]):
        raise RuntimeError(
            f"Tracked low-Q component is not fully inside its closed component at "
            f"t={data.snapshot}"
        )
    return closed_component, component_id, int(closed_count)


def build_hybrid_mask(
    data: SliceData, lowq_component: np.ndarray, logq_cut: float
) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Construct the solid row-wise hybrid envelope used by the endpoint masks."""

    closed_component, closed_id, closed_count = select_closed_component(
        data, lowq_component
    )
    geometric_span = np.zeros_like(lowq_component, dtype=bool)
    selected_rows = np.flatnonzero(np.any(lowq_component, axis=1))
    for row in selected_rows:
        lowq_columns = np.flatnonzero(lowq_component[row])
        closed_columns = np.flatnonzero(closed_component[row])
        if closed_columns.size == 0:
            raise RuntimeError(
                f"No axis-connected closed pixels on selected row {row} at "
                f"t={data.snapshot}"
            )
        geometric_span[
            row, int(closed_columns.min()) : int(lowq_columns.max()) + 1
        ] = True

    closed_only = geometric_span & closed_component
    # This precisely preserves the approved endpoint convention: the outer
    # row-wise bounds do not move, while every fully enclosed internal hole is
    # made solid, including rboundary != 11 holes at t=95.
    mask = binary_fill_holes(closed_only) & geometric_span
    bottom_to_bottom = mask & (data.rboundary == 11)
    filled_internal_holes = mask & (data.rboundary != 11)
    if not np.all(lowq_component <= mask):
        raise RuntimeError(f"Hybrid mask lost tracked low-Q pixels at t={data.snapshot}")

    rows, columns = np.where(mask)
    if rows.size == 0:
        raise RuntimeError(f"Hybrid mask is empty at t={data.snapshot}")
    dx = float(np.median(np.diff(data.x_mm)))
    dh = float(np.median(np.diff(data.height_mm)))
    left_boundary_matches = 0
    right_highq_matches = 0
    right_highq_or_closed_join_matches = 0
    for row in np.flatnonzero(np.any(mask, axis=1)):
        mask_columns = np.flatnonzero(mask[row])
        left = int(mask_columns.min())
        right = int(mask_columns.max())
        if left == 0 or not closed_component[row, left - 1]:
            left_boundary_matches += 1
        if right + 1 < mask.shape[1]:
            next_is_highq = bool(
                np.isfinite(data.logq[row, right + 1])
                and data.logq[row, right + 1] >= logq_cut
            )
            next_is_closed_join = not bool(closed_component[row, right + 1])
            if next_is_highq:
                right_highq_matches += 1
            if next_is_highq or next_is_closed_join:
                right_highq_or_closed_join_matches += 1

    _, final_component_count = label(mask, structure=FOUR_CONNECTED)
    residual_hole_count = int(np.count_nonzero(binary_fill_holes(mask) & ~mask))

    statistics: dict[str, Any] = {
        "closed_component_id": closed_id,
        "closed_component_count": closed_count,
        "lowq_component_pixel_count": int(np.count_nonzero(lowq_component)),
        "geometric_span_pixel_count": int(np.count_nonzero(geometric_span)),
        "bottom_to_bottom_closed_pixel_count": int(
            np.count_nonzero(bottom_to_bottom)
        ),
        "filled_internal_hole_pixel_count": int(
            np.count_nonzero(filled_internal_holes)
        ),
        "mask_pixel_count": int(rows.size),
        "area_Mm2": float(rows.size * abs(dx * dh)),
        "bbox_x_min_Mm": float(data.x_mm[columns].min()),
        "bbox_x_max_Mm": float(data.x_mm[columns].max()),
        "bbox_height_min_Mm": float(data.height_mm[rows].min()),
        "bbox_height_max_Mm": float(data.height_mm[rows].max()),
        "selected_height_row_count": int(np.unique(rows).size),
        "left_rows_matching_closed_boundary": int(left_boundary_matches),
        "right_rows_matching_high_Q_ridge": int(right_highq_matches),
        "right_rows_matching_high_Q_or_closed_join": int(
            right_highq_or_closed_join_matches
        ),
        "final_mask_component_count": int(final_component_count),
        "residual_internal_hole_pixel_count": residual_hole_count,
        "centroid_x_Mm": float(np.mean(data.x_mm[columns])),
        "centroid_height_Mm": float(np.mean(data.height_mm[rows])),
        "touches_lowest_Q_row": bool(np.any(mask[0, :])),
        "touches_top_Q_row": bool(np.any(mask[-1, :])),
        "touches_left_Q_edge": bool(np.any(mask[:, 0])),
        "touches_right_Q_edge": bool(np.any(mask[:, -1])),
    }
    return mask, bottom_to_bottom, filled_internal_holes, statistics


def marker_grid_index(
    data: SliceData, marker_x_height_mm: tuple[float, float]
) -> tuple[int, int]:
    """Map a physical ``(x, height)`` marker to its nearest grid cell."""

    marker_x, marker_height = marker_x_height_mm
    ix = int(np.argmin(np.abs(data.x_mm - marker_x)))
    iz = int(np.argmin(np.abs(data.height_mm - marker_height)))
    return iz, ix


def choose_watershed_side(
    data: SliceData,
    logq_cut: float,
    left_marker_x_height_mm: tuple[float, float],
    right_marker_x_height_mm: tuple[float, float],
) -> WatershedChoice:
    """Partition one closed dome along its high-Q ridge using two markers."""

    left_index = marker_grid_index(data, left_marker_x_height_mm)
    right_index = marker_grid_index(data, right_marker_x_height_mm)
    if left_index == right_index:
        raise ValueError("The two watershed markers map to the same grid cell")

    closed_labels, _ = label(data.rboundary == 11, structure=FOUR_CONNECTED)
    left_closed_id = int(closed_labels[left_index])
    right_closed_id = int(closed_labels[right_index])
    if left_closed_id == 0 or right_closed_id == 0:
        raise ValueError(
            f"Watershed marker outside rboundary==11 at t={data.snapshot}: "
            f"left_id={left_closed_id}, right_id={right_closed_id}"
        )
    if left_closed_id != right_closed_id:
        raise ValueError(
            f"Watershed markers lie in different closed components at "
            f"t={data.snapshot}: {left_closed_id} != {right_closed_id}"
        )

    native_closed = closed_labels == left_closed_id
    filled_closed = binary_fill_holes(native_closed)
    finite = np.isfinite(data.logq)
    finite_values = data.logq[finite & filled_closed]
    if finite_values.size == 0:
        raise RuntimeError(
            f"No finite log10Q values in the closed dome at t={data.snapshot}"
        )
    floor = float(np.min(finite_values))
    ceiling = float(max(np.max(finite_values), logq_cut + 1.0))
    elevation = np.nan_to_num(
        data.logq,
        nan=ceiling,
        posinf=ceiling,
        neginf=floor,
    )
    elevation = np.clip(elevation, floor, ceiling)

    markers = np.zeros(data.logq.shape, dtype=np.int16)
    markers[left_index] = 1
    markers[right_index] = 2
    partition = watershed(
        elevation,
        markers=markers,
        connectivity=1,
        mask=filled_closed,
    )
    selected_partition = partition == 1
    if not selected_partition[left_index] or selected_partition[right_index]:
        raise RuntimeError(f"Invalid watershed marker assignment at t={data.snapshot}")

    selected_lowq = (
        selected_partition
        & native_closed
        & finite
        & (data.logq < logq_cut)
    )
    if not selected_lowq[left_index]:
        raise ValueError(
            f"Left marker is not a finite low-Q closed pixel at t={data.snapshot}: "
            f"rboundary={data.rboundary[left_index]}, "
            f"log10Q={data.logq[left_index]}"
        )

    return WatershedChoice(
        selected_partition=selected_partition,
        selected_lowq=selected_lowq,
        native_closed_component=native_closed,
        filled_closed_component=filled_closed,
        left_marker_index=left_index,
        right_marker_index=right_index,
        left_marker_grid_Mm=(
            float(data.x_mm[left_index[1]]),
            float(data.height_mm[left_index[0]]),
        ),
        right_marker_grid_Mm=(
            float(data.x_mm[right_index[1]]),
            float(data.height_mm[right_index[0]]),
        ),
    )


def build_watershed_hybrid_mask(
    data: SliceData,
    choice: WatershedChoice,
    logq_cut: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, dict[str, Any]]:
    """Build one connected, solid mask from the watershed-selected Q side."""

    lowq_region = choice.selected_lowq
    geometric_span = np.zeros_like(lowq_region, dtype=bool)
    selected_rows = np.flatnonzero(np.any(lowq_region, axis=1))
    for row in selected_rows:
        lowq_columns = np.flatnonzero(lowq_region[row])
        closed_columns = np.flatnonzero(choice.filled_closed_component[row])
        if closed_columns.size == 0:
            raise RuntimeError(
                f"No closed-dome pixels on selected row {row} at t={data.snapshot}"
            )
        geometric_span[
            row, int(closed_columns.min()) : int(lowq_columns.max()) + 1
        ] = True

    closed_only = geometric_span & choice.native_closed_component
    candidate = binary_fill_holes(closed_only) & geometric_span
    candidate_labels, candidate_count = label(
        candidate, structure=FOUR_CONNECTED
    )
    selected_id = int(candidate_labels[choice.left_marker_index])
    if selected_id == 0:
        raise RuntimeError(
            f"The solid-envelope candidate lost the left marker at t={data.snapshot}"
        )
    mask = candidate_labels == selected_id
    mask = (
        binary_fill_holes(mask)
        & geometric_span
        & choice.filled_closed_component
    )
    _, final_component_count = label(mask, structure=FOUR_CONNECTED)
    if final_component_count != 1:
        raise RuntimeError(
            f"Expected one connected final mask at t={data.snapshot}; "
            f"found {final_component_count}"
        )

    selected_lowq_within_mask = lowq_region & mask
    bottom_to_bottom = mask & choice.native_closed_component
    filled_internal_holes = mask & ~choice.native_closed_component
    rows, columns = np.where(mask)
    if rows.size == 0:
        raise RuntimeError(f"Hybrid mask is empty at t={data.snapshot}")

    dx = float(np.median(np.diff(data.x_mm)))
    dh = float(np.median(np.diff(data.height_mm)))
    left_boundary_matches = 0
    right_highq_matches = 0
    right_partition_matches = 0
    right_highq_or_closed_join_matches = 0
    for row in np.flatnonzero(np.any(mask, axis=1)):
        mask_columns = np.flatnonzero(mask[row])
        left = int(mask_columns.min())
        right = int(mask_columns.max())
        if left == 0 or not choice.filled_closed_component[row, left - 1]:
            left_boundary_matches += 1
        if right + 1 < mask.shape[1]:
            next_is_highq = bool(
                np.isfinite(data.logq[row, right + 1])
                and data.logq[row, right + 1] >= logq_cut
            )
            next_is_other_partition = bool(
                choice.filled_closed_component[row, right + 1]
                and not choice.selected_partition[row, right + 1]
            )
            next_is_closed_join = not bool(
                choice.filled_closed_component[row, right + 1]
            )
            if next_is_highq:
                right_highq_matches += 1
            if next_is_other_partition:
                right_partition_matches += 1
            if next_is_highq or next_is_closed_join:
                right_highq_or_closed_join_matches += 1

    closed_labels, closed_count = label(
        data.rboundary == 11, structure=FOUR_CONNECTED
    )
    residual_hole_count = int(np.count_nonzero(binary_fill_holes(mask) & ~mask))
    statistics: dict[str, Any] = {
        "closed_component_id": int(closed_labels[choice.left_marker_index]),
        "closed_component_count": int(closed_count),
        "watershed_selected_partition_pixel_count": int(
            np.count_nonzero(choice.selected_partition)
        ),
        "watershed_other_partition_pixel_count": int(
            np.count_nonzero(
                choice.filled_closed_component & ~choice.selected_partition
            )
        ),
        "selected_lowq_pixel_count": int(np.count_nonzero(lowq_region)),
        "selected_lowq_within_final_mask_pixel_count": int(
            np.count_nonzero(selected_lowq_within_mask)
        ),
        "discarded_detached_lowq_pixel_count": int(
            np.count_nonzero(lowq_region & ~mask)
        ),
        "preselection_component_count": int(candidate_count),
        "geometric_span_pixel_count": int(np.count_nonzero(geometric_span)),
        "bottom_to_bottom_closed_pixel_count": int(
            np.count_nonzero(bottom_to_bottom)
        ),
        "filled_internal_hole_pixel_count": int(
            np.count_nonzero(filled_internal_holes)
        ),
        "mask_pixel_count": int(rows.size),
        "area_Mm2": float(rows.size * abs(dx * dh)),
        "bbox_x_min_Mm": float(data.x_mm[columns].min()),
        "bbox_x_max_Mm": float(data.x_mm[columns].max()),
        "bbox_height_min_Mm": float(data.height_mm[rows].min()),
        "bbox_height_max_Mm": float(data.height_mm[rows].max()),
        "selected_height_row_count": int(np.unique(rows).size),
        "left_rows_matching_closed_boundary": int(left_boundary_matches),
        "right_rows_matching_high_Q_ridge": int(right_highq_matches),
        "right_rows_matching_watershed_partition": int(right_partition_matches),
        "right_rows_matching_high_Q_or_closed_join": int(
            right_highq_or_closed_join_matches
        ),
        "final_mask_component_count": int(final_component_count),
        "residual_internal_hole_pixel_count": residual_hole_count,
        "centroid_x_Mm": float(np.mean(data.x_mm[columns])),
        "centroid_height_Mm": float(np.mean(data.height_mm[rows])),
        "left_marker_grid_x_Mm": choice.left_marker_grid_Mm[0],
        "left_marker_grid_height_Mm": choice.left_marker_grid_Mm[1],
        "right_marker_grid_x_Mm": choice.right_marker_grid_Mm[0],
        "right_marker_grid_height_Mm": choice.right_marker_grid_Mm[1],
        "touches_lowest_Q_row": bool(np.any(mask[0, :])),
        "touches_top_Q_row": bool(np.any(mask[-1, :])),
        "touches_left_Q_edge": bool(np.any(mask[:, 0])),
        "touches_right_Q_edge": bool(np.any(mask[:, -1])),
    }
    return (
        mask,
        selected_lowq_within_mask,
        bottom_to_bottom,
        filled_internal_holes,
        statistics,
    )


def find_shared_merge_window(
    forward: dict[int, ComponentChoice],
    reverse: dict[int, ComponentChoice],
    snapshots: list[int],
) -> list[int]:
    """Find the interior run where endpoint-anchored tracks share one component.

    Before the run the forward track identifies the selected side, after the
    run the reverse track does.  Equality in the middle is the many-to-one
    merge that makes component identity alone ambiguous.
    """

    equal = [
        np.array_equal(forward[snapshot].mask, reverse[snapshot].mask)
        for snapshot in snapshots
    ]
    runs: list[list[int]] = []
    current: list[int] = []
    for snapshot, is_equal in zip(snapshots, equal):
        if is_equal:
            current.append(snapshot)
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)

    interior = [
        run
        for run in runs
        if run[0] != snapshots[0] and run[-1] != snapshots[-1]
    ]
    if len(interior) != 1:
        raise RuntimeError(
            "Expected exactly one interior forward/reverse component-merge "
            f"window; found {interior}"
        )
    window = interior[0]
    indices = [snapshots.index(snapshot) for snapshot in window]
    expected = list(range(indices[0], indices[-1] + 1))
    if indices != expected:
        raise RuntimeError(f"Non-contiguous component-merge window: {window}")
    return window


def right_boundary_indices(mask: np.ndarray) -> np.ndarray:
    """Return ``(iz, ix)`` for the rightmost selected cell on every mask row."""

    rows = np.flatnonzero(np.any(mask, axis=1))
    return np.asarray(
        [(int(row), int(np.flatnonzero(mask[row]).max())) for row in rows],
        dtype=np.int32,
    )


def summarize_final_mask(
    data: SliceData,
    mask: np.ndarray,
    selected_lowq: np.ndarray,
    *,
    logq_cut: float,
    selection_mode: str,
    component_choice: ComponentChoice | None,
    seam_state: SeamState | None,
    left_marker_x_height_mm: tuple[float, float],
) -> dict[str, Any]:
    """Create one uniform QA record for component- and seam-selected masks."""

    mask = np.asarray(mask, dtype=bool)
    selected_lowq = np.asarray(selected_lowq, dtype=bool)
    if mask.shape != data.logq.shape or selected_lowq.shape != mask.shape:
        raise ValueError("Mask and selected low-Q arrays must match the Q grid")
    rows, columns = np.where(mask)
    if rows.size == 0:
        raise RuntimeError(f"Hybrid mask is empty at t={data.snapshot}")

    left_marker_index = marker_grid_index(data, left_marker_x_height_mm)
    closed_labels, closed_count = label(
        data.rboundary == 11, structure=FOUR_CONNECTED
    )
    closed_id = int(closed_labels[left_marker_index])
    if closed_id == 0:
        raise RuntimeError(f"Left marker is outside the closed dome at t={data.snapshot}")
    native_closed = closed_labels == closed_id
    filled_closed = binary_fill_holes(native_closed)
    if not mask[left_marker_index]:
        raise RuntimeError(f"Final mask does not contain the left marker at t={data.snapshot}")

    boundaries = right_boundary_indices(mask)
    highq = np.isfinite(data.logq) & (data.logq >= logq_cut)
    dx = float(np.median(np.diff(data.x_mm)))
    dh = float(np.median(np.diff(data.height_mm)))
    highq_distance_mm = distance_transform_edt(
        ~highq, sampling=(abs(dh), abs(dx))
    )
    boundary_distances = highq_distance_mm[
        boundaries[:, 0], boundaries[:, 1]
    ]

    left_matches = 0
    right_highq_matches = 0
    right_highq_or_closed_join_matches = 0
    for row, right in boundaries:
        row = int(row)
        right = int(right)
        left = int(np.flatnonzero(mask[row]).min())
        if left == 0 or not filled_closed[row, left - 1]:
            left_matches += 1
        if right + 1 < mask.shape[1]:
            next_is_highq = bool(highq[row, right + 1])
            next_is_closed_join = not bool(filled_closed[row, right + 1])
            if next_is_highq:
                right_highq_matches += 1
            if next_is_highq or next_is_closed_join:
                right_highq_or_closed_join_matches += 1

    _, eligible_component_count = label(
        native_closed & np.isfinite(data.logq) & (data.logq < logq_cut),
        structure=FOUR_CONNECTED,
    )
    _, final_component_count = label(mask, structure=FOUR_CONNECTED)
    residual_hole_count = int(np.count_nonzero(binary_fill_holes(mask) & ~mask))
    if final_component_count != 1 or residual_hole_count != 0:
        raise RuntimeError(
            f"Final mask QA failed at t={data.snapshot}: "
            f"components={final_component_count}, holes={residual_hole_count}"
        )

    seam_audit = seam_state.audit if seam_state is not None else {}
    return {
        "selection_mode": selection_mode,
        "source_component_id": (
            component_choice.component_id if component_choice is not None else None
        ),
        "source_component_pixel_count": (
            component_choice.pixel_count if component_choice is not None else None
        ),
        "eligible_component_count": int(eligible_component_count),
        "closed_component_id": closed_id,
        "closed_component_count": int(closed_count),
        "selected_lowq_pixel_count": int(np.count_nonzero(selected_lowq)),
        "bottom_to_bottom_closed_pixel_count": int(
            np.count_nonzero(mask & (data.rboundary == 11))
        ),
        "filled_internal_hole_pixel_count": int(
            np.count_nonzero(mask & (data.rboundary != 11))
        ),
        "mask_pixel_count": int(rows.size),
        "area_Mm2": float(rows.size * abs(dx * dh)),
        "bbox_x_min_Mm": float(data.x_mm[columns].min()),
        "bbox_x_max_Mm": float(data.x_mm[columns].max()),
        "bbox_height_min_Mm": float(data.height_mm[rows].min()),
        "bbox_height_max_Mm": float(data.height_mm[rows].max()),
        "selected_height_row_count": int(np.unique(rows).size),
        "left_rows_matching_closed_boundary": int(left_matches),
        "right_rows_matching_high_Q_ridge": int(right_highq_matches),
        "right_rows_matching_high_Q_or_closed_join": int(
            right_highq_or_closed_join_matches
        ),
        "right_boundary_median_distance_to_high_Q_Mm": float(
            np.median(boundary_distances)
        ),
        "right_boundary_max_distance_to_high_Q_Mm": float(
            np.max(boundary_distances)
        ),
        "final_mask_component_count": int(final_component_count),
        "residual_internal_hole_pixel_count": residual_hole_count,
        "centroid_x_Mm": float(np.mean(data.x_mm[columns])),
        "centroid_height_Mm": float(np.mean(data.height_mm[rows])),
        "seam_join_row": seam_state.join_row if seam_state is not None else None,
        "seam_join_height_Mm": seam_audit.get("join_height_Mm"),
        "seam_join_side": seam_state.join_side if seam_state is not None else None,
        "seam_path_node_count": seam_audit.get("path_node_count"),
        "seam_maximum_gap_rows": seam_audit.get("maximum_gap_rows"),
        "seam_maximum_gap_Mm": seam_audit.get("maximum_gap_Mm"),
        "seam_objective": seam_state.objective if seam_state is not None else None,
        "seam_path_mean_cost": (
            seam_state.path_mean_cost if seam_state is not None else None
        ),
        "touches_lowest_Q_row": bool(np.any(mask[0, :])),
        "touches_top_Q_row": bool(np.any(mask[-1, :])),
        "touches_left_Q_edge": bool(np.any(mask[:, 0])),
        "touches_right_Q_edge": bool(np.any(mask[:, -1])),
    }


def endpoint_reference_path(reference_dir: Path, snapshot: int, cut: float) -> Path:
    return reference_dir / f"y0_axis_envelope_logq{cut:g}_t{snapshot:04d}.npz"


def verify_endpoint_mask(
    mask: np.ndarray,
    *,
    snapshot: int,
    reference_dir: Path,
    logq_cut: float,
) -> dict[str, Any]:
    """Require pixel-for-pixel equality with the approved endpoint mask."""

    path = endpoint_reference_path(reference_dir, snapshot, logq_cut)
    if not path.is_file():
        raise FileNotFoundError(f"Required endpoint reference mask is missing: {path}")
    with np.load(path, allow_pickle=False) as archive:
        if "mask" not in archive.files:
            raise KeyError(f"Endpoint reference has no 'mask' array: {path}")
        reference = np.asarray(archive["mask"], dtype=bool)
    if reference.shape != mask.shape:
        raise ValueError(
            f"Endpoint reference shape mismatch at t={snapshot}: "
            f"{reference.shape} != {mask.shape}"
        )
    differing = int(np.count_nonzero(reference ^ mask))
    if differing:
        raise RuntimeError(
            f"Endpoint mask regression failed at t={snapshot}: {differing} "
            f"pixels differ from {path}"
        )
    return {
        "reference": str(path.resolve()),
        "pixel_equal": True,
        "differing_pixel_count": 0,
    }


def seed_coordinates(
    mask: np.ndarray,
    x_mm: np.ndarray,
    height_mm: np.ndarray,
    fixed_y: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return selected centers in plotting and AMRVAC coordinates."""

    rows, columns = np.where(mask)
    x_selected = x_mm[columns]
    height_selected = height_mm[rows]
    x_height = np.column_stack((x_selected, height_selected)).astype(float)
    xyz = np.column_stack(
        (
            x_selected / LENGTH_UNIT_MM,
            np.full(rows.size, fixed_y, dtype=float),
            height_selected / LENGTH_UNIT_MM + BOTTOM_Z_CODE,
        )
    )
    return x_height, xyz


def save_mask_npz(
    path: Path,
    *,
    data: SliceData,
    mask: np.ndarray,
    selected_lowq_region: np.ndarray,
    selection_mode: str,
    seam_state: SeamState | None,
    bottom_to_bottom: np.ndarray,
    filled_internal_holes: np.ndarray,
    fixed_y: float,
    delta: float,
    logq_cut: float,
) -> None:
    selected_x_height, selected_xyz = seed_coordinates(
        mask, data.x_mm, data.height_mm, fixed_y
    )
    selected_rows, selected_columns = np.where(mask)
    boundary_indices = right_boundary_indices(mask)
    boundary_x_height = np.column_stack(
        (
            data.x_mm[boundary_indices[:, 1]],
            data.height_mm[boundary_indices[:, 0]],
        )
    ).astype(np.float64)
    seam_indices = (
        np.empty((0, 2), dtype=np.int32)
        if seam_state is None
        else np.asarray(seam_state.path_nodes_iz_ix, dtype=np.int32)
    )
    seam_x_height = (
        np.empty((0, 2), dtype=np.float64)
        if seam_indices.size == 0
        else np.column_stack(
            (
                data.x_mm[seam_indices[:, 1]],
                data.height_mm[seam_indices[:, 0]],
            )
        ).astype(np.float64)
    )
    np.savez_compressed(
        path,
        mask=mask.astype(np.uint8),
        selected_lowq_region=selected_lowq_region.astype(np.uint8),
        tracked_lowq_component=selected_lowq_region.astype(np.uint8),
        bottom_to_bottom_closed_within_mask=bottom_to_bottom.astype(np.uint8),
        filled_internal_holes=filled_internal_holes.astype(np.uint8),
        x_Mm=data.x_mm.astype(np.float64),
        height_Mm=data.height_mm.astype(np.float64),
        selected_seed_x_height_Mm=selected_x_height,
        selected_seed_xyz_amrvac=selected_xyz,
        selected_seed_grid_indices_iz_ix=np.column_stack(
            (selected_rows, selected_columns)
        ).astype(np.int32),
        selected_seed_log10Q=data.logq[mask].astype(np.float64),
        selected_seed_rboundary=data.rboundary[mask].astype(np.int16),
        selected_seed_is_bottom_to_bottom=bottom_to_bottom[mask].astype(np.uint8),
        selected_seed_is_filled_internal_hole=filled_internal_holes[mask].astype(
            np.uint8
        ),
        right_boundary_grid_indices_iz_ix=boundary_indices,
        right_boundary_x_height_Mm=boundary_x_height,
        q4_seam_path_grid_indices_iz_ix=seam_indices,
        q4_seam_path_x_height_Mm=seam_x_height,
        q4_seam_join_row=np.asarray(
            -1 if seam_state is None else seam_state.join_row, dtype=np.int32
        ),
        q4_seam_join_side=np.asarray(
            "not_applicable" if seam_state is None else seam_state.join_side
        ),
        selected_seed_x_height_Mm_columns=np.asarray(("x_Mm", "height_Mm")),
        selected_seed_xyz_amrvac_columns=np.asarray(
            ("x", "y", "z_absolute")
        ),
        snapshot=np.asarray(data.snapshot, dtype=np.int32),
        time_minutes=np.asarray(
            data.snapshot * TIME_MINUTES_PER_SNAPSHOT, dtype=np.float64
        ),
        source=np.asarray(str(data.source)),
        log10Q_cut=np.asarray(logq_cut, dtype=np.float64),
        fixed_y_amrvac=np.asarray(fixed_y, dtype=np.float64),
        fastqsl2_delta_original_cells=np.asarray(delta, dtype=np.float64),
        length_unit_Mm=np.asarray(LENGTH_UNIT_MM, dtype=np.float64),
        bottom_z_amrvac=np.asarray(BOTTOM_Z_CODE, dtype=np.float64),
        connectivity=np.asarray(4, dtype=np.int8),
        selection_method=np.asarray(selection_mode),
        definition=np.asarray(
            "solid connected row-wise hybrid envelope: left edge from the "
            "rboundary==11 dome; right edge from the same central log10Q=cut "
            "ridge tracked through time; enclosed holes filled"
        ),
    )


def render_frame(
    data: SliceData,
    mask: np.ndarray,
    *,
    closed_component: np.ndarray,
    logq_cut: float,
    view_limits: tuple[float, float, float, float],
    dpi: int,
) -> Image.Image:
    """Render one white-background evolution frame as RGB."""

    figure, axis = plt.subplots(figsize=(7.2, 5.8), constrained_layout=True)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")

    region_layer = np.ma.masked_where(~mask, np.ones(mask.shape, dtype=float))
    axis.pcolormesh(
        data.x_mm,
        data.height_mm,
        region_layer,
        shading="nearest",
        cmap=ListedColormap([MASK_COLOR]),
        vmin=0.0,
        vmax=1.0,
        alpha=1.0,
        rasterized=True,
        zorder=1,
    )
    axis.contour(
        data.x_mm,
        data.height_mm,
        mask.astype(float),
        levels=[0.5],
        colors=[MASK_EDGE_COLOR],
        linewidths=1.15,
        zorder=2,
    )
    # Hide small internal rboundary holes: the requested diagnostic is the
    # outer purple closed/open boundary while the selected region stays solid.
    outer_closed = closed_component
    axis.contour(
        data.x_mm,
        data.height_mm,
        outer_closed.astype(float),
        levels=[0.5],
        colors=[CLOSED_COLOR],
        linewidths=2.0,
        zorder=5,
    )
    finite_logq = data.logq[np.isfinite(data.logq)]
    if finite_logq.size and finite_logq.min() <= logq_cut <= finite_logq.max():
        axis.contour(
            data.x_mm,
            data.height_mm,
            data.logq,
            levels=[logq_cut],
            colors=["white"],
            linewidths=3.8,
            zorder=3.5,
        )
        axis.contour(
            data.x_mm,
            data.height_mm,
            data.logq,
            levels=[logq_cut],
            colors=[Q_COLOR],
            linewidths=2.0,
            zorder=4,
        )

    handles = [
        Patch(
            facecolor=MASK_COLOR,
            edgecolor=MASK_EDGE_COLOR,
            label="Selected cross-section",
        ),
        Line2D(
            [0], [0], color=CLOSED_COLOR, linewidth=2.2,
            label=r"Outer $r_{\rm boundary}=11$ boundary",
        ),
        Line2D(
            [0], [0], color=Q_COLOR, linewidth=2.2,
            label=rf"$\log_{{10}}Q={logq_cut:g}$",
        ),
    ]
    axis.legend(
        handles=handles,
        loc="upper left",
        frameon=True,
        framealpha=0.93,
        facecolor="white",
        edgecolor="#CCCCCC",
        fontsize=11.5,
    )
    axis.set_title(
        f"t = {data.snapshot * TIME_MINUTES_PER_SNAPSHOT:.1f} min",
        fontsize=18,
        pad=8,
    )
    axis.set_xlabel("Solar-X [Mm]", fontsize=14)
    axis.set_ylabel("Height above bottom boundary [Mm]", fontsize=14)
    axis.set_xlim(view_limits[:2])
    axis.set_ylim(view_limits[2:])
    axis.set_aspect("equal")
    axis.tick_params(direction="in", top=True, right=True, labelsize=12)

    buffer = BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi, facecolor="white")
    plt.close(figure)
    buffer.seek(0)
    with Image.open(buffer) as image:
        return image.convert("RGB")


def save_gif(
    path: Path,
    frames: list[Image.Image],
    *,
    fps: float,
    endpoint_pause_seconds: float,
) -> None:
    if not frames:
        raise ValueError("Cannot save an empty GIF")
    base_duration_ms = max(1, int(round(1000.0 / fps)))
    endpoint_duration_ms = max(
        base_duration_ms, int(round(1000.0 * endpoint_pause_seconds))
    )

    # One shared palette avoids temporal color shimmer between GIF frames.
    adaptive = getattr(Image, "Palette", Image).ADAPTIVE
    dither_none = getattr(Image, "Dither", Image).NONE
    first = frames[0].convert("P", palette=adaptive, colors=64)
    paletted = [first]
    paletted.extend(
        frame.quantize(palette=first, dither=dither_none) for frame in frames[1:]
    )
    durations = [base_duration_ms] * len(paletted)
    durations[0] = endpoint_duration_ms
    durations[-1] = endpoint_duration_ms
    first.save(
        path,
        save_all=True,
        append_images=paletted[1:],
        duration=durations,
        loop=0,
        disposal=2,
        optimize=False,
    )


def json_safe(value: Any) -> Any:
    if isinstance(value, Path):
        return str(value.resolve())
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    raise TypeError(f"Cannot encode {type(value).__name__} as JSON")


def write_csv(path: Path, records: list[dict[str, Any]]) -> None:
    if not records:
        raise ValueError("No quality records to write")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--magnetic-dir", type=Path, default=MAGNETIC_DIR)
    parser.add_argument("--slice-dir", type=Path, default=SLICE_DIR)
    parser.add_argument("--axis-feature-dir", type=Path, default=AXIS_FEATURE_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument(
        "--endpoint-reference-dir", type=Path, default=ENDPOINT_REFERENCE_DIR
    )
    parser.add_argument("--start", type=int, default=START_SNAPSHOT)
    parser.add_argument("--stop", type=int, default=STOP_SNAPSHOT)
    parser.add_argument("--step", type=int, default=1)
    parser.add_argument("--fixed-y", type=float, default=FIXED_Y_CODE)
    parser.add_argument("--delta", type=float, default=DELTA_INDEX)
    parser.add_argument("--logq-cut", type=float, default=LOGQ_CUT)
    parser.add_argument("--threads", type=int, default=THREADS)
    parser.add_argument("--fps", type=float, default=8.0)
    parser.add_argument("--endpoint-pause", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=130)
    args = parser.parse_args()

    if args.start != START_SNAPSHOT or args.stop != STOP_SNAPSHOT:
        parser.error(
            "Endpoint regression currently requires --start 10 and --stop 95"
        )
    if args.step <= 0:
        parser.error("--step must be positive")
    if not np.isfinite(args.fixed_y) or args.fixed_y != 0.0:
        parser.error("This workflow is defined for exact --fixed-y 0")
    if not np.isfinite(args.delta) or args.delta <= 0.0:
        parser.error("--delta must be positive and finite")
    if not np.isfinite(args.logq_cut) or args.logq_cut <= 0.0:
        parser.error("--logq-cut must be positive and finite")
    if args.logq_cut != LOGQ_CUT:
        parser.error(
            "Endpoint regression masks are defined for --logq-cut 4; use 4"
        )
    if args.threads <= 0:
        parser.error("--threads must be positive")
    if not np.isfinite(args.fps) or args.fps <= 0.0:
        parser.error("--fps must be positive and finite")
    if not np.isfinite(args.endpoint_pause) or args.endpoint_pause < 0.0:
        parser.error("--endpoint-pause must be non-negative and finite")
    if args.dpi <= 0:
        parser.error("--dpi must be positive")

    snapshots = list(range(args.start, args.stop + 1, args.step))
    if snapshots[-1] != args.stop:
        parser.error("--step must land exactly on snapshot 95")
    args.slice_dir.mkdir(parents=True, exist_ok=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    sources: dict[int, Path] = {}
    cache_status: dict[int, str] = {}
    print(f"Preparing {len(snapshots)} fixed-y slices...", flush=True)
    for snapshot in snapshots:
        source, status = ensure_slice(
            snapshot,
            magnetic_dir=args.magnetic_dir,
            slice_dir=args.slice_dir,
            axis_feature_dir=args.axis_feature_dir,
            fixed_y=args.fixed_y,
            delta=args.delta,
            threads=args.threads,
        )
        sources[snapshot] = source
        cache_status[snapshot] = status

    reference_grid = load_slice(sources[snapshots[0]], snapshots[0])
    for snapshot in snapshots[1:]:
        candidate_grid = load_slice(sources[snapshot], snapshot)
        if not np.array_equal(candidate_grid.x_mm, reference_grid.x_mm) or not np.array_equal(
            candidate_grid.height_mm, reference_grid.height_mm
        ):
            raise ValueError(
                f"Fixed-y grid at t={snapshot} differs from t={snapshots[0]}"
            )

    print("Tracking endpoint-anchored low-Q components...", flush=True)
    forward = forward_track(sources, snapshots, args.logq_cut)
    reverse = reverse_track(sources, snapshots, args.logq_cut)
    merge_window = find_shared_merge_window(forward, reverse, snapshots)
    merge_start = merge_window[0]
    merge_stop = merge_window[-1]
    print(
        f"Central components merge at t={merge_start}..{merge_stop}; "
        "tracking the same Q4 seam through that window...",
        flush=True,
    )

    masks: dict[int, FrameMask] = {}
    for snapshot in snapshots:
        if snapshot in merge_window:
            continue
        choice = forward[snapshot] if snapshot < merge_start else reverse[snapshot]
        selection_mode = (
            "forward_axis_anchored_component"
            if snapshot < merge_start
            else "reverse_axis_anchored_component"
        )
        data = load_slice(sources[snapshot], snapshot)
        mask, bottom_to_bottom, filled_holes, _ = build_hybrid_mask(
            data, choice.mask, args.logq_cut
        )
        native_closed, _, _ = select_closed_component(data, choice.mask)
        statistics = summarize_final_mask(
            data,
            mask,
            choice.mask,
            logq_cut=args.logq_cut,
            selection_mode=selection_mode,
            component_choice=choice,
            seam_state=None,
            left_marker_x_height_mm=LEFT_MARKER_X_HEIGHT_MM,
        )
        masks[snapshot] = FrameMask(
            mask=mask,
            selected_lowq=choice.mask,
            bottom_to_bottom=bottom_to_bottom,
            filled_internal_holes=filled_holes,
            closed_component=binary_fill_holes(native_closed),
            selection_mode=selection_mode,
            seam_state=None,
            statistics=statistics,
        )

    # Evolve the seam state from the first approved endpoint even while the
    # low-Q components are still distinct.  This retains the identity and
    # join-height history of the same Q4 branch before it enters the merge.
    first_snapshot = snapshots[0]
    first_data = load_slice(sources[first_snapshot], first_snapshot)
    seam_state = seam_from_mask(
        first_data,
        masks[first_snapshot].mask,
        logq_cut=args.logq_cut,
    )
    for snapshot in snapshots[1 : snapshots.index(merge_stop) + 1]:
        data = load_slice(sources[snapshot], snapshot)
        seam_state = track_q4_seam(
            data,
            seam_state,
            logq_cut=args.logq_cut,
        )
        if snapshot not in merge_window:
            continue
        mask, seam_state = solid_mask_from_seam(data, seam_state)
        selected_lowq = (
            mask
            & (data.rboundary == 11)
            & np.isfinite(data.logq)
            & (data.logq < args.logq_cut)
        )
        bottom_to_bottom = mask & (data.rboundary == 11)
        filled_holes = mask & (data.rboundary != 11)
        native_closed, _, _ = select_closed_component(data, selected_lowq)
        selection_mode = "dynamic_programming_same_Q4_seam"
        statistics = summarize_final_mask(
            data,
            mask,
            selected_lowq,
            logq_cut=args.logq_cut,
            selection_mode=selection_mode,
            component_choice=None,
            seam_state=seam_state,
            left_marker_x_height_mm=LEFT_MARKER_X_HEIGHT_MM,
        )
        masks[snapshot] = FrameMask(
            mask=mask,
            selected_lowq=selected_lowq,
            bottom_to_bottom=bottom_to_bottom,
            filled_internal_holes=filled_holes,
            closed_component=binary_fill_holes(native_closed),
            selection_mode=selection_mode,
            seam_state=seam_state,
            statistics=statistics,
        )

    endpoint_checks = {
        snapshot: verify_endpoint_mask(
            masks[snapshot].mask,
            snapshot=snapshot,
            reference_dir=args.endpoint_reference_dir,
            logq_cut=args.logq_cut,
        )
        for snapshot in (snapshots[0], snapshots[-1])
    }

    records: list[dict[str, Any]] = []
    frames: list[Image.Image] = []
    mask_paths: list[Path] = []
    frame_paths: list[Path] = []
    frames_dir = args.output_dir / "frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    previous_final_mask: np.ndarray | None = None
    previous_area: float | None = None
    print("Saving masks and rendering GIF frames...", flush=True)
    for snapshot in snapshots:
        data = load_slice(sources[snapshot], snapshot)
        frame_mask = masks[snapshot]
        mask = frame_mask.mask
        selected_lowq = frame_mask.selected_lowq
        bottom_to_bottom = frame_mask.bottom_to_bottom
        filled_holes = frame_mask.filled_internal_holes
        statistics = frame_mask.statistics
        if previous_final_mask is None:
            mask_iou_previous = None
        else:
            mask_intersection = int(
                np.count_nonzero(mask & previous_final_mask)
            )
            mask_union = int(np.count_nonzero(mask | previous_final_mask))
            mask_iou_previous = (
                mask_intersection / mask_union if mask_union else 1.0
            )
        area_change_fraction_previous = (
            None
            if previous_area is None
            else (statistics["area_Mm2"] - previous_area) / previous_area
        )
        topology_change_flag = bool(
            mask_iou_previous is not None
            and (
                mask_iou_previous < 0.5
                or abs(float(area_change_fraction_previous)) > 0.5
            )
        )

        mask_path = args.output_dir / (
            f"y0_hybrid_mask_logq{args.logq_cut:g}_t{snapshot:04d}.npz"
        )
        save_mask_npz(
            mask_path,
            data=data,
            mask=mask,
            selected_lowq_region=selected_lowq,
            selection_mode=frame_mask.selection_mode,
            seam_state=frame_mask.seam_state,
            bottom_to_bottom=bottom_to_bottom,
            filled_internal_holes=filled_holes,
            fixed_y=args.fixed_y,
            delta=args.delta,
            logq_cut=args.logq_cut,
        )
        mask_paths.append(mask_path.resolve())

        axis_point = (
            data.axis_hits_mm[0].tolist()
            if data.axis_hits_mm.shape == (1, 2)
            else [None, None]
        )
        record: dict[str, Any] = {
            "snapshot": snapshot,
            "time_minutes": snapshot * TIME_MINUTES_PER_SNAPSHOT,
            "slice_source": str(data.source),
            "slice_status": cache_status[snapshot],
            "mask_npz": str(mask_path.resolve()),
            "log10Q_cut": args.logq_cut,
            "axis_x_Mm": axis_point[0],
            "axis_height_Mm": axis_point[1],
            "selection_method": frame_mask.selection_mode,
            "mask_iou_previous": mask_iou_previous,
            "area_change_fraction_previous": area_change_fraction_previous,
            "topology_change_flag": topology_change_flag,
            "endpoint_reference_equal": (
                endpoint_checks[snapshot]["pixel_equal"]
                if snapshot in endpoint_checks
                else None
            ),
            **statistics,
        }
        frame = render_frame(
            data,
            mask,
            closed_component=frame_mask.closed_component,
            logq_cut=args.logq_cut,
            view_limits=VIEW_LIMITS_MM,
            dpi=args.dpi,
        )
        frame_path = frames_dir / f"y0_hybrid_mask_logq4_t{snapshot:04d}.png"
        frame.save(frame_path, format="PNG")
        record["frame_png"] = str(frame_path.resolve())
        records.append(record)
        frames.append(frame)
        frame_paths.append(frame_path.resolve())
        previous_final_mask = mask
        previous_area = float(statistics["area_Mm2"])

    stem = (
        f"y0_hybrid_mask_evolution_logq{args.logq_cut:g}_"
        f"t{snapshots[0]:04d}_t{snapshots[-1]:04d}"
    )
    gif_path = args.output_dir / f"{stem}.gif"
    csv_path = args.output_dir / f"{stem}_quality.csv"
    json_path = args.output_dir / f"{stem}_quality.json"
    mask_stack_path = args.output_dir / f"{stem}_mask_stack.npz"
    np.savez_compressed(
        mask_stack_path,
        masks=np.stack([masks[snapshot].mask for snapshot in snapshots]).astype(
            np.uint8
        ),
        snapshots=np.asarray(snapshots, dtype=np.int32),
        time_minutes=np.asarray(snapshots, dtype=np.float64)
        * TIME_MINUTES_PER_SNAPSHOT,
        x_Mm=reference_grid.x_mm.astype(np.float64),
        height_Mm=reference_grid.height_mm.astype(np.float64),
        fixed_y_amrvac=np.asarray(args.fixed_y, dtype=np.float64),
        log10Q_cut=np.asarray(args.logq_cut, dtype=np.float64),
        merge_window_snapshots=np.asarray(merge_window, dtype=np.int32),
        length_unit_Mm=np.asarray(LENGTH_UNIT_MM, dtype=np.float64),
        bottom_z_amrvac=np.asarray(BOTTOM_Z_CODE, dtype=np.float64),
        coordinate_conversion=np.asarray(
            "x_amrvac=x_Mm/10; y_amrvac=0; "
            "z_absolute_amrvac=height_Mm/10+0.1"
        ),
    )
    save_gif(
        gif_path,
        frames,
        fps=args.fps,
        endpoint_pause_seconds=args.endpoint_pause,
    )
    write_csv(csv_path, records)

    topology_change_snapshots = [
        int(record["snapshot"])
        for record in records
        if record["topology_change_flag"]
    ]
    metadata = {
        "implementation": {
            "script": str(Path(__file__).resolve()),
            "seam_tracker_module": str(
                (Path(__file__).resolve().parent / "q4_seam_tracking.py")
            ),
            "reproduction_command": (
                "/opt/anaconda3/bin/python "
                f"{Path(__file__).resolve()}"
            ),
        },
        "method": (
            "real-axis-anchored component tracking outside the merge window; "
            "dynamic-programming tracking of the same bottom-rooted "
            "log10Q=4 seam inside the merge window; solid row-wise fill from "
            "the current rboundary==11 left boundary to that seam"
        ),
        "settings": {
            "snapshots": snapshots,
            "fixed_y_amrvac": args.fixed_y,
            "fastqsl2_delta_original_cells": args.delta,
            "log10Q_cut": args.logq_cut,
            "threads_for_missing_slices": args.threads,
            "view_limits_Mm": list(VIEW_LIMITS_MM),
            "fps": args.fps,
            "endpoint_pause_seconds": args.endpoint_pause,
            "time_minutes_per_snapshot": TIME_MINUTES_PER_SNAPSHOT,
            "length_unit_Mm": LENGTH_UNIT_MM,
            "bottom_z_amrvac": BOTTOM_Z_CODE,
        },
        "tracking": {
            "forward_anchor": "real magnetic-axis intersection at snapshot 10",
            "reverse_anchor": "real magnetic-axis intersection at snapshot 95",
            "merge_window_snapshots": merge_window,
            "before_merge": "forward axis-anchored low-Q component",
            "during_merge": (
                "dynamic-programming path through contemporary row-wise "
                "low-to-high log10Q=4 crossings"
            ),
            "after_merge": "reverse axis-anchored low-Q component",
            "seam_parameters": {
                "maximum_gap_rows": 24,
                "maximum_gap_Mm": 5.625,
                "start_gate_Mm": 4.0,
                "temporal_gate_Mm": 7.0,
                "minimum_join_height_Mm": 5.0,
                "join_boundary_distance_pixels": 2.5,
            },
            "temporal_interpolation_used": False,
            "topology_change_flag_rule": (
                "previous-mask IoU below 0.5 or absolute fractional area "
                "change above 0.5"
            ),
            "topology_change_snapshots": topology_change_snapshots,
        },
        "endpoint_regression": endpoint_checks,
        "saved_seed_coordinates": {
            "selected_seed_x_height_Mm": ["x_Mm", "height_Mm"],
            "selected_seed_xyz_amrvac": ["x", "y", "z_absolute"],
            "conversion": (
                "x=x_Mm/10, y=0, z_absolute=height_Mm/10+0.1"
            ),
            "per_seed_audit_arrays": [
                "selected_seed_grid_indices_iz_ix",
                "selected_seed_log10Q",
                "selected_seed_rboundary",
                "selected_seed_is_bottom_to_bottom",
                "selected_seed_is_filled_internal_hole",
                "right_boundary_grid_indices_iz_ix",
                "right_boundary_x_height_Mm",
                "q4_seam_path_grid_indices_iz_ix",
                "q4_seam_path_x_height_Mm",
            ],
        },
        "outputs": {
            "gif": str(gif_path.resolve()),
            "quality_csv": str(csv_path.resolve()),
            "quality_json": str(json_path.resolve()),
            "mask_stack_npz": str(mask_stack_path.resolve()),
            "mask_npz_files": [str(path) for path in mask_paths],
            "frame_png_files": [str(path) for path in frame_paths],
        },
        "frames": records,
    }
    json_path.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False, default=json_safe) + "\n",
        encoding="utf-8",
    )

    print(gif_path.resolve())
    print(mask_stack_path.resolve())
    print(csv_path.resolve())
    print(json_path.resolve())


if __name__ == "__main__":
    main()

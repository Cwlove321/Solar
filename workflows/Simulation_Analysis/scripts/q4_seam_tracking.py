"""Dynamic-programming tracker for a central ``log10(Q)=4`` seam.

The module tracks the low-to-high Q crossing rooted at the bottom boundary,
then fills the solid region between the current ``rboundary == 11`` left edge
and that seam.  It is used in a merge window where connected-component
tracking alone cannot distinguish the central seam from a large right-hand
branch.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Mapping, Sequence

import numpy as np
from scipy.ndimage import binary_fill_holes, distance_transform_edt, label


FOUR_CONNECTED = np.asarray(
    [[0, 1, 0], [1, 1, 1], [0, 1, 0]], dtype=np.uint8
)


@dataclass(frozen=True)
class SeamState:
    """One tracked central-Q seam and the information needed by the next frame."""

    seam_x_by_row_Mm: dict[int, float]
    join_row: int
    join_side: str
    path_nodes_iz_ix: np.ndarray
    objective: float
    path_mean_cost: float
    audit: dict[str, Any]


def _field(data: Any, name: str) -> np.ndarray:
    """Read an array from either an attribute-based object or a mapping."""

    value = data[name] if isinstance(data, Mapping) else getattr(data, name)
    return np.asarray(value)


def _snapshot(data: Any) -> int | None:
    if isinstance(data, Mapping):
        value = data.get("snapshot")
    else:
        value = getattr(data, "snapshot", None)
    return None if value is None else int(value)


def _grid(data: Any) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    logq = np.asarray(_field(data, "logq"), dtype=float)
    rboundary = np.asarray(_field(data, "rboundary"), dtype=np.int16)
    x_mm = np.asarray(_field(data, "x_mm"), dtype=float)
    height_mm = np.asarray(_field(data, "height_mm"), dtype=float)
    if x_mm.ndim == 2:
        x_mm = x_mm[0]
    if height_mm.ndim == 2:
        height_mm = height_mm[:, 0]
    if logq.shape != rboundary.shape or logq.shape != (height_mm.size, x_mm.size):
        raise ValueError("Inconsistent logq, rboundary, x, and height shapes")
    if x_mm.size < 2 or height_mm.size < 2:
        raise ValueError("The fixed-y grid must have at least two cells per axis")
    return logq, rboundary, x_mm, height_mm


def _interpolate_previous(previous: Mapping[int, float], row: int) -> float:
    rows = np.asarray(sorted(previous), dtype=int)
    if rows.size == 0:
        raise ValueError("The previous seam is empty")
    values = np.asarray([previous[int(index)] for index in rows], dtype=float)
    return float(
        np.interp(row, rows, values, left=values[0], right=values[-1])
    )


def seam_from_mask(
    data: Any,
    mask: np.ndarray,
    *,
    logq_cut: float = 4.0,
) -> SeamState:
    """Initialize a seam from the approved t=10 solid endpoint mask."""

    logq, _, x_mm, height_mm = _grid(data)
    selected = np.asarray(mask, dtype=bool)
    if selected.shape != logq.shape:
        raise ValueError("The initial mask shape differs from the Q slice")
    rows = np.flatnonzero(np.any(selected, axis=1))
    if rows.size == 0:
        raise ValueError("The initial mask is empty")

    right_by_row: dict[int, float] = {}
    q_adjacent_rows: list[int] = []
    path_nodes: list[tuple[int, int]] = []
    for row in rows:
        right = int(np.flatnonzero(selected[row]).max())
        right_by_row[int(row)] = float(x_mm[right])
        if (
            right + 1 < logq.shape[1]
            and np.isfinite(logq[row, right + 1])
            and logq[row, right + 1] >= logq_cut
        ):
            q_adjacent_rows.append(int(row))
            path_nodes.append((int(row), right))
    if not q_adjacent_rows:
        raise ValueError("The initial mask has no right edge adjacent to the Q ridge")
    join_row = max(q_adjacent_rows)
    seam = {
        row: x
        for row, x in right_by_row.items()
        if row <= join_row
    }
    path = np.asarray(
        [(row, int(np.searchsorted(x_mm, x, side="right") - 1))
         for row, x in seam.items()],
        dtype=np.int32,
    )
    return SeamState(
        seam_x_by_row_Mm=seam,
        join_row=join_row,
        join_side="reference",
        path_nodes_iz_ix=path,
        objective=0.0,
        path_mean_cost=0.0,
        audit={
            "snapshot": _snapshot(data),
            "source": "approved mask",
            "join_height_Mm": float(height_mm[join_row]),
            "selected_height_row_count": int(rows.size),
        },
    )


def track_q4_seam(
    data: Any,
    previous: SeamState | Mapping[int, float],
    *,
    previous_join_row: int | None = None,
    logq_cut: float = 4.0,
    x_limits_Mm: tuple[float, float] = (-8.0, 47.0),
    maximum_height_Mm: float = 35.5,
    start_rows: int = 4,
    start_gate_Mm: float = 4.0,
    maximum_gap_rows: int = 24,
    slope_intercept_Mm: float = 1.0,
    maximum_slope: float = 2.5,
    temporal_gate_Mm: float = 7.0,
    temporal_gate_extra_rows: int = 24,
    temporal_scale_Mm: float = 2.0,
    slope_weight: float = 0.06,
    skipped_row_weight: float = 0.35,
    minimum_join_height_Mm: float = 5.0,
    join_boundary_distance_pixels: float = 2.5,
    join_height_scale_Mm: float = 4.0,
) -> SeamState:
    """Track one bottom-rooted Q4 seam using row-wise dynamic programming."""

    logq, rboundary, x_mm, height_mm = _grid(data)
    if isinstance(previous, SeamState):
        previous_seam = previous.seam_x_by_row_Mm
        prior_join_row = previous.join_row
    else:
        previous_seam = {int(row): float(value) for row, value in previous.items()}
        if previous_join_row is None:
            raise ValueError("previous_join_row is required for a seam mapping")
        prior_join_row = int(previous_join_row)

    closed = binary_fill_holes(rboundary == 11)
    distance_to_exterior = distance_transform_edt(closed)
    row_limit = int(np.searchsorted(height_mm, maximum_height_Mm, side="right"))
    row_limit = max(1, min(row_limit, logq.shape[0]))

    candidates: dict[int, list[int]] = {}
    for row in range(row_limit):
        qrow = logq[row]
        columns = np.flatnonzero(
            closed[row, :-1]
            & np.isfinite(qrow[:-1])
            & (qrow[:-1] < logq_cut)
            & np.isfinite(qrow[1:])
            & (qrow[1:] >= logq_cut)
        )
        candidates[row] = [
            int(column)
            for column in columns
            if x_limits_Mm[0] < x_mm[column] < x_limits_Mm[1]
        ]

    previous_bottom_x = _interpolate_previous(previous_seam, 0)
    cost: dict[tuple[int, int], float] = {}
    predecessor: dict[tuple[int, int], tuple[int, int] | None] = {}
    node_count: dict[tuple[int, int], int] = {}
    for row in range(min(start_rows, row_limit)):
        for column in candidates[row]:
            if abs(x_mm[column] - previous_bottom_x) > start_gate_Mm:
                continue
            key = (row, column)
            cost[key] = (
                ((x_mm[column] - previous_bottom_x) / 1.2) ** 2
                + 0.2 * row
            )
            predecessor[key] = None
            node_count[key] = 1

    if not cost:
        raise RuntimeError(
            f"No admissible bottom Q4 root at snapshot {_snapshot(data)}"
        )

    for row in range(1, row_limit):
        for column in candidates[row]:
            x_current = float(x_mm[column])
            temporal_target = _interpolate_previous(previous_seam, row)
            if (
                row <= prior_join_row + temporal_gate_extra_rows
                and abs(x_current - temporal_target) > temporal_gate_Mm
            ):
                continue

            best: tuple[float, tuple[int, int]] | None = None
            for earlier_row in range(max(0, row - maximum_gap_rows), row):
                for earlier_column in candidates[earlier_row]:
                    earlier_key = (earlier_row, earlier_column)
                    if earlier_key not in cost:
                        continue
                    delta_height = float(height_mm[row] - height_mm[earlier_row])
                    delta_x = float(x_current - x_mm[earlier_column])
                    if abs(delta_x) > slope_intercept_Mm + maximum_slope * delta_height:
                        continue
                    temporal_cost = 0.0
                    if row <= prior_join_row + temporal_gate_extra_rows:
                        temporal_cost = (
                            (x_current - temporal_target) / temporal_scale_Mm
                        ) ** 2
                    candidate_cost = (
                        cost[earlier_key]
                        + temporal_cost
                        + slope_weight * (delta_x / delta_height) ** 2
                        + skipped_row_weight * (row - earlier_row - 1)
                    )
                    if best is None or candidate_cost < best[0]:
                        best = (candidate_cost, earlier_key)
            if best is not None:
                key = (row, column)
                cost[key] = best[0]
                predecessor[key] = best[1]
                node_count[key] = node_count[best[1]] + 1

    goals: list[tuple[float, tuple[int, int], float]] = []
    median_dh = float(np.median(np.diff(height_mm)))
    for key, cumulative_cost in cost.items():
        row, column = key
        if height_mm[row] < minimum_join_height_Mm:
            continue
        if distance_to_exterior[row, column] > join_boundary_distance_pixels:
            continue
        mean_cost = cumulative_cost / node_count[key]
        join_height_penalty = (
            ((row - prior_join_row) * median_dh) / join_height_scale_Mm
        ) ** 2
        goals.append((mean_cost + join_height_penalty, key, mean_cost))
    if not goals:
        raise RuntimeError(
            f"No admissible Q4-to-purple-boundary join at snapshot {_snapshot(data)}"
        )

    objective, key, mean_cost = min(goals)
    reverse_path: list[tuple[int, int]] = []
    while key is not None:
        reverse_path.append(key)
        key = predecessor[key]
    path = np.asarray(reverse_path[::-1], dtype=np.int32)
    path_rows = path[:, 0]
    path_x = x_mm[path[:, 1]]
    join_row = int(path_rows[-1])
    all_rows = np.arange(join_row + 1, dtype=int)
    interpolated_x = np.interp(all_rows, path_rows, path_x)
    seam = {
        int(row): float(value)
        for row, value in zip(all_rows, interpolated_x)
    }

    gaps = np.diff(path_rows)
    return SeamState(
        seam_x_by_row_Mm=seam,
        join_row=join_row,
        join_side="unclassified",
        path_nodes_iz_ix=path,
        objective=float(objective),
        path_mean_cost=float(mean_cost),
        audit={
            "snapshot": _snapshot(data),
            "candidate_node_count": int(sum(map(len, candidates.values()))),
            "path_node_count": int(path.shape[0]),
            "maximum_gap_rows": int(gaps.max()) if gaps.size else 0,
            "maximum_gap_Mm": float(gaps.max() * median_dh) if gaps.size else 0.0,
            "join_height_Mm": float(height_mm[join_row]),
            "bottom_x_Mm": float(interpolated_x[0]),
            "join_x_Mm": float(interpolated_x[-1]),
        },
    )


def solid_mask_from_seam(data: Any, seam: SeamState) -> tuple[np.ndarray, SeamState]:
    """Fill a solid mask from the purple left edge to the tracked Q4 seam."""

    _, rboundary, x_mm, height_mm = _grid(data)
    closed = binary_fill_holes(rboundary == 11)
    mask = np.zeros_like(closed, dtype=bool)
    for row, seam_x in seam.seam_x_by_row_Mm.items():
        closed_columns = np.flatnonzero(closed[row])
        if closed_columns.size == 0:
            raise RuntimeError(f"No closed pixels on tracked row {row}")
        right = int(np.searchsorted(x_mm, seam_x, side="right") - 1)
        if right < closed_columns.min():
            raise RuntimeError(f"The tracked seam lies left of the closed region on row {row}")
        mask[row, int(closed_columns.min()) : right + 1] = True

    join_row = seam.join_row
    join_x = seam.seam_x_by_row_Mm[join_row]
    join_column = int(np.searchsorted(x_mm, join_x, side="right") - 1)
    closed_columns = np.flatnonzero(closed[join_row])
    intervals = np.split(
        closed_columns, np.flatnonzero(np.diff(closed_columns) > 1) + 1
    )
    interval = min(intervals, key=lambda values: np.min(np.abs(values - join_column)))
    join_side = (
        "right"
        if int(interval[-1]) - join_column < join_column - int(interval[0])
        else "left"
    )

    if join_side == "right":
        previous_columns = np.flatnonzero(mask[join_row])
        for row in range(join_row + 1, closed.shape[0]):
            closed_columns = np.flatnonzero(closed[row])
            if closed_columns.size == 0:
                break
            intervals = np.split(
                closed_columns, np.flatnonzero(np.diff(closed_columns) > 1) + 1
            )

            def interval_distance(values: np.ndarray) -> int:
                if values[-1] < previous_columns[0]:
                    return int(previous_columns[0] - values[-1])
                if values[0] > previous_columns[-1]:
                    return int(values[0] - previous_columns[-1])
                return 0

            interval = min(intervals, key=interval_distance)
            mask[row, int(interval[0]) : int(interval[-1]) + 1] = True
            previous_columns = interval

    mask = binary_fill_holes(mask)
    component_count = int(label(mask, structure=FOUR_CONNECTED)[1])
    residual_holes = int(np.count_nonzero(binary_fill_holes(mask) & ~mask))
    if component_count != 1 or residual_holes != 0:
        raise RuntimeError(
            f"Solid-mask QA failed: components={component_count}, holes={residual_holes}"
        )

    dx = float(np.median(np.diff(x_mm)))
    dh = float(np.median(np.diff(height_mm)))
    updated_audit = dict(seam.audit)
    updated_audit.update(
        {
            "join_side": join_side,
            "mask_pixel_count": int(np.count_nonzero(mask)),
            "area_Mm2": float(np.count_nonzero(mask) * abs(dx * dh)),
            "mask_component_count": component_count,
            "residual_hole_pixel_count": residual_holes,
        }
    )
    updated = SeamState(
        seam_x_by_row_Mm=seam.seam_x_by_row_Mm,
        join_row=seam.join_row,
        join_side=join_side,
        path_nodes_iz_ix=seam.path_nodes_iz_ix,
        objective=seam.objective,
        path_mean_cost=seam.path_mean_cost,
        audit=updated_audit,
    )
    return mask, updated


def track_sequence(
    datasets: Sequence[Any],
    initial_mask: np.ndarray,
    *,
    logq_cut: float = 4.0,
) -> dict[int, tuple[SeamState, np.ndarray]]:
    """Track a sequence beginning with an approved endpoint mask."""

    if not datasets:
        raise ValueError("No datasets were supplied")
    first = datasets[0]
    state = seam_from_mask(first, initial_mask, logq_cut=logq_cut)
    results: dict[int, tuple[SeamState, np.ndarray]] = {
        int(_snapshot(first) if _snapshot(first) is not None else 0):
        (state, np.asarray(initial_mask, dtype=bool))
    }
    for index, data in enumerate(datasets[1:], start=1):
        state = track_q4_seam(data, state, logq_cut=logq_cut)
        mask, state = solid_mask_from_seam(data, state)
        snapshot = _snapshot(data)
        results[int(snapshot if snapshot is not None else index)] = (state, mask)
    return results

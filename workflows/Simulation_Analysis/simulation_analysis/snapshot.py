"""AMR-aware snapshot inspection and coordinate-aligned slice extraction."""

from __future__ import annotations

import gc
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .geometry import AXES, axis_centers, axis_edges, axis_id, block_centers, block_edge_indices
from .physics import derive_state, gamma_from_dataset, resolve_state_fields


@dataclass
class AxisSlice:
    values: dict[str, np.ndarray]
    level: np.ndarray
    coordinates: tuple[np.ndarray, np.ndarray]
    coordinate_names: tuple[str, str]
    fixed_axis: str
    requested_coordinate: float
    coordinate: float
    angular_unit: str
    target_level: int
    gamma: float
    time: float
    source: str


@dataclass
class PlotPlane:
    x: np.ndarray
    y: np.ndarray
    x_name: str
    y_name: str
    values: dict[str, np.ndarray]
    level: np.ndarray
    projected_bx: np.ndarray
    projected_by: np.ndarray
    fixed_label: str
    time: float
    source: str


def _open_dataset(path: str | Path, ghost_width: int = 0):
    from simesh.amrvac import open_dataset

    return open_dataset(str(Path(path).expanduser().resolve()), ghost_width=ghost_width)


def _native_coordinate(axis: str, value: float, angular_unit: str) -> float:
    if axis == "r" or angular_unit == "radian":
        return float(value)
    if angular_unit != "degree":
        raise ValueError("angular_unit must be degree or radian")
    return float(np.radians(value))


def _display_coordinate(axis: str, value: float, angular_unit: str) -> float:
    if axis == "r" or angular_unit == "radian":
        return float(value)
    return float(np.degrees(value))


def _slice_windows(dataset, fixed_axis: str, level: int, config: dict) -> dict[int, tuple[int, int]]:
    windows = {}
    remaining = [dimension for dimension in range(3) if dimension != axis_id(fixed_axis)]
    radius_unit = float(config["grid"]["radius_unit_Mm"])
    surface = float(config["grid"]["surface_radius_code"])
    max_height = config["slice"].get("max_height_Mm")
    angular_bounds = config["slice"].get("angular_bounds_degree")

    for dimension in remaining:
        centers = axis_centers(dataset, dimension, level, config)
        keep = np.ones(centers.size, dtype=bool)
        if dimension == 0 and max_height is not None:
            keep &= centers <= surface + float(max_height) / radius_unit
        if angular_bounds is not None:
            lower, upper = sorted(float(value) for value in angular_bounds)
            if fixed_axis == "phi" and dimension == 1:
                latitude = 90.0 - np.degrees(centers)
                keep &= (latitude >= lower) & (latitude <= upper)
            elif fixed_axis == "theta" and dimension == 2:
                longitude = np.degrees(centers)
                keep &= (longitude >= lower) & (longitude <= upper)
        indices = np.flatnonzero(keep)
        if indices.size == 0:
            raise ValueError(f"Requested bounds contain no {AXES[dimension]} cells")
        windows[dimension] = (int(indices[0]), int(indices[-1]) + 1)
    return windows


def extract_axis_slice(
    path: str | Path,
    config: dict,
    *,
    axis: str | None = None,
    coordinate: float | None = None,
    angular_unit: str | None = None,
    target_level: int | None = None,
) -> AxisSlice:
    """Extract a piecewise-constant finest-grid slice from spherical AMR data."""
    fixed_axis = (axis or config["slice"]["axis"]).lower()
    fixed_id = axis_id(fixed_axis)
    requested = float(
        config["slice"]["coordinate"] if coordinate is None else coordinate
    )
    angle_unit = angular_unit or config["slice"]["angular_unit"]
    requested_native = _native_coordinate(fixed_axis, requested, angle_unit)

    dataset = _open_dataset(path)
    if int(dataset.ndim) != 3 or "spherical" not in str(dataset.geometry).lower():
        raise ValueError("extract_axis_slice requires a 3D spherical AMRVAC DAT")
    resolved = resolve_state_fields(dataset, config)
    original_indices = [resolved[key][1] for key in resolved]
    dataset.load_data(field_indices=original_indices)
    loaded_column = {logical: position for position, logical in enumerate(resolved)}
    levels = np.asarray(dataset.tree_info[0], dtype=int)
    finest = int(levels.max())
    selected_level = target_level or config["slice"].get("target_level") or finest
    selected_level = int(selected_level)
    if selected_level < finest:
        raise ValueError(
            f"target_level={selected_level} is below the finest leaf level {finest}"
        )

    fixed_edges = axis_edges(dataset, fixed_id, selected_level, config)
    if requested_native < fixed_edges[0] or requested_native > fixed_edges[-1]:
        raise ValueError(
            f"Requested {fixed_axis}={requested_native} lies outside "
            f"[{fixed_edges[0]}, {fixed_edges[-1]}]"
        )
    fixed_index = int(np.searchsorted(fixed_edges, requested_native, side="right") - 1)
    fixed_index = min(max(fixed_index, 0), fixed_edges.size - 2)
    actual_native = float(0.5 * (fixed_edges[fixed_index] + fixed_edges[fixed_index + 1]))

    remaining = [dimension for dimension in range(3) if dimension != fixed_id]
    windows = _slice_windows(dataset, fixed_axis, selected_level, config)
    output_shape = tuple(windows[dimension][1] - windows[dimension][0] for dimension in remaining)
    raw = {
        logical: np.full(output_shape, np.nan, dtype=np.float64)
        for logical in resolved
    }
    level_map = np.zeros(output_shape, dtype=np.int8)

    for block in range(int(dataset.nleafs)):
        level = int(levels[block])
        repeat = 2 ** (selected_level - level)
        start, stop = block_edge_indices(dataset, block, level)
        target_start = start * repeat
        target_stop = stop * repeat
        if not (target_start[fixed_id] <= fixed_index < target_stop[fixed_id]):
            continue
        local_fixed = (fixed_index - target_start[fixed_id]) // repeat

        intersections = []
        for dimension in remaining:
            window_start, window_stop = windows[dimension]
            lower = max(int(target_start[dimension]), window_start)
            upper = min(int(target_stop[dimension]), window_stop)
            if lower >= upper:
                break
            intersections.append((dimension, lower, upper, window_start))
        else:
            output_slices = tuple(
                slice(lower - window_start, upper - window_start)
                for _, lower, upper, window_start in intersections
            )
            for logical, column in loaded_column.items():
                block_values = np.take(
                    dataset.data[block, column], int(local_fixed), axis=fixed_id
                )
                fine = block_values
                for dimension_position in range(2):
                    fine = np.repeat(fine, repeat, axis=dimension_position)
                fine_slices = tuple(
                    slice(lower - int(target_start[dimension]), upper - int(target_start[dimension]))
                    for dimension, lower, upper, _ in intersections
                )
                raw[logical][output_slices] = fine[fine_slices]
            level_map[output_slices] = level

    if np.any(level_map == 0):
        missing = int(np.count_nonzero(level_map == 0))
        raise RuntimeError(f"Slice assembly left {missing} cells unfilled")

    gamma = gamma_from_dataset(dataset, config)
    derived = derive_state(raw, gamma, config)
    coordinates = tuple(
        axis_centers(dataset, dimension, selected_level, config)[
            windows[dimension][0] : windows[dimension][1]
        ]
        for dimension in remaining
    )
    result = AxisSlice(
        values=derived,
        level=level_map,
        coordinates=coordinates,
        coordinate_names=tuple(AXES[dimension] for dimension in remaining),
        fixed_axis=fixed_axis,
        requested_coordinate=requested,
        coordinate=_display_coordinate(fixed_axis, actual_native, angle_unit),
        angular_unit=angle_unit,
        target_level=selected_level,
        gamma=gamma,
        time=float(dataset.metadata.get("time", 0.0)),
        source=str(Path(path).resolve()),
    )
    del dataset
    gc.collect()
    return result


def to_plot_plane(section: AxisSlice, config: dict) -> PlotPlane:
    """Convert native spherical coordinates to a locally physical 2D plane."""
    radius_unit = float(config["grid"]["radius_unit_Mm"])
    surface = float(config["grid"]["surface_radius_code"])
    reference = config["slice"].get("reference_angle_degree")
    values = {name: np.asarray(value) for name, value in section.values.items()}
    level = np.asarray(section.level)

    if section.fixed_axis == "phi":
        r, theta = section.coordinates
        latitude = 90.0 - np.degrees(theta)
        if reference is None:
            reference = 0.5 * (float(latitude[0]) + float(latitude[-1]))
        x = np.radians(latitude[::-1] - float(reference)) * surface * radius_unit
        y = (r - surface) * radius_unit
        values = {name: value[:, ::-1] for name, value in values.items()}
        level = level[:, ::-1]
        bx = -values["b2"]
        by = values["b1"]
        x_name = f"Surface arc distance from latitude {reference:.2f} deg [Mm]"
        fixed_label = f"phi={section.coordinate:.3f} {section.angular_unit}"
    elif section.fixed_axis == "theta":
        r, phi = section.coordinates
        longitude = np.degrees(phi)
        if reference is None:
            reference = 0.5 * (float(longitude[0]) + float(longitude[-1]))
        theta_native = (
            np.radians(section.coordinate)
            if section.angular_unit == "degree"
            else section.coordinate
        )
        x = (
            np.radians(longitude - float(reference))
            * surface
            * radius_unit
            * np.sin(theta_native)
        )
        y = (r - surface) * radius_unit
        bx = values["b3"]
        by = values["b1"]
        x_name = f"Surface arc distance from longitude {reference:.2f} deg [Mm]"
        fixed_label = f"theta={section.coordinate:.3f} {section.angular_unit}"
    else:
        theta, phi = section.coordinates
        latitude = 90.0 - np.degrees(theta)
        longitude = np.degrees(phi)
        x = longitude
        y = latitude[::-1]
        values = {name: value[::-1, :] for name, value in values.items()}
        level = level[::-1, :]
        bx = values["b3"]
        by = -values["b2"]
        x_name = "Longitude [deg]"
        fixed_label = f"r={section.coordinate:.6f}"

    return PlotPlane(
        x=np.asarray(x),
        y=np.asarray(y),
        x_name=x_name,
        y_name="Height above surface [Mm]" if section.fixed_axis != "r" else "Latitude [deg]",
        values=values,
        level=level,
        projected_bx=np.asarray(bx),
        projected_by=np.asarray(by),
        fixed_label=fixed_label,
        time=section.time,
        source=section.source,
    )


def inspect_snapshot(
    path: str | Path,
    config: dict,
    *,
    variable: str = "temperature_K",
    extreme: str = "max",
) -> dict:
    """Summarize a DAT and locate one derived-variable extremum."""
    if extreme not in {"min", "max", "absmax"}:
        raise ValueError("extreme must be min, max, or absmax")
    dataset = _open_dataset(path)
    resolved = resolve_state_fields(dataset, config)
    dataset.load_data(field_indices=[resolved[key][1] for key in resolved])
    gamma = gamma_from_dataset(dataset, config)
    levels = np.asarray(dataset.tree_info[0], dtype=int)
    ranges: dict[str, dict[str, float | int]] = {}
    best_value = None
    best_payload = None

    for block in range(int(dataset.nleafs)):
        raw = {
            logical: dataset.data[block, position]
            for position, logical in enumerate(resolved)
        }
        state = derive_state(raw, gamma, config)
        if variable not in state and variable != "amr_level":
            raise KeyError(f"Unknown variable {variable}; available={sorted(state)}")
        common = (
            "rho_code", "number_density_cm3", "temperature_K", "pressure_code",
            "bmag_G", "vmag_km_s", "sound_speed_km_s", "alfven_speed_km_s",
            "alfven_reduced_km_s", "beta"
        )
        for name in common:
            finite = np.asarray(state[name])[np.isfinite(state[name])]
            if finite.size == 0:
                continue
            entry = ranges.setdefault(name, {"min": float("inf"), "max": float("-inf"), "finite": 0})
            entry["min"] = min(float(entry["min"]), float(finite.min()))
            entry["max"] = max(float(entry["max"]), float(finite.max()))
            entry["finite"] = int(entry["finite"]) + int(finite.size)

        target = (
            np.full_like(raw["rho"], levels[block], dtype=float)
            if variable == "amr_level"
            else np.asarray(state[variable])
        )
        score = np.abs(target) if extreme == "absmax" else target
        flat = int(np.nanargmin(score) if extreme == "min" else np.nanargmax(score))
        local = tuple(int(value) for value in np.unravel_index(flat, target.shape))
        candidate_score = float(score[local])
        better = best_value is None or (
            candidate_score < best_value if extreme == "min" else candidate_score > best_value
        )
        if better:
            r, theta, phi = block_centers(dataset, block, int(levels[block]), config)
            position = (float(r[local[0]]), float(theta[local[1]]), float(phi[local[2]]))
            context = {}
            for name in common:
                context[name] = float(np.asarray(state[name])[local])
            best_value = candidate_score
            best_payload = {
                "value": float(target[local]),
                "block": block,
                "level": int(levels[block]),
                "local_index": list(local),
                "r_code": position[0],
                "height_Mm": (
                    position[0] - float(config["grid"]["surface_radius_code"])
                ) * float(config["grid"]["radius_unit_Mm"]),
                "theta_degree": float(np.degrees(position[1])),
                "latitude_degree": float(90.0 - np.degrees(position[1])),
                "phi_degree": float(np.degrees(position[2])),
                "context": context,
            }

    counts = {str(level): int(np.count_nonzero(levels == level)) for level in sorted(set(levels))}
    report = {
        "source": str(Path(path).expanduser().resolve()),
        "time_code": float(dataset.metadata.get("time", 0.0)),
        "time_s": float(dataset.metadata.get("time", 0.0)) * float(config["units"]["time_s"]),
        "geometry": str(dataset.geometry),
        "domain_nx": [int(value) for value in dataset.domain_nx],
        "block_nx": [int(value) for value in dataset.block_nx],
        "physical_domain": np.asarray(dataset.physical_domain, dtype=float).tolist(),
        "fields": [str(value) for value in dataset.wnames],
        "nleafs": int(dataset.nleafs),
        "leaf_levels": counts,
        "gamma": gamma,
        "ranges": ranges,
        "extremum": {"variable": variable, "kind": extreme, **(best_payload or {})},
    }
    del dataset
    gc.collect()
    return report

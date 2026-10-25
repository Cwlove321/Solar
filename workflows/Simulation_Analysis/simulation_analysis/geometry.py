"""Spherical AMR grid geometry, including the mature-code radial stretch."""

from __future__ import annotations

import numpy as np


AXES = ("r", "theta", "phi")


def axis_id(axis: str) -> int:
    try:
        return AXES.index(axis.lower())
    except ValueError as exc:
        raise ValueError("axis must be r, theta, or phi") from exc


def axis_edges(dataset, dimension: int, level: int, config: dict) -> np.ndarray:
    """Return exact cell edges for one AMR level."""
    factor = 2 ** (int(level) - 1)
    lower, upper = (np.asarray(value, dtype=float) for value in dataset.physical_domain)
    count = int(dataset.domain_nx[dimension]) * factor
    if dimension != 0 or config["grid"]["radial_type"] == "uniform":
        return np.linspace(lower[dimension], upper[dimension], count + 1)

    q_base = float(config["grid"]["qstretch_base"])
    q = q_base ** (1.0 / factor)
    if np.isclose(q, 1.0):
        return np.linspace(lower[dimension], upper[dimension], count + 1)
    dr0 = (upper[dimension] - lower[dimension]) * (q - 1.0) / (q**count - 1.0)
    index = np.arange(count + 1, dtype=float)
    return lower[dimension] + dr0 * (q**index - 1.0) / (q - 1.0)


def axis_centers(dataset, dimension: int, level: int, config: dict) -> np.ndarray:
    edges = axis_edges(dataset, dimension, level, config)
    return 0.5 * (edges[:-1] + edges[1:])


def axis_edges_at_indices(
    dataset, dimension: int, level: int, indices: np.ndarray, config: dict
) -> np.ndarray:
    """Evaluate grid edges at arbitrary logical indices, including ghost indices."""
    factor = 2 ** (int(level) - 1)
    lower, upper = (np.asarray(value, dtype=float) for value in dataset.physical_domain)
    count = int(dataset.domain_nx[dimension]) * factor
    index = np.asarray(indices, dtype=float)
    if dimension != 0 or config["grid"]["radial_type"] == "uniform":
        spacing = (upper[dimension] - lower[dimension]) / count
        return lower[dimension] + index * spacing
    q_base = float(config["grid"]["qstretch_base"])
    q = q_base ** (1.0 / factor)
    if np.isclose(q, 1.0):
        spacing = (upper[dimension] - lower[dimension]) / count
        return lower[dimension] + index * spacing
    dr0 = (upper[dimension] - lower[dimension]) * (q - 1.0) / (q**count - 1.0)
    return lower[dimension] + dr0 * (q**index - 1.0) / (q - 1.0)


def block_edge_indices(dataset, block: int, level: int) -> tuple[np.ndarray, np.ndarray]:
    """Return inclusive-start/exclusive-stop logical indices for a leaf block."""
    logical = np.asarray(dataset.tree_info[1][block], dtype=int)
    block_nx = np.asarray(dataset.block_nx, dtype=int)
    start = (logical - 1) * block_nx
    return start, start + block_nx


def block_centers(dataset, block: int, level: int, config: dict) -> tuple[np.ndarray, ...]:
    start, stop = block_edge_indices(dataset, block, level)
    output = []
    for dimension in range(3):
        edges = axis_edges(dataset, dimension, level, config)
        local = edges[start[dimension] : stop[dimension] + 1]
        output.append(0.5 * (local[:-1] + local[1:]))
    return tuple(output)


def padded_block_centers(
    dataset, block: int, level: int, ghost_width: int, config: dict
) -> tuple[np.ndarray, ...]:
    """Return centers for an AMR block plus exchanged ghost cells."""
    start, stop = block_edge_indices(dataset, block, level)
    output = []
    for dimension in range(3):
        index = np.arange(
            start[dimension] - ghost_width,
            stop[dimension] + ghost_width + 1,
        )
        edges = axis_edges_at_indices(dataset, dimension, level, index, config)
        output.append(0.5 * (edges[:-1] + edges[1:]))
    return tuple(output)


def block_volume(dataset, block: int, level: int, config: dict) -> np.ndarray:
    """Exact spherical leaf-cell volume in code-length cubed."""
    start, stop = block_edge_indices(dataset, block, level)
    r = axis_edges(dataset, 0, level, config)[start[0] : stop[0] + 1]
    theta = axis_edges(dataset, 1, level, config)[start[1] : stop[1] + 1]
    phi = axis_edges(dataset, 2, level, config)[start[2] : stop[2] + 1]
    radial = (r[1:] ** 3 - r[:-1] ** 3) / 3.0
    polar = np.cos(theta[:-1]) - np.cos(theta[1:])
    azimuthal = np.diff(phi)
    return radial[:, None, None] * polar[None, :, None] * azimuthal[None, None, :]


def coordinate_from_index(
    dataset, block: int, level: int, local_index: tuple[int, int, int], config: dict
) -> tuple[float, float, float]:
    centers = block_centers(dataset, block, level, config)
    return tuple(float(centers[d][local_index[d]]) for d in range(3))

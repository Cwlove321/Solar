#!/usr/bin/env python3
"""Locate isolated magnetic nulls in a 3-D rectilinear magnetic field.

The input may be an ``.npz`` file containing three arrays (``Bx``, ``By``,
``Bz`` by default) with shape ``(nx, ny, nz)``, or the QSL/FastQSL binary
format containing an int32 ``(nx, ny, nz)`` header followed by three float64
arrays stored in ``(z, y, x)`` order. Optional one-dimensional ``x``, ``y``,
and ``z`` arrays in an NPZ give the physical grid-point coordinates. If
coordinates are absent, index coordinates or ``--origin``/``--spacing`` are
used.

The search has four stages:

1. vectorized component-sign screening on the eight vertices of every cell;
2. a subdivided solid-angle Poincare degree diagnostic on the cell surface;
3. multi-start, damped Newton iteration of the trilinear field, constrained to
   the candidate cell and using its analytic Jacobian;
4. magnetic-residual validation, global de-duplication, and eigenanalysis of
   the physical-coordinate Jacobian.

The default searches every sign-compatible cell and records the Poincare
degree as a diagnostic. The optional ``--candidate-method degree`` filter is
faster when there are many candidates, but assumes isolated, non-degenerate
nulls with non-zero net topological degree in each cell.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import dataclass
from itertools import product
from pathlib import Path
from typing import Iterable

import numpy as np


# Each face is ordered counter-clockwise when viewed from outside the cell.
# Splitting the six faces gives twelve outward-oriented surface triangles.
_FACES = (
    ((0, 0, 0), (0, 0, 1), (0, 1, 1), (0, 1, 0)),  # -x
    ((1, 0, 0), (1, 1, 0), (1, 1, 1), (1, 0, 1)),  # +x
    ((0, 0, 0), (1, 0, 0), (1, 0, 1), (0, 0, 1)),  # -y
    ((0, 1, 0), (0, 1, 1), (1, 1, 1), (1, 1, 0)),  # +y
    ((0, 0, 0), (0, 1, 0), (1, 1, 0), (1, 0, 0)),  # -z
    ((0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)),  # +z
)
def _plain_array(values: np.ndarray) -> list:
    """Return a JSON-safe nested list, replacing non-finite values with null."""
    array = np.asarray(values)
    if np.iscomplexobj(array):
        raise TypeError("Use separate real and imaginary arrays for JSON output")
    safe = array.astype(object)
    safe[~np.isfinite(array)] = None
    return safe.tolist()


def _finite_or_none(value: float) -> float | None:
    number = float(value)
    return number if np.isfinite(number) else None


def _canonical_vector(vector: np.ndarray) -> np.ndarray:
    """Normalize a real vector and choose a reproducible (but non-physical) sign."""
    real = np.asarray(np.real_if_close(vector), dtype=float)
    norm = float(np.linalg.norm(real))
    if not np.isfinite(norm) or norm == 0.0:
        return np.full(3, np.nan)
    real = real / norm
    pivot = int(np.argmax(np.abs(real)))
    return -real if real[pivot] < 0.0 else real


@dataclass(frozen=True)
class MagneticNull:
    """One validated magnetic null and its local linear topology."""

    cell: tuple[int, int, int]
    index_position: np.ndarray
    position: np.ndarray
    residual: float
    relative_residual: float
    poincare_degree: float
    jacobian: np.ndarray
    eigenvalues: np.ndarray
    eigenvectors: np.ndarray
    null_type: str
    determinant: float
    jacobian_condition_number: float
    divergence: float
    relative_divergence: float
    spine_vector: np.ndarray
    fan_basis: np.ndarray

    def to_dict(self) -> dict:
        return {
            "cell": list(self.cell),
            "index_position": _plain_array(self.index_position),
            "position": {
                name: float(value)
                for name, value in zip(("x", "y", "z"), self.position)
            },
            "residual": self.residual,
            "relative_residual": self.relative_residual,
            "poincare_degree": _finite_or_none(self.poincare_degree),
            "type": self.null_type,
            "determinant": self.determinant,
            "jacobian_condition_number": _finite_or_none(self.jacobian_condition_number),
            "divergence": self.divergence,
            "relative_divergence": self.relative_divergence,
            "jacobian": _plain_array(self.jacobian),
            "eigenvalues": [
                {"real": float(value.real), "imag": float(value.imag)}
                for value in self.eigenvalues
            ],
            # Eigenvectors are columns, matching numpy.linalg.eig.
            "eigenvectors_real": _plain_array(self.eigenvectors.real),
            "eigenvectors_imag": _plain_array(self.eigenvectors.imag),
            "spine_vector": _plain_array(self.spine_vector),
            # Two real vectors spanning the local fan plane.
            "fan_basis": _plain_array(self.fan_basis),
        }


@dataclass(frozen=True)
class NullSearchResult:
    """Nulls plus candidate counts needed to audit a search."""

    nulls: tuple[MagneticNull, ...]
    sign_candidates: int
    searched_candidates: int
    degree_rejected: int
    degree_undefined: int
    degenerate_rejected: int


def _validate_inputs(
    bx: np.ndarray,
    by: np.ndarray,
    bz: np.ndarray,
    x: np.ndarray | None,
    y: np.ndarray | None,
    z: np.ndarray | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, tuple[np.ndarray, np.ndarray, np.ndarray]]:
    fields = tuple(np.asarray(field, dtype=float) for field in (bx, by, bz))
    if any(field.ndim != 3 for field in fields):
        raise ValueError("Bx, By, and Bz must all be three-dimensional arrays")
    if fields[0].shape != fields[1].shape or fields[0].shape != fields[2].shape:
        raise ValueError("Bx, By, and Bz must have identical shapes")
    if any(size < 2 for size in fields[0].shape):
        raise ValueError("Every magnetic-field dimension must contain at least two grid points")
    if not all(np.all(np.isfinite(field)) for field in fields):
        raise ValueError("Magnetic-field arrays contain NaN or infinite values")

    supplied = (x, y, z)
    coordinates = []
    for name, values, size in zip(("x", "y", "z"), supplied, fields[0].shape):
        coordinate = np.arange(size, dtype=float) if values is None else np.asarray(values, dtype=float)
        if coordinate.ndim != 1 or coordinate.size != size:
            raise ValueError(f"{name} must be a one-dimensional array of length {size}")
        if not np.all(np.isfinite(coordinate)):
            raise ValueError(f"{name} contains NaN or infinite values")
        if not np.all(np.diff(coordinate) > 0.0):
            raise ValueError(f"{name} coordinates must be strictly increasing")
        coordinates.append(coordinate)
    return fields[0], fields[1], fields[2], tuple(coordinates)


def _component_sign_mask(field: np.ndarray, atol: float) -> np.ndarray:
    shape = tuple(size - 1 for size in field.shape)
    lower = np.full(shape, np.inf, dtype=float)
    upper = np.full(shape, -np.inf, dtype=float)
    for di, dj, dk in product((0, 1), repeat=3):
        corner = field[
            di : di + shape[0],
            dj : dj + shape[1],
            dk : dk + shape[2],
        ]
        np.minimum(lower, corner, out=lower)
        np.maximum(upper, corner, out=upper)
    return (lower <= atol) & (upper >= -atol)


def component_sign_candidates(
    bx: np.ndarray,
    by: np.ndarray,
    bz: np.ndarray,
    *,
    atol: float = 0.0,
    boundary_cells: int = 0,
) -> np.ndarray:
    """Return cells where every field component brackets zero at its vertices."""
    if atol < 0.0:
        raise ValueError("atol must be non-negative")
    if boundary_cells < 0:
        raise ValueError("boundary_cells must be non-negative")
    mask = (
        _component_sign_mask(np.asarray(bx), atol)
        & _component_sign_mask(np.asarray(by), atol)
        & _component_sign_mask(np.asarray(bz), atol)
    )
    if boundary_cells:
        for axis, size in enumerate(mask.shape):
            if 2 * boundary_cells >= size:
                return np.empty((0, 3), dtype=int)
            selector = [slice(None)] * 3
            selector[axis] = slice(0, boundary_cells)
            mask[tuple(selector)] = False
            selector[axis] = slice(size - boundary_cells, size)
            mask[tuple(selector)] = False
    return np.argwhere(mask)


def _cell_field(bx: np.ndarray, by: np.ndarray, bz: np.ndarray, cell: Iterable[int]) -> np.ndarray:
    i, j, k = (int(value) for value in cell)
    return np.stack(
        (
            bx[i : i + 2, j : j + 2, k : k + 2],
            by[i : i + 2, j : j + 2, k : k + 2],
            bz[i : i + 2, j : j + 2, k : k + 2],
        ),
        axis=-1,
    )


def poincare_degree(
    cell_field: np.ndarray,
    *,
    zero_tol: float = 0.0,
    subdivisions: int = 2,
) -> float:
    """Calculate the magnetic-field topological degree on one cell surface.

    ``cell_field`` has shape ``(2, 2, 2, 3)``. A NaN result means the degree
    is undefined because a sampled surface vertex has nearly zero field.
    """
    values = np.asarray(cell_field, dtype=float)
    if values.shape != (2, 2, 2, 3):
        raise ValueError("cell_field must have shape (2, 2, 2, 3)")
    if zero_tol < 0.0:
        raise ValueError("zero_tol must be non-negative")
    if subdivisions < 1:
        raise ValueError("subdivisions must be at least 1")

    def field_at(local: np.ndarray) -> np.ndarray:
        u, v, w = local
        weights = ((1.0 - u, u), (1.0 - v, v), (1.0 - w, w))
        field = np.zeros(3, dtype=float)
        for i, j, k in product((0, 1), repeat=3):
            field += values[i, j, k] * weights[0][i] * weights[1][j] * weights[2][k]
        return field

    solid_angle = 0.0
    coordinates = np.linspace(0.0, 1.0, subdivisions + 1)
    for face in _FACES:
        p0, p1, p2, p3 = (np.asarray(point, dtype=float) for point in face)
        surface = np.empty((subdivisions + 1, subdivisions + 1, 3), dtype=float)
        for row, u in enumerate(coordinates):
            for column, v in enumerate(coordinates):
                local = (
                    (1.0 - u) * (1.0 - v) * p0
                    + u * (1.0 - v) * p1
                    + u * v * p2
                    + (1.0 - u) * v * p3
                )
                surface[row, column] = field_at(local)
        norms = np.linalg.norm(surface, axis=-1)
        if np.any(norms <= zero_tol):
            return float("nan")
        surface /= norms[..., None]
        for row in range(subdivisions):
            for column in range(subdivisions):
                a = surface[row, column]
                b = surface[row + 1, column]
                c = surface[row + 1, column + 1]
                d = surface[row, column + 1]
                for first, second, third in ((a, b, c), (a, c, d)):
                    numerator = float(np.dot(first, np.cross(second, third)))
                    denominator = float(
                        1.0
                        + np.dot(first, second)
                        + np.dot(second, third)
                        + np.dot(third, first)
                    )
                    solid_angle += 2.0 * np.arctan2(numerator, denominator)
    return float(solid_angle / (4.0 * np.pi))


def trilinear_field_and_jacobian(cell_field: np.ndarray, local: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Evaluate B and dB/d(local coordinate) inside one cell analytically."""
    values = np.asarray(cell_field, dtype=float)
    u, v, w = np.asarray(local, dtype=float)
    wu = np.array((1.0 - u, u))
    wv = np.array((1.0 - v, v))
    ww = np.array((1.0 - w, w))

    field = np.zeros(3, dtype=float)
    for i, j, k in product((0, 1), repeat=3):
        field += values[i, j, k] * wu[i] * wv[j] * ww[k]

    jacobian = np.empty((3, 3), dtype=float)
    jacobian[:, 0] = sum(
        (values[1, j, k] - values[0, j, k]) * wv[j] * ww[k]
        for j, k in product((0, 1), repeat=2)
    )
    jacobian[:, 1] = sum(
        (values[i, 1, k] - values[i, 0, k]) * wu[i] * ww[k]
        for i, k in product((0, 1), repeat=2)
    )
    jacobian[:, 2] = sum(
        (values[i, j, 1] - values[i, j, 0]) * wu[i] * wv[j]
        for i, j in product((0, 1), repeat=2)
    )
    return field, jacobian


def _newton_root(
    cell_field: np.ndarray,
    seed: np.ndarray,
    *,
    residual_limit: float,
    step_tol: float,
    max_iter: int,
    minimum_damping: float = 2.0**-20,
) -> tuple[np.ndarray, float] | None:
    local = np.asarray(seed, dtype=float).copy()
    for _ in range(max_iter):
        field, jacobian = trilinear_field_and_jacobian(cell_field, local)
        residual = float(np.linalg.norm(field))
        if residual <= residual_limit:
            return np.clip(local, 0.0, 1.0), residual
        try:
            step = np.linalg.solve(jacobian, -field)
        except np.linalg.LinAlgError:
            return None
        if not np.all(np.isfinite(step)):
            return None

        damping = 1.0
        accepted = False
        while damping >= minimum_damping:
            trial = local + damping * step
            if np.all(trial >= 0.0) and np.all(trial <= 1.0):
                trial_field, _ = trilinear_field_and_jacobian(cell_field, trial)
                trial_residual = float(np.linalg.norm(trial_field))
                if trial_residual < residual or trial_residual <= residual_limit:
                    local = trial
                    accepted = True
                    break
            damping *= 0.5
        if not accepted:
            return None
        if float(np.linalg.norm(damping * step)) <= step_tol:
            field, _ = trilinear_field_and_jacobian(cell_field, local)
            residual = float(np.linalg.norm(field))
            return (local, residual) if residual <= residual_limit else None

    field, _ = trilinear_field_and_jacobian(cell_field, local)
    residual = float(np.linalg.norm(field))
    return (local, residual) if residual <= residual_limit else None


def _seed_points(seed_grid: int) -> list[np.ndarray]:
    if seed_grid < 1:
        raise ValueError("seed_grid must be at least 1")
    if seed_grid == 1:
        return [np.full(3, 0.5)]
    values = np.linspace(0.15, 0.85, seed_grid)
    seeds = [np.asarray(item, dtype=float) for item in product(values, repeat=3)]
    center = np.full(3, 0.5)
    if not any(np.array_equal(seed, center) for seed in seeds):
        seeds.insert(0, center)
    return seeds


def _orthonormal_fan(first: np.ndarray, second: np.ndarray) -> np.ndarray:
    matrix = np.column_stack((np.asarray(first, dtype=float), np.asarray(second, dtype=float)))
    if not np.all(np.isfinite(matrix)) or np.linalg.matrix_rank(matrix) < 2:
        return np.full((2, 3), np.nan)
    basis, _ = np.linalg.qr(matrix)
    return np.stack((_canonical_vector(basis[:, 0]), _canonical_vector(basis[:, 1])))


def _classify_jacobian(
    jacobian: np.ndarray,
    *,
    eigen_tol: float,
) -> tuple[np.ndarray, np.ndarray, str, np.ndarray, np.ndarray]:
    eigenvalues, eigenvectors = np.linalg.eig(jacobian)
    scale = max(float(np.max(np.abs(eigenvalues))), np.finfo(float).tiny)
    tolerance = eigen_tol * scale
    real_mask = np.abs(eigenvalues.imag) <= tolerance
    real_indices = np.flatnonzero(real_mask)

    spine = np.full(3, np.nan)
    fan = np.full((2, 3), np.nan)
    if real_indices.size == 3:
        real_parts = eigenvalues.real
        positive = np.flatnonzero(real_parts > tolerance)
        negative = np.flatnonzero(real_parts < -tolerance)
        if positive.size == 2 and negative.size == 1:
            spine_index = int(negative[0])
            fan_indices = positive
            kind = "radial_fan_positive_spine_negative"
        elif negative.size == 2 and positive.size == 1:
            spine_index = int(positive[0])
            fan_indices = negative
            kind = "radial_fan_negative_spine_positive"
        else:
            return eigenvalues, eigenvectors, "degenerate_radial", spine, fan
        spine = _canonical_vector(eigenvectors[:, spine_index])
        fan = np.stack([_canonical_vector(eigenvectors[:, index]) for index in fan_indices])
        return eigenvalues, eigenvectors, kind, spine, fan

    if real_indices.size == 1:
        spine_index = int(real_indices[0])
        complex_indices = np.flatnonzero(~real_mask)
        preferred = complex_indices[eigenvalues[complex_indices].imag > 0.0]
        fan_index = int(preferred[0] if preferred.size else complex_indices[0])
        spine_growth = float(eigenvalues[spine_index].real)
        fan_growth = float(eigenvalues[fan_index].real)
        if fan_growth > tolerance and spine_growth < -tolerance:
            kind = "spiral_fan_positive_spine_negative"
        elif fan_growth < -tolerance and spine_growth > tolerance:
            kind = "spiral_fan_negative_spine_positive"
        else:
            kind = "degenerate_spiral"
        spine = _canonical_vector(eigenvectors[:, spine_index])
        complex_vector = eigenvectors[:, fan_index]
        fan = _orthonormal_fan(complex_vector.real, complex_vector.imag)
        return eigenvalues, eigenvectors, kind, spine, fan

    return eigenvalues, eigenvectors, "degenerate_mixed", spine, fan


def _make_null(
    cell: tuple[int, int, int],
    local: np.ndarray,
    residual: float,
    scale: float,
    degree: float,
    cell_field: np.ndarray,
    coordinates: tuple[np.ndarray, np.ndarray, np.ndarray],
    eigen_tol: float,
) -> MagneticNull:
    index_position = np.asarray(cell, dtype=float) + local
    position = np.array(
        [
            axis[index] + local_coordinate * (axis[index + 1] - axis[index])
            for axis, index, local_coordinate in zip(coordinates, cell, local)
        ]
    )
    _, local_jacobian = trilinear_field_and_jacobian(cell_field, local)
    widths = np.array([axis[index + 1] - axis[index] for axis, index in zip(coordinates, cell)])
    jacobian = local_jacobian / widths[np.newaxis, :]
    eigenvalues, eigenvectors, kind, spine, fan = _classify_jacobian(
        jacobian, eigen_tol=eigen_tol
    )
    norm = float(np.linalg.norm(jacobian, ord="fro"))
    divergence = float(np.trace(jacobian))
    return MagneticNull(
        cell=cell,
        index_position=index_position,
        position=position,
        residual=residual,
        relative_residual=residual / max(scale, np.finfo(float).tiny),
        poincare_degree=degree,
        jacobian=jacobian,
        eigenvalues=eigenvalues,
        eigenvectors=eigenvectors,
        null_type=kind,
        determinant=float(np.linalg.det(jacobian)),
        jacobian_condition_number=float(np.linalg.cond(jacobian)),
        divergence=divergence,
        relative_divergence=abs(divergence) / max(norm, np.finfo(float).tiny),
        spine_vector=spine,
        fan_basis=fan,
    )


def _deduplicate(nulls: list[MagneticNull], tolerance: float) -> tuple[MagneticNull, ...]:
    if tolerance < 0.0:
        raise ValueError("dedup_tol must be non-negative")
    unique: list[MagneticNull] = []
    for candidate in sorted(nulls, key=lambda item: item.residual):
        if not any(
            np.linalg.norm(candidate.index_position - accepted.index_position) <= tolerance
            for accepted in unique
        ):
            unique.append(candidate)
    return tuple(sorted(unique, key=lambda item: tuple(item.index_position)))


def find_magnetic_nulls(
    bx: np.ndarray,
    by: np.ndarray,
    bz: np.ndarray,
    *,
    x: np.ndarray | None = None,
    y: np.ndarray | None = None,
    z: np.ndarray | None = None,
    candidate_method: str = "sign",
    sign_atol: float = 0.0,
    degree_threshold: float = 0.5,
    degree_zero_tol: float = 0.0,
    degree_subdivisions: int = 2,
    field_atol: float = 0.0,
    field_rtol: float = 1.0e-10,
    step_tol: float = 1.0e-10,
    max_iter: int = 50,
    seed_grid: int = 2,
    dedup_tol: float = 1.0e-5,
    eigen_tol: float = 1.0e-8,
    jacobian_rcond: float = 1.0e-10,
    keep_degenerate: bool = False,
    boundary_cells: int = 0,
) -> NullSearchResult:
    """Find magnetic nulls and return auditable candidate statistics.

    Positions are returned both in continuous array-index coordinates and in
    the supplied rectilinear physical coordinates. ``dedup_tol`` is measured
    in array-index units, while field tolerances use the input field units.
    """
    bx, by, bz, coordinates = _validate_inputs(bx, by, bz, x, y, z)
    if candidate_method not in {"degree", "sign"}:
        raise ValueError("candidate_method must be 'degree' or 'sign'")
    if not 0.0 <= degree_threshold <= 1.0:
        raise ValueError("degree_threshold must lie between zero and one")
    if min(
        degree_zero_tol,
        field_atol,
        field_rtol,
        step_tol,
        eigen_tol,
        jacobian_rcond,
    ) < 0.0:
        raise ValueError("all tolerances must be non-negative")
    if max_iter < 1:
        raise ValueError("max_iter must be at least 1")
    if degree_subdivisions < 1:
        raise ValueError("degree_subdivisions must be at least 1")

    candidates = component_sign_candidates(
        bx, by, bz, atol=sign_atol, boundary_cells=boundary_cells
    )
    seeds = _seed_points(seed_grid)
    found: list[MagneticNull] = []
    searched = 0
    rejected = 0
    undefined = 0
    degenerate_rejected = 0

    for raw_cell in candidates:
        cell = tuple(int(value) for value in raw_cell)
        values = _cell_field(bx, by, bz, cell)
        degree = poincare_degree(
            values,
            zero_tol=degree_zero_tol,
            subdivisions=degree_subdivisions,
        )
        if not np.isfinite(degree):
            undefined += 1
        if candidate_method == "degree" and np.isfinite(degree) and abs(degree) < degree_threshold:
            rejected += 1
            continue
        searched += 1

        vertex_scale = float(np.max(np.linalg.norm(values, axis=-1)))
        residual_limit = field_atol + field_rtol * vertex_scale
        local_roots: list[tuple[np.ndarray, float]] = []
        for seed in seeds:
            root = _newton_root(
                values,
                seed,
                residual_limit=residual_limit,
                step_tol=step_tol,
                max_iter=max_iter,
            )
            if root is None:
                continue
            local, residual = root
            if any(np.linalg.norm(local - previous[0]) <= dedup_tol for previous in local_roots):
                continue
            local_roots.append((local, residual))

        for local, residual in local_roots:
            null = _make_null(
                cell,
                local,
                residual,
                vertex_scale,
                degree,
                values,
                coordinates,
                eigen_tol,
            )
            reciprocal_condition = 1.0 / null.jacobian_condition_number
            if (
                not keep_degenerate
                and (
                    not np.isfinite(reciprocal_condition)
                    or reciprocal_condition < jacobian_rcond
                )
            ):
                degenerate_rejected += 1
                continue
            found.append(null)

    return NullSearchResult(
        nulls=_deduplicate(found, dedup_tol),
        sign_candidates=int(len(candidates)),
        searched_candidates=searched,
        degree_rejected=rejected,
        degree_undefined=undefined,
        degenerate_rejected=degenerate_rejected,
    )


def _coordinates_from_archive(
    archive: np.lib.npyio.NpzFile,
    shape: tuple[int, int, int],
    keys: tuple[str, str, str],
    origin: tuple[float, float, float],
    spacing: tuple[float, float, float] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    coordinates = []
    for axis, size, key, start, step in zip(
        ("x", "y", "z"), shape, keys, origin, spacing or (1.0, 1.0, 1.0)
    ):
        if key in archive:
            coordinate = np.asarray(archive[key], dtype=float)
        else:
            coordinate = start + step * np.arange(size, dtype=float)
        if coordinate.ndim != 1 or coordinate.size != size:
            raise ValueError(f"coordinate key {key!r} for {axis} must have length {size}")
        coordinates.append(coordinate)
    return tuple(coordinates)


def _load_qsl_binary(path: Path, axis_order: str) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Read the three-component binary format used by the QSL workflows."""
    with path.open("rb") as stream:
        dimensions = np.fromfile(stream, dtype=np.int32, count=3)
        if dimensions.size != 3 or np.any(dimensions < 2):
            raise ValueError(f"Invalid QSL binary dimensions in {path}")
        nx, ny, nz = (int(value) for value in dimensions)
        count = nx * ny * nz
        fields = []
        for component in ("Bx", "By", "Bz"):
            values = np.fromfile(stream, dtype=np.float64, count=count)
            if values.size != count:
                raise ValueError(f"Incomplete {component} array in {path}")
            if axis_order == "zyx":
                field = values.reshape((nz, ny, nx)).transpose(2, 1, 0)
            elif axis_order == "xyz":
                field = values.reshape((nx, ny, nz))
            else:
                raise ValueError("binary axis order must be 'zyx' or 'xyz'")
            fields.append(field)
        if stream.read(1):
            raise ValueError(f"Unexpected trailing data in {path}")
    return fields[0], fields[1], fields[2]


def _save_npz(path: Path, result: NullSearchResult) -> None:
    nulls = result.nulls
    count = len(nulls)
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        path,
        position=np.stack([item.position for item in nulls]) if count else np.empty((0, 3)),
        index_position=(
            np.stack([item.index_position for item in nulls]) if count else np.empty((0, 3))
        ),
        cell=np.asarray([item.cell for item in nulls], dtype=int).reshape(count, 3),
        residual=np.asarray([item.residual for item in nulls]),
        relative_residual=np.asarray([item.relative_residual for item in nulls]),
        poincare_degree=np.asarray([item.poincare_degree for item in nulls]),
        null_type=np.asarray([item.null_type for item in nulls]),
        jacobian=(
            np.stack([item.jacobian for item in nulls]) if count else np.empty((0, 3, 3))
        ),
        eigenvalues=(
            np.stack([item.eigenvalues for item in nulls]) if count else np.empty((0, 3), complex)
        ),
        eigenvectors=(
            np.stack([item.eigenvectors for item in nulls])
            if count
            else np.empty((0, 3, 3), complex)
        ),
        spine_vector=(
            np.stack([item.spine_vector for item in nulls]) if count else np.empty((0, 3))
        ),
        fan_basis=(
            np.stack([item.fan_basis for item in nulls]) if count else np.empty((0, 2, 3))
        ),
    )


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="NPZ or QSL binary rectilinear 3-D field")
    parser.add_argument("--output", type=Path, help="output JSON; default: INPUT_nulls.json")
    parser.add_argument("--save-npz", action="store_true", help="also write a compact NPZ result")
    parser.add_argument(
        "--input-format",
        choices=("auto", "npz", "qsl-binary"),
        default="auto",
        help="auto selects NPZ for .npz and QSL binary for .binary",
    )
    parser.add_argument(
        "--binary-axis-order",
        choices=("zyx", "xyz"),
        default="zyx",
        help="array storage order after the QSL binary header",
    )
    parser.add_argument("--bx-key", default="Bx")
    parser.add_argument("--by-key", default="By")
    parser.add_argument("--bz-key", default="Bz")
    parser.add_argument("--x-key", default="x")
    parser.add_argument("--y-key", default="y")
    parser.add_argument("--z-key", default="z")
    parser.add_argument("--origin", nargs=3, type=float, default=(0.0, 0.0, 0.0), metavar=("X0", "Y0", "Z0"))
    parser.add_argument("--spacing", nargs=3, type=float, metavar=("DX", "DY", "DZ"))
    parser.add_argument("--candidate-method", choices=("degree", "sign"), default="sign")
    parser.add_argument("--sign-atol", type=float, default=0.0)
    parser.add_argument("--degree-threshold", type=float, default=0.5)
    parser.add_argument("--degree-zero-tol", type=float, default=0.0)
    parser.add_argument("--degree-subdivisions", type=int, default=2)
    parser.add_argument("--field-atol", type=float, default=0.0)
    parser.add_argument("--field-rtol", type=float, default=1.0e-10)
    parser.add_argument("--step-tol", type=float, default=1.0e-10)
    parser.add_argument("--max-iter", type=int, default=50)
    parser.add_argument("--seed-grid", type=int, default=2)
    parser.add_argument("--dedup-tol", type=float, default=1.0e-5)
    parser.add_argument("--eigen-tol", type=float, default=1.0e-8)
    parser.add_argument("--jacobian-rcond", type=float, default=1.0e-10)
    parser.add_argument(
        "--keep-degenerate",
        action="store_true",
        help="retain roots whose physical-coordinate Jacobian is singular or ill-conditioned",
    )
    parser.add_argument("--boundary-cells", type=int, default=0)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    input_format = args.input_format
    if input_format == "auto":
        if args.input.suffix.lower() == ".npz":
            input_format = "npz"
        elif args.input.suffix.lower() == ".binary":
            input_format = "qsl-binary"
        else:
            raise SystemExit("Cannot infer input format; use --input-format")

    if input_format == "npz":
        with np.load(args.input, allow_pickle=False) as archive:
            missing = [key for key in (args.bx_key, args.by_key, args.bz_key) if key not in archive]
            if missing:
                raise SystemExit(f"Missing magnetic-field keys: {', '.join(missing)}")
            bx = np.asarray(archive[args.bx_key], dtype=float)
            by = np.asarray(archive[args.by_key], dtype=float)
            bz = np.asarray(archive[args.bz_key], dtype=float)
            coordinates = _coordinates_from_archive(
                archive,
                bx.shape,
                (args.x_key, args.y_key, args.z_key),
                tuple(args.origin),
                tuple(args.spacing) if args.spacing is not None else None,
            )
    else:
        bx, by, bz = _load_qsl_binary(args.input, args.binary_axis_order)
        steps = tuple(args.spacing) if args.spacing is not None else (1.0, 1.0, 1.0)
        coordinates = tuple(
            start + step * np.arange(size, dtype=float)
            for start, step, size in zip(args.origin, steps, bx.shape)
        )

    result = find_magnetic_nulls(
        bx,
        by,
        bz,
        x=coordinates[0],
        y=coordinates[1],
        z=coordinates[2],
        candidate_method=args.candidate_method,
        sign_atol=args.sign_atol,
        degree_threshold=args.degree_threshold,
        degree_zero_tol=args.degree_zero_tol,
        degree_subdivisions=args.degree_subdivisions,
        field_atol=args.field_atol,
        field_rtol=args.field_rtol,
        step_tol=args.step_tol,
        max_iter=args.max_iter,
        seed_grid=args.seed_grid,
        dedup_tol=args.dedup_tol,
        eigen_tol=args.eigen_tol,
        jacobian_rcond=args.jacobian_rcond,
        keep_degenerate=args.keep_degenerate,
        boundary_cells=args.boundary_cells,
    )
    output = args.output or args.input.with_name(f"{args.input.stem}_nulls.json")
    output.parent.mkdir(parents=True, exist_ok=True)
    report = {
        "source": str(args.input.resolve()),
        "input_format": input_format,
        "binary_axis_order": args.binary_axis_order if input_format == "qsl-binary" else None,
        "shape": list(bx.shape),
        "coordinate_ranges": {
            name: [float(axis[0]), float(axis[-1])]
            for name, axis in zip(("x", "y", "z"), coordinates)
        },
        "algorithm": {
            "candidate_method": args.candidate_method,
            "interpolation": "trilinear",
            "root_solver": "cell-constrained damped Newton with analytic Jacobian",
            "sign_atol": args.sign_atol,
            "degree_threshold": args.degree_threshold,
            "degree_zero_tol": args.degree_zero_tol,
            "degree_subdivisions": args.degree_subdivisions,
            "field_atol": args.field_atol,
            "field_rtol": args.field_rtol,
            "step_tol": args.step_tol,
            "max_iter": args.max_iter,
            "seed_grid": args.seed_grid,
            "dedup_tol_index_units": args.dedup_tol,
            "eigen_tol": args.eigen_tol,
            "jacobian_rcond": args.jacobian_rcond,
            "keep_degenerate": args.keep_degenerate,
            "boundary_cells": args.boundary_cells,
        },
        "candidate_counts": {
            "component_sign": result.sign_candidates,
            "searched": result.searched_candidates,
            "degree_rejected": result.degree_rejected,
            "degree_undefined": result.degree_undefined,
            "degenerate_rejected": result.degenerate_rejected,
        },
        "null_count": len(result.nulls),
        "nulls": [item.to_dict() for item in result.nulls],
    }
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    if args.save_npz:
        _save_npz(output.with_suffix(".npz"), result)
    print(json.dumps({"output": str(output.resolve()), "null_count": len(result.nulls)}, indent=2))


if __name__ == "__main__":
    main()

"""Helicity and energy flux helpers."""

from __future__ import annotations

import numpy as np


def h_flux(
    vector_potential: list[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray],
    magnetic_field: list[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray],
    velocity: list[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray],
    dx: float = 1.0,
) -> float:
    """Calculate magnetic helicity flux through the six faces of a 3D box.

    The component order follows the original research scripts: `[component3,
    component2, component1]` on arrays indexed as `(z, y, x)`.
    """

    flux = 0.0
    for i in range(3):
        for n in (0, -1):
            ap1 = np.take(vector_potential[(i + 2) % 3], indices=n, axis=i)
            ap2 = np.take(vector_potential[(i + 1) % 3], indices=n, axis=i)
            b1 = np.take(magnetic_field[(i + 2) % 3], indices=n, axis=i)
            b2 = np.take(magnetic_field[(i + 1) % 3], indices=n, axis=i)
            v1 = np.take(velocity[(i + 2) % 3], indices=n, axis=i)
            v2 = np.take(velocity[(i + 1) % 3], indices=n, axis=i)
            b3 = np.take(magnetic_field[i], indices=n, axis=i)
            v3 = np.take(velocity[i], indices=n, axis=i)
            flux += 2.0 * np.sum((ap1 * b1 + ap2 * b2) * v3 - (ap1 * v1 + ap2 * v2) * b3) * (-1) ** n
    return float(flux * dx**2)


def e_flux(
    magnetic_field: list[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray],
    velocity: list[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray],
    dx: float = 1.0,
) -> float:
    """Calculate Poynting flux through the six faces of a 3D box."""

    flux = 0.0
    for i in range(3):
        for n in (0, -1):
            b1 = np.take(magnetic_field[(i + 2) % 3], indices=n, axis=i)
            b2 = np.take(magnetic_field[(i + 1) % 3], indices=n, axis=i)
            v1 = np.take(velocity[(i + 2) % 3], indices=n, axis=i)
            v2 = np.take(velocity[(i + 1) % 3], indices=n, axis=i)
            b3 = np.take(magnetic_field[i], indices=n, axis=i)
            v3 = np.take(velocity[i], indices=n, axis=i)
            flux += np.sum((b1**2 + b2**2) * v3 - (b1 * v1 + b2 * v2) * b3) * (-1) ** n / (4.0 * np.pi)
    return float(flux * dx**2)


def h_transpose(
    potential_field: list[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray],
    magnetic_field: list[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray],
    velocity: list[np.ndarray] | tuple[np.ndarray, np.ndarray, np.ndarray],
    dx: float = 1.0,
) -> float:
    """Calculate helicity transfer from ``Hj`` to ``Hpj``."""

    bp3, bp2, bp1 = potential_field
    b3, b2, b1 = magnetic_field
    v3, v2, v1 = velocity

    c1 = v2 * b3 - v3 * b2
    c2 = b3 * b1 - v1 * b3
    c3 = v1 * b2 - v2 * b1

    return float(-2.0 * volume_integral(c1 * bp1 + c2 * bp2 + c3 * bp3, dx=dx))


def volume_integral(array: np.ndarray, dx: float = 1.0) -> float:
    """Integrate a 3D cell-centered array with boundary half-cell weights."""

    array = np.asarray(array)
    if array.ndim != 3:
        raise ValueError("volume_integral expects a 3D array")
    if min(array.shape) < 2:
        raise ValueError("array must have at least two cells along every axis")

    corners = [
        array[0, 0, 0],
        array[0, 0, -1],
        array[0, -1, 0],
        array[0, -1, -1],
        array[-1, 0, 0],
        array[-1, 0, -1],
        array[-1, -1, 0],
        array[-1, -1, -1],
    ]
    edges = [
        array[0, 0, 1:-1],
        array[0, -1, 1:-1],
        array[-1, 0, 1:-1],
        array[-1, -1, 1:-1],
        array[0, 1:-1, 0],
        array[0, 1:-1, -1],
        array[-1, 1:-1, 0],
        array[-1, 1:-1, -1],
        array[1:-1, 0, 0],
        array[1:-1, 0, -1],
        array[1:-1, -1, 0],
        array[1:-1, -1, -1],
    ]
    faces = [
        array[0, 1:-1, 1:-1],
        array[-1, 1:-1, 1:-1],
        array[1:-1, 0, 1:-1],
        array[1:-1, -1, 1:-1],
        array[1:-1, 1:-1, 0],
        array[1:-1, 1:-1, -1],
    ]

    total = 0.125 * np.sum(corners)
    total += 0.25 * sum(np.sum(edge) for edge in edges)
    total += 0.5 * sum(np.sum(face) for face in faces)
    total += np.sum(array[1:-1, 1:-1, 1:-1])
    return float(total * dx**3)


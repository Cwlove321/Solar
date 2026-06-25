"""Binary I/O helpers for simulation magnetic-field cubes."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class MagneticCube:
    """A 3D magnetic-field cube stored as arrays indexed `(z, y, x)`."""

    b1: np.ndarray
    b2: np.ndarray
    b3: np.ndarray
    dx: float | None = None

    def __post_init__(self) -> None:
        if self.b1.shape != self.b2.shape or self.b1.shape != self.b3.shape:
            raise ValueError("b1, b2, and b3 must have the same shape")
        if self.b1.ndim != 3:
            raise ValueError("magnetic-field components must be 3D arrays")


def write_magnetic_cube(path: str | Path, cube: MagneticCube, include_dx: bool = True) -> None:
    """Write a magnetic-field cube to the project binary format."""

    path = Path(path)
    nz, ny, nx = cube.b1.shape
    with path.open("wb") as stream:
        if include_dx:
            if cube.dx is None:
                raise ValueError("cube.dx is required when include_dx=True")
            np.array(cube.dx, dtype=np.double).tofile(stream)
        np.array(nx, dtype=np.int32).tofile(stream)
        np.array(ny, dtype=np.int32).tofile(stream)
        np.array(nz, dtype=np.int32).tofile(stream)
        np.asarray(cube.b1, dtype=np.double).tofile(stream)
        np.asarray(cube.b2, dtype=np.double).tofile(stream)
        np.asarray(cube.b3, dtype=np.double).tofile(stream)


def read_magnetic_cube(path: str | Path, has_dx: bool = True) -> MagneticCube:
    """Read a magnetic-field cube written by :func:`write_magnetic_cube`."""

    path = Path(path)
    with path.open("rb") as stream:
        dx = float(np.fromfile(stream, dtype=np.double, count=1)[0]) if has_dx else None
        nx, ny, nz = np.fromfile(stream, dtype=np.int32, count=3)
        count = int(nx * ny * nz)
        b1 = np.fromfile(stream, dtype=np.double, count=count).reshape((nz, ny, nx))
        b2 = np.fromfile(stream, dtype=np.double, count=count).reshape((nz, ny, nx))
        b3 = np.fromfile(stream, dtype=np.double, count=count).reshape((nz, ny, nx))
    return MagneticCube(b1=b1, b2=b2, b3=b3, dx=dx)


def init_boundary(path: str | Path, b1: np.ndarray, b2: np.ndarray, b3: np.ndarray) -> float:
    """Write six boundary faces for arrays indexed as `(z, y, x)`.

    Returns the total signed boundary flux, matching the diagnostic value used
    in the original scripts.
    """

    b1 = np.asarray(b1)
    b2 = np.asarray(b2)
    b3 = np.asarray(b3)
    if b1.shape != b2.shape or b1.shape != b3.shape:
        raise ValueError("b1, b2, and b3 must have the same shape")
    if b1.ndim != 3:
        raise ValueError("boundary export expects 3D arrays")

    faces = [
        -b1[:, :, 0],
        b1[:, :, -1],
        -b2[:, 0, :],
        b2[:, -1, :],
        -b3[0, :, :],
        b3[-1, :, :],
    ]

    total_flux = 0.0
    with Path(path).open("wb") as stream:
        for face in faces:
            np.asarray(face, dtype=np.double).tofile(stream)
            total_flux += float(np.sum(face))
    return total_flux


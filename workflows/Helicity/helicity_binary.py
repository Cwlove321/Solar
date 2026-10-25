"""Binary contract shared by the AMRVAC preprocessor and Helicity solver.

Python fields use ``(z, y, x)`` C-order arrays.  Their byte stream is identical
to a Fortran ``(x, y, z)`` array because ``x`` is the fastest-varying index in
both layouts.
"""

from __future__ import annotations

from collections.abc import Iterator, Sequence
from contextlib import contextmanager
import os
from pathlib import Path
import tempfile

import numpy as np
from numpy.typing import NDArray


FLOAT_DTYPE = np.dtype("<f8")
INT_DTYPE = np.dtype("<i4")
FloatArray = NDArray[np.float64]
BoundaryFaces = tuple[FloatArray, FloatArray, FloatArray, FloatArray, FloatArray, FloatArray]


def _validated_components(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
) -> tuple[FloatArray, FloatArray, FloatArray]:
    """Return same-shaped, three-dimensional, C-contiguous float64 fields."""
    components = tuple(
        np.asarray(component, dtype=FLOAT_DTYPE, order="C")
        for component in (b1, b2, b3)
    )
    shape = components[0].shape
    if len(shape) != 3:
        raise ValueError(f"磁场数组必须是三维 (z,y,x)，实际形状为 {shape}")
    if any(component.shape != shape for component in components[1:]):
        raise ValueError("b1、b2、b3 的形状必须一致")
    if any(size < 1 for size in shape):
        raise ValueError(f"磁场数组各维长度必须大于零，实际形状为 {shape}")
    return components


@contextmanager
def _atomic_binary_writer(filename: str | Path) -> Iterator[object]:
    """Write beside the destination and replace it only after a clean close."""
    destination = Path(filename)
    destination.parent.mkdir(parents=True, exist_ok=True)
    temporary_name: str | None = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="wb",
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as stream:
            temporary_name = stream.name
            yield stream
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary_name, destination)
        temporary_name = None
    finally:
        if temporary_name is not None:
            Path(temporary_name).unlink(missing_ok=True)


def boundary_faces(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
) -> BoundaryFaces:
    """Return outward-normal fields in MKL order: ax, bx, ay, by, az, bz."""
    return _boundary_faces_from_components(_validated_components(b1, b2, b3))


def _boundary_faces_from_components(
    components: tuple[FloatArray, FloatArray, FloatArray],
) -> BoundaryFaces:
    bx_field, by_field, bz_field = components
    return (
        -bx_field[:, :, 0],
        bx_field[:, :, -1],
        -by_field[:, 0, :],
        by_field[:, -1, :],
        -bz_field[0, :, :],
        bz_field[-1, :, :],
    )


def _write_boundary_faces(filename: str | Path, faces: BoundaryFaces) -> float:
    with _atomic_binary_writer(filename) as stream:
        for face in faces:
            np.asarray(face, dtype=FLOAT_DTYPE, order="C").tofile(stream)
    return float(sum(np.sum(face, dtype=np.float64) for face in faces))


def _write_magnetic_components(
    filename: str | Path,
    spacing: float,
    components: tuple[FloatArray, FloatArray, FloatArray],
) -> None:
    nz, ny, nx = components[0].shape
    with _atomic_binary_writer(filename) as stream:
        np.asarray(spacing, dtype=FLOAT_DTYPE).tofile(stream)
        np.asarray((nx, ny, nz), dtype=INT_DTYPE).tofile(stream)
        for component in components:
            component.tofile(stream)


def write_boundary_data(
    filename: str | Path,
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
) -> float:
    """Write six boundary faces and return their unweighted discrete sum."""
    faces = _boundary_faces_from_components(_validated_components(b1, b2, b3))
    return _write_boundary_faces(filename, faces)


def write_magnetic_data(
    filename: str | Path,
    spacing: float,
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
) -> None:
    """Write ``spacing, nx, ny, nz, Bx, By, Bz`` as a Fortran stream."""
    components = _validated_components(b1, b2, b3)
    _write_magnetic_components(filename, spacing, components)


def write_helicity_inputs(
    magnetic_filename: str | Path,
    boundary_filename: str | Path,
    spacing: float,
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
) -> float:
    """Validate once and write the paired magnetic and boundary input files."""
    components = _validated_components(b1, b2, b3)
    faces = _boundary_faces_from_components(components)
    _write_magnetic_components(magnetic_filename, spacing, components)
    return _write_boundary_faces(boundary_filename, faces)


def _read_exact(stream: object, dtype: np.dtype, count: int, label: str) -> np.ndarray:
    values = np.fromfile(stream, dtype=dtype, count=count)
    if values.size != count:
        raise ValueError(f"{label}不完整：预期 {count} 个值，实际 {values.size} 个")
    return values


def _reject_trailing_bytes(stream: object, filename: str | Path) -> None:
    if stream.read(1):
        raise ValueError(f"文件包含未定义的尾部数据：{filename}")


def read_magnetic_data(filename: str | Path) -> tuple[float, FloatArray, FloatArray, FloatArray]:
    """Read and strictly validate a magnetic-input stream."""
    with Path(filename).open("rb") as stream:
        spacing = _read_exact(stream, FLOAT_DTYPE, 1, "网格间距")
        dimensions = _read_exact(stream, INT_DTYPE, 3, "网格尺寸")
        nx, ny, nz = (int(value) for value in dimensions)
        if min(nx, ny, nz) < 1:
            raise ValueError(f"非法网格尺寸：{(nx, ny, nz)}")
        count = nx * ny * nz
        components = [
            _read_exact(stream, FLOAT_DTYPE, count, f"磁场分量 {index}")
            for index in range(1, 4)
        ]
        _reject_trailing_bytes(stream, filename)
    fields = tuple(component.reshape(nz, ny, nx) for component in components)
    return float(spacing[0]), fields[0], fields[1], fields[2]


def read_boundary_data(filename: str | Path, shape: Sequence[int]) -> BoundaryFaces:
    """Read six boundary faces for a magnetic field with ``(z,y,x)`` shape."""
    if len(shape) != 3:
        raise ValueError(f"shape 必须是 (z,y,x)，实际为 {tuple(shape)}")
    nz, ny, nx = (int(value) for value in shape)
    if min(nx, ny, nz) < 1:
        raise ValueError(f"非法网格尺寸：{(nx, ny, nz)}")
    face_shapes = ((nz, ny), (nz, ny), (nz, nx), (nz, nx), (ny, nx), (ny, nx))
    faces: list[FloatArray] = []
    with Path(filename).open("rb") as stream:
        for index, face_shape in enumerate(face_shapes, start=1):
            count = int(np.prod(face_shape))
            face = _read_exact(stream, FLOAT_DTYPE, count, f"边界面 {index}")
            faces.append(face.reshape(face_shape))
        _reject_trailing_bytes(stream, filename)
    return tuple(faces)  # type: ignore[return-value]

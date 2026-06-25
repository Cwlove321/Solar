"""Vector calculus helpers for gridded magnetic fields."""

from __future__ import annotations

from typing import Literal

import numpy as np

ArrayLike = np.ndarray
AxisOrder = Literal["zyx", "xyz"]


def curl(
    b1: ArrayLike,
    b2: ArrayLike,
    b3: ArrayLike,
    dx: float = 1.0,
    dy: float = 1.0,
    dz: float = 1.0,
    axis_order: AxisOrder = "zyx",
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return curl of a vector field.

    Parameters
    ----------
    b1, b2, b3:
        Vector components on the same 3D grid.
    dx, dy, dz:
        Grid spacing along x, y, and z.
    axis_order:
        `"zyx"` means arrays are indexed as `(z, y, x)`, matching most of the
        existing research scripts. `"xyz"` means arrays are indexed as `(x, y, z)`.
    """

    b1 = np.asarray(b1)
    b2 = np.asarray(b2)
    b3 = np.asarray(b3)
    if b1.shape != b2.shape or b1.shape != b3.shape:
        raise ValueError("b1, b2, and b3 must have the same shape")
    if b1.ndim != 3:
        raise ValueError("curl expects 3D arrays")

    if axis_order == "zyx":
        curl_x = np.gradient(b3, axis=1) / dy - np.gradient(b2, axis=0) / dz
        curl_y = np.gradient(b1, axis=0) / dz - np.gradient(b3, axis=2) / dx
        curl_z = np.gradient(b2, axis=2) / dx - np.gradient(b1, axis=1) / dy
    elif axis_order == "xyz":
        curl_x = np.gradient(b3, axis=1) / dy - np.gradient(b2, axis=2) / dz
        curl_y = np.gradient(b1, axis=2) / dz - np.gradient(b3, axis=0) / dx
        curl_z = np.gradient(b2, axis=0) / dx - np.gradient(b1, axis=1) / dy
    else:
        raise ValueError("axis_order must be 'zyx' or 'xyz'")

    return curl_x, curl_y, curl_z


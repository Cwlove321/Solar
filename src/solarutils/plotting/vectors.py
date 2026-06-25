"""Vector-field plotting helpers."""

from __future__ import annotations

import numpy as np
import matplotlib.pyplot as plt


def plot_signed_vectors(
    bz: np.ndarray,
    vx: np.ndarray,
    vy: np.ndarray,
    interval: int,
    scale: float,
    ax=None,
    bzmax: float = 1000.0,
    threshold: float = 0.1,
):
    """Plot transverse vectors colored by the sign of the background field."""

    if ax is None:
        _, ax = plt.subplots(dpi=200)

    bz = np.asarray(bz)
    vx = np.asarray(vx)
    vy = np.asarray(vy)
    if bz.shape != vx.shape or bz.shape != vy.shape:
        raise ValueError("bz, vx, and vy must have the same 2D shape")

    ny, nx = bz.shape
    x, y = np.meshgrid(np.arange(nx), np.arange(ny))
    ax.imshow(bz, cmap="gray", origin="lower", vmax=bzmax, vmin=-bzmax)

    mask = np.sqrt(vx**2 + vy**2) * np.sign(bz)

    vx_pos = vx.copy()
    vy_pos = vy.copy()
    vx_pos[mask < threshold] = np.nan
    vy_pos[mask < threshold] = np.nan
    quiver = ax.quiver(
        x[::interval, ::interval],
        y[::interval, ::interval],
        vx_pos[::interval, ::interval],
        vy_pos[::interval, ::interval],
        scale=scale,
        color="skyblue",
        width=0.002,
        angles="xy",
    )

    vx_neg = vx.copy()
    vy_neg = vy.copy()
    vx_neg[mask > -threshold] = np.nan
    vy_neg[mask > -threshold] = np.nan
    ax.quiver(
        x[::interval, ::interval],
        y[::interval, ::interval],
        vx_neg[::interval, ::interval],
        vy_neg[::interval, ::interval],
        scale=scale,
        color="tab:pink",
        width=0.002,
        angles="xy",
    )

    return quiver


"""Consistent plot styles for simulation diagnostics."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, LogNorm, Normalize, TwoSlopeNorm
from scipy.interpolate import RegularGridInterpolator

from .snapshot import PlotPlane


FIELD_STYLES = {
    "temperature_K": ("Temperature [K]", "inferno", "log", "temperature_K"),
    "number_density_cm3": (r"Number density [cm$^{-3}$]", "viridis", "log", "number_density_cm3"),
    "rho_code": ("Density [code units]", "viridis", "log", None),
    "pressure_code": ("Pressure [code units]", "magma", "log", None),
    "entropy_proxy": (r"$p/rho^gamma$ [code units]", "cividis", "log", None),
    "br_G": (r"$B_r$ [G]", "coolwarm", "symmetric", "magnetic_field_G"),
    "btheta_G": (r"$B_theta$ [G]", "coolwarm", "symmetric", "magnetic_field_G"),
    "bphi_G": (r"$B_phi$ [G]", "coolwarm", "symmetric", "magnetic_field_G"),
    "bmag_G": (r"$|B|$ [G]", "magma", "log", None),
    "vr_km_s": (r"$v_r$ [km s$^{-1}$]", "coolwarm", "symmetric", "velocity_km_s"),
    "vtheta_km_s": (r"$v_theta$ [km s$^{-1}$]", "coolwarm", "symmetric", "velocity_km_s"),
    "vphi_km_s": (r"$v_phi$ [km s$^{-1}$]", "coolwarm", "symmetric", "velocity_km_s"),
    "vmag_km_s": (r"$|v|$ [km s$^{-1}$]", "plasma", "linear", None),
    "sound_speed_km_s": (r"$c_s$ [km s$^{-1}$]", "plasma", "log", None),
    "alfven_speed_km_s": (r"$v_A$ [km s$^{-1}$]", "plasma", "log", None),
    "alfven_reduced_km_s": (r"Reduced $v_A$ [km s$^{-1}$]", "plasma", "log", None),
    "fast_speed_upper_km_s": (r"$(c_s^2+v_A^2)^{1/2}$ [km s$^{-1}$]", "plasma", "log", None),
    "beta": (r"Plasma $beta$", "cividis", "log", None),
    "amr_level": ("AMR level", "viridis", "level", None),
}


def _finite_limits(values: np.ndarray, mode: str) -> tuple[float, float]:
    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return 0.0, 1.0
    if mode == "log":
        positive = finite[finite > 0.0]
        if positive.size == 0:
            return 1.0e-30, 1.0
        return float(np.percentile(positive, 1.0)), float(np.percentile(positive, 99.0))
    if mode == "symmetric":
        limit = max(float(np.percentile(np.abs(finite), 99.0)), 1.0e-30)
        return -limit, limit
    lower, upper = np.percentile(finite, [1.0, 99.0])
    if lower == upper:
        upper = lower + max(abs(lower), 1.0) * 1.0e-6
    return float(lower), float(upper)


def field_style(field: str, values: np.ndarray, config: dict):
    label, cmap, mode, config_range = FIELD_STYLES.get(
        field, (field, "viridis", "linear", None)
    )
    limits = config["plot"].get(config_range) if config_range else None
    vmin, vmax = tuple(float(value) for value in limits) if limits else _finite_limits(values, mode)
    if mode == "log":
        norm = LogNorm(max(vmin, 1.0e-30), max(vmax, max(vmin, 1.0e-30) * 1.0001))
    elif mode == "symmetric":
        limit = max(abs(vmin), abs(vmax), 1.0e-30)
        norm = TwoSlopeNorm(vcenter=0.0, vmin=-limit, vmax=limit)
    elif mode == "level":
        low = int(np.nanmin(values))
        high = int(np.nanmax(values))
        norm = BoundaryNorm(np.arange(low - 0.5, high + 1.5), high - low + 1)
    else:
        norm = Normalize(vmin, vmax)
    return label, cmap, norm


def add_fieldlines(axis, plane: PlotPlane, *, color: str = "white") -> None:
    """Overlay the magnetic field projected into the displayed plane."""
    nx = min(180, plane.x.size)
    ny = min(140, plane.y.size)
    x = np.linspace(float(plane.x[0]), float(plane.x[-1]), nx)
    y = np.linspace(float(plane.y[0]), float(plane.y[-1]), ny)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    points = np.column_stack((yy.ravel(), xx.ravel()))
    bx = RegularGridInterpolator(
        (plane.y, plane.x), plane.projected_bx, bounds_error=False, fill_value=np.nan
    )(points).reshape(ny, nx)
    by = RegularGridInterpolator(
        (plane.y, plane.x), plane.projected_by, bounds_error=False, fill_value=np.nan
    )(points).reshape(ny, nx)
    bx = np.ma.masked_invalid(bx)
    by = np.ma.masked_invalid(by)
    axis.streamplot(
        x,
        y,
        bx,
        by,
        color=color,
        density=1.0,
        linewidth=0.45,
        arrowsize=0.55,
    )


def draw_scalar(axis, plane: PlotPlane, values: np.ndarray, *, cmap, norm):
    """Draw cell-centered values without pretending a stretched axis is uniform."""
    return axis.pcolormesh(
        plane.x,
        plane.y,
        values,
        shading="nearest",
        cmap=cmap,
        norm=norm,
        rasterized=True,
    )


def plot_plane_panels(
    plane: PlotPlane,
    fields: list[str],
    config: dict,
    output: str | Path,
    *,
    fieldlines: bool = False,
    columns: int = 2,
) -> Path:
    rows = int(math.ceil(len(fields) / columns))
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(7.2 * columns, 5.2 * rows),
        constrained_layout=True,
        squeeze=False,
    )
    for panel, field in zip(axes.flat, fields):
        values = plane.level if field == "amr_level" else plane.values[field]
        label, cmap, norm = field_style(field, values, config)
        image = draw_scalar(panel, plane, values, cmap=cmap, norm=norm)
        panel.set_aspect(
            "equal" if config["slice"]["equal_physical_scale"] else "auto",
            adjustable="box",
        )
        if fieldlines and field != "amr_level":
            add_fieldlines(panel, plane)
        if field == "amr_level":
            levels = np.arange(int(np.nanmin(values)) + 0.5, int(np.nanmax(values)) + 0.5)
            if levels.size:
                panel.contour(plane.x, plane.y, values, levels=levels, colors="white", linewidths=0.55)
        panel.set_title(label)
        panel.set_xlabel(plane.x_name)
        panel.set_ylabel(plane.y_name)
        panel.grid(alpha=0.10, linewidth=0.35)
        figure.colorbar(image, ax=panel, label=label, shrink=0.9)
    for panel in axes.flat[len(fields) :]:
        panel.remove()
    figure.suptitle(f"{plane.fixed_label}, t={plane.time:.6g}")
    destination = Path(output).expanduser().resolve()
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=int(config["plot"]["dpi"]), bbox_inches="tight")
    plt.close(figure)
    return destination

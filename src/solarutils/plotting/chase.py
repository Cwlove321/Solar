"""Plotting helpers for CHASE spectral maps."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import numpy as np

from solarutils.obs.chase import chase_doppler_velocity_map, chase_map_at_velocity
from solarutils.plotting.maps import map_extent_arcsec


def plot_chase_channel(
    chase_map,
    *,
    extent: tuple[float, float, float, float] | None = None,
    ax=None,
    output: str | Path | None = None,
    cmap: str = "gray",
    vmin: float | None = None,
    vmax: float | None = None,
    major_tick: float = 100,
    title: str | None = None,
    figsize: tuple[float, float] = (5.8, 4.0),
    dpi: int = 180,
):
    """Plot one CHASE wavelength/velocity channel."""

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure
    if extent is None:
        extent = map_extent_arcsec(chase_map)
    if vmax is None:
        finite = np.asarray(chase_map.data, dtype=float)
        vmax = float(np.nanpercentile(finite, 99.5))
    image = ax.imshow(chase_map.data, origin="lower", extent=extent, cmap=cmap, vmin=vmin, vmax=vmax)
    if title is None:
        velocity = chase_map.meta.get("doppler", None)
        time = chase_map.date.value[-12:-4] if hasattr(chase_map, "date") else ""
        if velocity is None:
            title = f"CHASE Halpha {time} UT".strip()
        else:
            title = f"CHASE Halpha {velocity:g} km/s {time} UT".strip()
    ax.set_title(title, color="black", loc="center")
    _format_solar_axes(ax, major_tick=major_tick)
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")
    return image


def plot_chase_file_channel(
    path: str | Path,
    velocity_kms: float,
    *,
    extent: tuple[float, float, float, float] | None = None,
    output: str | Path | None = None,
    **kwargs,
):
    """Load and plot one CHASE file at a Doppler velocity channel."""

    smap = chase_map_at_velocity(path, velocity_kms, extent=extent)
    image = plot_chase_channel(smap, extent=extent, output=output, **kwargs)
    return smap, image


def plot_chase_doppler_velocity(
    velocity_map,
    *,
    extent: tuple[float, float, float, float] | None = None,
    ax=None,
    output: str | Path | None = None,
    vmax: float = 15.0,
    cmap: str = "RdBu_r",
    colorbar: bool = True,
    major_tick: float = 100,
    title: str | None = None,
    figsize: tuple[float, float] = (5.8, 4.0),
    dpi: int = 180,
):
    """Plot a CHASE Doppler velocity map in km/s."""

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure
    if extent is None:
        extent = map_extent_arcsec(velocity_map)
    image = ax.imshow(velocity_map.data, origin="lower", extent=extent, cmap=cmap, vmin=-vmax, vmax=vmax)
    if title is None:
        time = velocity_map.date.value[-12:-4] if hasattr(velocity_map, "date") else ""
        title = f"CHASE Doppler velocity {time} UT".strip()
    ax.set_title(title, color="black", loc="center")
    _format_solar_axes(ax, major_tick=major_tick)
    cbar = None
    if colorbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="4%", pad=0.08)
        cbar = fig.colorbar(image, cax=cax)
        cbar.set_label("Velocity (km/s)")
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")
    return {"figure": fig, "axis": ax, "image": image, "colorbar": cbar}


def plot_chase_file_doppler_velocity(
    path: str | Path,
    *,
    extent: tuple[float, float, float, float] | None = None,
    output: str | Path | None = None,
    **kwargs,
):
    """Calculate and plot a CHASE Doppler velocity map from one file."""

    velocity_map, line_center_map = chase_doppler_velocity_map(path, extent=extent)
    artists = plot_chase_doppler_velocity(velocity_map, extent=extent, output=output, **kwargs)
    return velocity_map, line_center_map, artists


def _format_solar_axes(ax, major_tick: float = 100) -> None:
    ax.set_xlabel("Solar-X (arcsec)")
    ax.set_ylabel("Solar-Y (arcsec)")
    ax.xaxis.set_major_locator(MultipleLocator(base=major_tick))
    ax.yaxis.set_major_locator(MultipleLocator(base=major_tick))
    ax.minorticks_on()
    ax.tick_params(which="major", direction="out", top=True, right=True)
    ax.tick_params(which="minor", direction="out", top=True, right=True)

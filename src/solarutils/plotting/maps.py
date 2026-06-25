"""Common plotting helpers for SunPy maps."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.patches import ConnectionPatch, Rectangle
from matplotlib.ticker import MultipleLocator
import numpy as np


def map_extent_arcsec(smap) -> tuple[float, float, float, float]:
    """Return a SunPy map extent as `(x_min, x_max, y_min, y_max)` in arcsec."""

    return (
        smap.bottom_left_coord.Tx.to_value("arcsec"),
        smap.top_right_coord.Tx.to_value("arcsec"),
        smap.bottom_left_coord.Ty.to_value("arcsec"),
        smap.top_right_coord.Ty.to_value("arcsec"),
    )


def plot_aia_map(
    smap,
    extent=None,
    ax=None,
    norm_min: float | None = None,
    norm_max: float | None = None,
    major_tick: float = 50,
):
    """Plot an AIA map with the map's native colormap and normalization."""

    if ax is None:
        _, ax = plt.subplots()
    if extent is None:
        extent = map_extent_arcsec(smap)

    norm = smap.plot_settings["norm"]
    norm.vmin = 1 if norm_min is None else norm_min
    if norm_max is not None:
        norm.vmax = norm_max

    image = ax.imshow(smap.data, norm=norm, extent=extent, origin="lower", cmap=smap.cmap)
    wave = int(smap.wavelength.value) if hasattr(smap, "wavelength") else ""
    time = smap.date.value[-12:-4]
    ax.set_title(f"AIA {wave} Å {time} UT".strip(), color="black", loc="center")
    _format_solar_axes(ax, major_tick=major_tick)
    return image


def plot_hmi_los(smap, extent=None, vmax: float = 500.0, ax=None, major_tick: float = 50):
    """Plot an HMI line-of-sight magnetogram."""

    if ax is None:
        _, ax = plt.subplots()
    if extent is None:
        extent = map_extent_arcsec(smap)

    image = ax.imshow(smap.data, vmax=vmax, vmin=-vmax, extent=extent, origin="lower", cmap=smap.cmap)
    time = smap.date.value[-12:-4]
    ax.set_title(f"HMI LOS {time} UT", color="black", loc="center")
    _format_solar_axes(ax, major_tick=major_tick)
    return image


def plot_hmi_on_aia(
    aia_map,
    hmi_map,
    *,
    extent: tuple[float, float, float, float] | None = None,
    ax=None,
    output: str | Path | None = None,
    aia_norm_min: float | None = None,
    aia_norm_max: float | None = None,
    hmi_threshold: float | None = 30.0,
    hmi_vmax: float = 100.0,
    hmi_alpha: float = 0.5,
    hmi_cmap: str = "RdBu_r",
    contour_levels: tuple[float, ...] | list[float] | None = None,
    contour_kwargs: dict | None = None,
    boxes: list[tuple[float, float, float, float]] | None = None,
    colorbar: bool = True,
    colorbar_label: str = "Magnetic Field (Gauss)",
    shade_threshold_on_colorbar: bool = True,
    major_tick: float = 100,
    title: str | None = None,
    figsize: tuple[float, float] = (6.8, 5.2),
    dpi: int = 180,
):
    """Overlay an HMI LOS magnetogram on an AIA map in solar-X/Y coordinates.

    The maps are assumed to share the same solar-coordinate footprint. This is
    the maintained version of the common notebook plot where AIA is the base
    image and weak HMI fields are masked to transparent before overlaying.
    """

    if ax is None:
        fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    else:
        fig = ax.figure
    if extent is None:
        extent = map_extent_arcsec(aia_map)

    aia_image = plot_aia_map(
        aia_map,
        extent=extent,
        ax=ax,
        norm_min=aia_norm_min,
        norm_max=aia_norm_max,
        major_tick=major_tick,
    )

    masked_hmi = mask_weak_hmi_field(hmi_map.data, threshold=hmi_threshold)
    cmap = mpl.colormaps[hmi_cmap].copy()
    cmap.set_bad(alpha=0.0)
    hmi_image = ax.imshow(
        masked_hmi,
        origin="lower",
        cmap=cmap,
        extent=extent,
        vmin=-hmi_vmax,
        vmax=hmi_vmax,
        alpha=hmi_alpha,
    )

    contours = None
    if contour_levels is not None:
        contour_style = {"alpha": 1.0, "cmap": "RdBu_r", "linewidths": 1.0}
        if contour_kwargs:
            contour_style.update(contour_kwargs)
        contours = ax.contour(
            hmi_map.data,
            levels=contour_levels,
            origin="lower",
            extent=extent,
            **contour_style,
        )

    box_artists = []
    if boxes:
        for box_extent in boxes:
            box_artists.append(draw_submap_box(ax, box_extent, edgecolor="white", linewidth=2.0))

    cbar = None
    if colorbar:
        from mpl_toolkits.axes_grid1 import make_axes_locatable

        divider = make_axes_locatable(ax)
        cax = divider.append_axes("right", size="5%", pad=0.1)
        cbar = fig.colorbar(hmi_image, cax=cax, extend="both")
        cbar.set_label(colorbar_label)
        if hmi_threshold is not None and shade_threshold_on_colorbar:
            threshold = abs(float(hmi_threshold))
            cbar.ax.axhspan(
                -threshold,
                threshold,
                facecolor="gray",
                alpha=0.3,
                hatch="////",
                edgecolor="black",
            )

    if title is None:
        wave = int(aia_map.wavelength.value) if hasattr(aia_map, "wavelength") else ""
        time = aia_map.date.value[-12:-4]
        if hmi_threshold is None:
            title = f"HMI overlaid on AIA {wave} Å {time} UT".strip()
        else:
            title = f"HMI (|B| > {hmi_threshold:g} G) overlaid on AIA {wave} Å {time} UT".strip()
    ax.set_title(title, color="black", loc="center")

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")

    return {
        "figure": fig,
        "axis": ax,
        "aia": aia_image,
        "hmi": hmi_image,
        "colorbar": cbar,
        "contours": contours,
        "boxes": box_artists,
    }


def mask_weak_hmi_field(data, threshold: float | None = 30.0):
    """Return HMI data with weak fields masked to NaN for transparent plotting."""

    masked = np.asarray(data, dtype=float).copy()
    if threshold is not None:
        masked[np.abs(masked) < abs(float(threshold))] = np.nan
    return masked


def draw_submap_box(ax, extent: tuple[float, float, float, float], **kwargs):
    """Draw the child-map footprint on a parent map axis."""

    style = {
        "edgecolor": "cyan",
        "facecolor": "none",
        "linewidth": 1.4,
        "linestyle": "--",
        "alpha": 0.95,
    }
    style.update(kwargs)
    x_min, x_max, y_min, y_max = extent
    box = Rectangle((x_min, y_min), x_max - x_min, y_max - y_min, **style)
    ax.add_patch(box)
    return box


def connect_submap_axes(
    parent_ax,
    child_ax,
    extent: tuple[float, float, float, float],
    color: str = "cyan",
    linewidth: float = 1.0,
    alpha: float = 0.9,
    linestyle: str = "--",
):
    """Connect a child-map footprint on the parent axis to the child axis."""

    x_min, x_max, y_min, y_max = extent
    connectors = []
    for parent_xy, child_xy in (
        ((x_max, y_max), (0, 1)),
        ((x_max, y_min), (0, 0)),
    ):
        connector = ConnectionPatch(
            xyA=parent_xy,
            coordsA=parent_ax.transData,
            xyB=child_xy,
            coordsB=child_ax.transAxes,
            color=color,
            linewidth=linewidth,
            linestyle=linestyle,
            alpha=alpha,
            clip_on=False,
        )
        child_ax.add_artist(connector)
        connectors.append(connector)
    return connectors


def plot_aia_context_submap(
    smap,
    submap,
    submap_extent: tuple[float, float, float, float],
    *,
    output: str | Path | None = None,
    norm_max: float | None = None,
    figsize: tuple[float, float] = (8.2, 4.3),
    dpi: int = 180,
):
    """Plot an AIA context parent map beside a child submap."""

    fig, (parent_ax, child_ax) = plt.subplots(1, 2, figsize=figsize, dpi=dpi)
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.14, top=0.88, wspace=0.22)

    plot_aia_map(smap, ax=parent_ax, norm_max=norm_max, major_tick=500)
    draw_submap_box(parent_ax, submap_extent)

    plot_aia_map(submap, ax=child_ax, norm_max=norm_max, major_tick=100)
    connect_submap_axes(parent_ax, child_ax, submap_extent)

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")
    return fig, (parent_ax, child_ax)


def _format_solar_axes(ax, major_tick: float = 50) -> None:
    ax.set_xlabel("Solar-X (arcsec)")
    ax.set_ylabel("Solar-Y (arcsec)")
    ax.yaxis.set_major_locator(MultipleLocator(base=major_tick))
    ax.xaxis.set_major_locator(MultipleLocator(base=major_tick))
    ax.minorticks_on()
    ax.tick_params(which="major", direction="out", top=True, right=True)
    ax.tick_params(which="minor", direction="out", top=True, right=True)

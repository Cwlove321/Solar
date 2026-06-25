"""Plotting helpers."""

from solarutils.plotting.chase import (
    plot_chase_channel,
    plot_chase_doppler_velocity,
    plot_chase_file_channel,
    plot_chase_file_doppler_velocity,
)
from solarutils.plotting.maps import (
    connect_submap_axes,
    draw_submap_box,
    map_extent_arcsec,
    mask_weak_hmi_field,
    plot_aia_context_submap,
    plot_aia_map,
    plot_hmi_on_aia,
    plot_hmi_los,
)
from solarutils.plotting.movie import MOVIE_QUALITY_PRESETS, make_movie_from_images
from solarutils.plotting.vectors import plot_signed_vectors

__all__ = [
    "connect_submap_axes",
    "draw_submap_box",
    "map_extent_arcsec",
    "mask_weak_hmi_field",
    "make_movie_from_images",
    "MOVIE_QUALITY_PRESETS",
    "plot_chase_channel",
    "plot_chase_doppler_velocity",
    "plot_chase_file_channel",
    "plot_chase_file_doppler_velocity",
    "plot_aia_context_submap",
    "plot_aia_map",
    "plot_hmi_on_aia",
    "plot_hmi_los",
    "plot_signed_vectors",
]

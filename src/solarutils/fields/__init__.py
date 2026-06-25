"""Magnetic-field analysis utilities."""

from solarutils.fields.flux import (
    VALID_FLUX_MODES,
    magnetic_flux_for_extent,
    magnetic_flux_from_data,
    magnetic_flux_pixel_area_cm2,
    magnetic_flux_series,
)
from solarutils.fields.helicity import e_flux, h_flux, h_transpose, volume_integral
from solarutils.fields.vector import curl

__all__ = [
    "VALID_FLUX_MODES",
    "curl",
    "e_flux",
    "h_flux",
    "h_transpose",
    "magnetic_flux_for_extent",
    "magnetic_flux_from_data",
    "magnetic_flux_pixel_area_cm2",
    "magnetic_flux_series",
    "volume_integral",
]

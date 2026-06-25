"""Magnetic flux helpers for photospheric magnetograms."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Literal

import astropy.units as u
import numpy as np

from solarutils.obs.maps import get_submap


FluxMode = Literal["unsigned", "positive", "negative"]
VALID_FLUX_MODES = ("unsigned", "positive", "negative")


def magnetic_flux_from_data(data, pixel_area_cm2: float, mode: FluxMode = "unsigned") -> float:
    """Calculate magnetic flux from a magnetogram array.

    Parameters
    ----------
    data
        Magnetic field in Gauss.
    pixel_area_cm2
        Pixel area in cm^2.
    mode
        `"unsigned"` integrates `abs(B)`, `"positive"` integrates only `B > 0`,
        and `"negative"` integrates only `B < 0` while keeping the negative sign.
    """

    if mode not in VALID_FLUX_MODES:
        raise ValueError(f"mode must be one of {VALID_FLUX_MODES}")

    data = np.asarray(data, dtype=float)
    finite = np.isfinite(data)
    if mode == "unsigned":
        selected = np.abs(data[finite])
    elif mode == "positive":
        selected = data[finite & (data > 0)]
    else:
        selected = data[finite & (data < 0)]
    return float(np.sum(selected) * pixel_area_cm2)

def magnetic_flux_pixel_area_cm2(smap) -> float:
    """Estimate a SunPy map pixel area in cm^2 from plate scale and distance."""

    try:
        cdelt1 = abs(smap.scale.axis1.to_value(u.rad / u.pix))
        cdelt2 = abs(smap.scale.axis2.to_value(u.rad / u.pix))
    except Exception:
        cdelt1 = abs(float(smap.meta["cdelt1"])) * u.arcsec.to(u.rad)
        cdelt2 = abs(float(smap.meta["cdelt2"])) * u.arcsec.to(u.rad)

    try:
        distance_cm = smap.dsun.to_value(u.cm)
    except Exception:
        distance_cm = float(smap.meta["dsun_obs"]) * u.m.to(u.cm)

    return float((distance_cm * cdelt1) * (distance_cm * cdelt2))


def magnetic_flux_for_extent(
    smap,
    extent: tuple[float, float, float, float] | None = None,
    *,
    modes: Iterable[FluxMode] = VALID_FLUX_MODES,
    pixel_area_cm2: float | None = None,
) -> dict[str, float]:
    """Calculate magnetic flux for one map extent.

    `extent` is `(x_min, x_max, y_min, y_max)` in arcsec. Pass `None` to use the
    whole input map.
    """

    region = smap if extent is None else get_submap(smap, extent)
    area = magnetic_flux_pixel_area_cm2(region) if pixel_area_cm2 is None else float(pixel_area_cm2)
    return {mode: magnetic_flux_from_data(region.data, area, mode=mode) for mode in modes}


def magnetic_flux_series(
    paths: Iterable[str | Path],
    extents: Iterable[tuple[float, float, float, float]],
    *,
    modes: Iterable[FluxMode] = VALID_FLUX_MODES,
    region_labels: Iterable[str] | None = None,
    output: str | Path | None = None,
) -> dict[str, np.ndarray]:
    """Calculate magnetic flux for a time sequence of magnetogram submaps.

    The optional `.npz` output contains `time_iso`, `time_mjd`, `paths`,
    `extents`, `region_labels`, `modes`, and `flux_mx` with shape
    `(time, region, mode)`.
    """

    from sunpy.map import Map

    paths = [Path(path) for path in paths]
    extents = [tuple(float(value) for value in extent) for extent in extents]
    modes = tuple(modes)
    _validate_flux_modes(modes)

    if not paths:
        raise ValueError("paths must contain at least one magnetogram")
    if not extents:
        raise ValueError("extents must contain at least one region")

    if region_labels is None:
        labels = [f"R{index + 1}" for index in range(len(extents))]
    else:
        labels = list(region_labels)
        if len(labels) != len(extents):
            raise ValueError("region_labels must have the same length as extents")

    flux = np.empty((len(paths), len(extents), len(modes)), dtype=float)
    time_iso = []
    time_mjd = []
    for time_index, path in enumerate(paths):
        smap = Map(path)
        time_iso.append(smap.date.isot)
        time_mjd.append(smap.date.mjd)
        for region_index, extent in enumerate(extents):
            values = magnetic_flux_for_extent(smap, extent, modes=modes)
            for mode_index, mode in enumerate(modes):
                flux[time_index, region_index, mode_index] = values[mode]

    result = {
        "time_iso": np.asarray(time_iso),
        "time_mjd": np.asarray(time_mjd, dtype=float),
        "paths": np.asarray([str(path) for path in paths]),
        "extents": np.asarray(extents, dtype=float),
        "region_labels": np.asarray(labels),
        "modes": np.asarray(modes),
        "flux_mx": flux,
    }

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(output, **result)

    return result


def _validate_flux_modes(modes: Iterable[str]) -> None:
    for mode in modes:
        if mode not in VALID_FLUX_MODES:
            raise ValueError(f"mode must be one of {VALID_FLUX_MODES}")

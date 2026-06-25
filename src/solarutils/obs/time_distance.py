"""Time-distance helpers for processed image sequences."""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from scipy.ndimage import map_coordinates


def sample_line_data(
    data,
    start_pixel: tuple[float, float],
    end_pixel: tuple[float, float],
    n_samples: int,
    *,
    order: int = 1,
    fill_value: float = np.nan,
) -> np.ndarray:
    """Sample a 2D image along a line in pixel coordinates.

    `start_pixel` and `end_pixel` are `(x, y)`. The returned array has length
    `n_samples` and includes both endpoints.
    """

    if n_samples < 2:
        raise ValueError("n_samples must be at least 2")

    x0, y0 = start_pixel
    x1, y1 = end_pixel
    x = np.linspace(float(x0), float(x1), n_samples)
    y = np.linspace(float(y0), float(y1), n_samples)
    return map_coordinates(
        np.asarray(data, dtype=float),
        [y, x],
        order=order,
        mode="constant",
        cval=fill_value,
    )


def sample_map_line(
    smap,
    start: tuple[float, float],
    end: tuple[float, float],
    n_samples: int,
    *,
    order: int = 1,
    fill_value: float = np.nan,
) -> np.ndarray:
    """Sample a SunPy map along a solar-X/Y line in arcsec."""

    start_coord = _skycoord_from_arcsec(start, smap)
    end_coord = _skycoord_from_arcsec(end, smap)
    x0, y0 = smap.wcs.world_to_pixel(start_coord)
    x1, y1 = smap.wcs.world_to_pixel(end_coord)
    return sample_line_data(
        smap.data,
        (float(x0), float(y0)),
        (float(x1), float(y1)),
        n_samples,
        order=order,
        fill_value=fill_value,
    )


def time_distance_series(
    paths: Iterable[str | Path],
    start: tuple[float, float],
    end: tuple[float, float],
    n_samples: int,
    *,
    order: int = 1,
    fill_value: float = np.nan,
    output: str | Path | None = None,
) -> dict[str, np.ndarray | float]:
    """Sample a processed map sequence along one fixed solar-X/Y line.

    The optional `.npz` output contains `time_iso`, `time_mjd`, `paths`,
    `start_arcsec`, `end_arcsec`, `distance_km`, `total_distance_km`, and
    `data` with shape `(time, distance)`.
    """

    from sunpy.map import Map

    paths = [Path(path) for path in paths]
    if not paths:
        raise ValueError("paths must contain at least one map")

    start = (float(start[0]), float(start[1]))
    end = (float(end[0]), float(end[1]))
    data = np.empty((len(paths), n_samples), dtype=float)
    time_iso = []
    time_mjd = []
    total_distance_km = None
    distance_km = None
    for index, path in enumerate(paths):
        smap = Map(path)
        if total_distance_km is None:
            total_distance_km = line_distance_km(smap, start, end)
            distance_km = np.linspace(0.0, total_distance_km, n_samples)
        data[index] = sample_map_line(
            smap,
            start,
            end,
            n_samples,
            order=order,
            fill_value=fill_value,
        )
        time_iso.append(smap.date.isot)
        time_mjd.append(smap.date.mjd)

    result = {
        "time_iso": np.asarray(time_iso),
        "time_mjd": np.asarray(time_mjd, dtype=float),
        "paths": np.asarray([str(path) for path in paths]),
        "start_arcsec": np.asarray(start, dtype=float),
        "end_arcsec": np.asarray(end, dtype=float),
        "distance_km": distance_km,
        "total_distance_km": total_distance_km,
        "data": data,
    }

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(output, **result)

    return result


def _skycoord_from_arcsec(point: tuple[float, float], smap) -> SkyCoord:
    return SkyCoord(point[0] * u.arcsec, point[1] * u.arcsec, frame=smap.coordinate_frame)


def line_distance_km(smap, start: tuple[float, float], end: tuple[float, float]) -> float:
    """Return the projected line length in km for two solar-X/Y arcsec points."""

    angular_distance = np.hypot(float(end[0]) - float(start[0]), float(end[1]) - float(start[1])) * u.arcsec
    return float((angular_distance.to_value(u.rad) * smap.dsun.to_value(u.km)))

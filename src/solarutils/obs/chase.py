"""CHASE spectral-map helpers."""

from __future__ import annotations

from pathlib import Path

import astropy.units as u
import numpy as np
from sunpy.coordinates import frames
from scipy.optimize import curve_fit
from sunpy.map import Map

from solarutils.io import read_chase
from solarutils.obs.maps import get_submap


SPEED_OF_LIGHT_KMS = 299792.458
SOLAR_RADIUS_KM = 695700.0
SNODGRASS_ROTATION_COEFFS_DEG_PER_DAY = (14.713, -2.396, -1.787)


def chase_wavelength_axis(header, n_wavelength: int | None = None) -> np.ndarray:
    """Return the CHASE wavelength axis in Angstrom."""

    if n_wavelength is None:
        n_wavelength = int(header["NAXIS3"])
    return float(header["CRVAL3"]) + np.arange(n_wavelength, dtype=float) * float(header["CDELT3"])


def chase_reference_line_center(
    data,
    header,
    wavelength=None,
    *,
    center_extent: tuple[float, float, float, float] = (-100.0, 100.0, -100.0, 100.0),
    window: tuple[float, float] | None = None,
) -> dict[str, float | np.ndarray]:
    """Fit the mean disk-center absorption profile and return the reference line center.

    The default reference region is the solar disk center plus/minus 100 arcsec.
    A Gaussian is fitted to `continuum - intensity`, so the fitted center is the
    absorption-line center.
    """

    wavelength = chase_wavelength_axis(header, data.shape[0]) if wavelength is None else np.asarray(wavelength)
    reference_map = Map(data[0], header)
    center_map = get_submap(reference_map, center_extent)
    x0, x1, y0, y1 = map(int, _map_pixel_bounds(reference_map, center_map))
    mean_profile = np.nanmean(data[:, y0:y1, x0:x1], axis=(1, 2))
    fitted_center = fit_absorption_gaussian_center(wavelength, mean_profile, window=window)
    return {
        "line_center": float(fitted_center),
        "wavelength": wavelength,
        "mean_profile": mean_profile,
        "center_extent": np.asarray(center_extent, dtype=float),
    }


def fit_absorption_gaussian_center(
    wavelength,
    profile,
    *,
    window: tuple[float, float] | None = None,
) -> float:
    """Fit a Gaussian absorption profile and return its center wavelength."""

    wavelength = np.asarray(wavelength, dtype=float)
    profile = np.asarray(profile, dtype=float)
    finite = np.isfinite(wavelength) & np.isfinite(profile)
    if window is not None:
        finite &= (wavelength >= min(window)) & (wavelength <= max(window))
    x = wavelength[finite]
    y = profile[finite]
    if x.size < 5:
        raise ValueError("At least five spectral samples are required for Gaussian fitting")

    continuum0 = float(np.nanmax(y))
    depth0 = max(float(continuum0 - np.nanmin(y)), 1.0)
    center0 = float(x[np.nanargmin(y)])
    sigma0 = max((x.max() - x.min()) / 8.0, abs(np.nanmedian(np.diff(x))) * 2.0)
    bounds = (
        [np.nanmin(y) - 10 * depth0, 0.0, x.min(), abs(np.nanmedian(np.diff(x))) / 4.0],
        [np.nanmax(y) + 10 * depth0, np.inf, x.max(), x.max() - x.min()],
    )
    try:
        params, _ = curve_fit(
            _absorption_gaussian,
            x,
            y,
            p0=(continuum0, depth0, center0, sigma0),
            bounds=bounds,
            maxfev=10000,
        )
        return float(params[2])
    except Exception:
        return centroid_line_center(x, y)


def centroid_line_center(wavelength, profile) -> float:
    """Return an absorption-line center from the intensity centroid."""

    wavelength = np.asarray(wavelength, dtype=float)
    profile = np.asarray(profile, dtype=float)
    continuum = np.nanmax(profile)
    weights = np.clip(continuum - profile, 0.0, None)
    total = np.nansum(weights)
    if not np.isfinite(total) or total <= 0:
        return float("nan")
    return float(np.nansum(wavelength * weights) / total)


def chase_channel_index_for_velocity(wavelength, line_center: float, velocity_kms: float) -> int:
    """Return the spectral index nearest to a Doppler velocity offset."""

    target = doppler_shifted_wavelength(line_center, velocity_kms)
    return nearest_wavelength_index(wavelength, target)


def nearest_wavelength_index(wavelength, target: float) -> int:
    """Return the index of the wavelength sample nearest to `target`."""

    wavelength = np.asarray(wavelength, dtype=float)
    return int(np.nanargmin(np.abs(wavelength - float(target))))


def doppler_shifted_wavelength(line_center: float, velocity_kms: float) -> float:
    """Convert a velocity offset to wavelength with the non-relativistic formula."""

    return float(line_center) * (1.0 + float(velocity_kms) / SPEED_OF_LIGHT_KMS)


def chase_map_at_velocity(
    path: str | Path,
    velocity_kms: float,
    *,
    extent: tuple[float, float, float, float] | None = None,
    reference_line_center: float | None = None,
    center_extent: tuple[float, float, float, float] = (-100.0, 100.0, -100.0, 100.0),
    line_window: tuple[float, float] | None = None,
):
    """Load one CHASE Doppler channel as a SunPy map.

    `velocity_kms=0` means the fitted reference line center. Positive velocity
    selects the red side, and negative velocity selects the blue side.
    """

    data, header, wavelength = read_chase(path)
    if reference_line_center is None:
        reference = chase_reference_line_center(
            data,
            header,
            wavelength,
            center_extent=center_extent,
            window=line_window,
        )
        reference_line_center = float(reference["line_center"])

    index = chase_channel_index_for_velocity(wavelength, reference_line_center, velocity_kms)
    smap = Map(data[index], header)
    if extent is not None:
        smap = get_submap(smap, extent)
    smap.meta["waveunit"] = "Angstrom"
    smap.meta["wavelnth"] = float(wavelength[index])
    smap.meta["doppler"] = float(velocity_kms)
    smap.meta["linecen"] = float(reference_line_center)
    return smap


def chase_doppler_velocity_map(
    path: str | Path,
    *,
    extent: tuple[float, float, float, float] | None = None,
    reference_line_center: float | None = None,
    center_extent: tuple[float, float, float, float] = (-100.0, 100.0, -100.0, 100.0),
    line_window: tuple[float, float] | None = None,
    method: str = "centroid",
    clip_absorption: bool = True,
    remove_solar_rotation: bool = False,
    rotation_reference_extent: tuple[float, float, float, float] | None = None,
):
    """Calculate a CHASE Doppler velocity map in km/s.

    The reference line center is fitted from the disk-center mean spectrum with
    a Gaussian. The first implemented pixel method is `centroid`.
    """

    if method != "centroid":
        raise ValueError("Only method='centroid' is currently implemented")

    data, header, wavelength = read_chase(path)
    if reference_line_center is None:
        reference = chase_reference_line_center(
            data,
            header,
            wavelength,
            center_extent=center_extent,
            window=line_window,
        )
        reference_line_center = float(reference["line_center"])

    if extent is None:
        spatial = data
        out_map = Map(data[nearest_wavelength_index(wavelength, reference_line_center)], header)
    else:
        reference_map = Map(data[nearest_wavelength_index(wavelength, reference_line_center)], header)
        submap = get_submap(reference_map, extent)
        x0, x1, y0, y1 = map(int, _map_pixel_bounds(reference_map, submap))
        spatial = data[:, y0:y1, x0:x1]
        out_map = Map(spatial[nearest_wavelength_index(wavelength, reference_line_center)], submap.meta.copy())

    selected_wavelength, selected_data = _spectral_window(wavelength, spatial, line_window, reference_line_center)
    line_center_map = centroid_line_center_map(selected_wavelength, selected_data, clip_absorption=clip_absorption)
    velocity = (line_center_map - reference_line_center) / reference_line_center * SPEED_OF_LIGHT_KMS
    if remove_solar_rotation:
        reference_map = Map(data[nearest_wavelength_index(wavelength, reference_line_center)], header)
        rotation_reference_extent = center_extent if rotation_reference_extent is None else rotation_reference_extent
        rotation_velocity = solar_rotation_los_velocity_map(out_map)
        reference_rotation = np.nanmean(
            solar_rotation_los_velocity_map(get_submap(reference_map, rotation_reference_extent))
        )
        velocity = velocity - (rotation_velocity - reference_rotation)

    meta = out_map.meta.copy()
    meta["bunit"] = "km / s"
    meta["linecen"] = float(reference_line_center)
    meta["rotcorr"] = bool(remove_solar_rotation)
    velocity_map = Map(velocity.astype(np.float32), meta)
    return velocity_map, line_center_map


def centroid_line_center_map(wavelength, data, *, clip_absorption: bool = True) -> np.ndarray:
    """Return a per-pixel absorption centroid map from a spectral cube."""

    wavelength = np.asarray(wavelength, dtype=float)
    data = np.asarray(data, dtype=float)
    continuum = np.nanmax(data, axis=0)
    weights = continuum[None, :, :] - data
    if clip_absorption:
        weights = np.clip(weights, 0.0, None)
    total = np.nansum(weights, axis=0)
    numerator = np.nansum(weights * wavelength[:, None, None], axis=0)
    center = np.full(total.shape, np.nan, dtype=float)
    valid = np.isfinite(total) & (total > 0)
    center[valid] = numerator[valid] / total[valid]
    return center


def solar_rotation_los_velocity_map(
    smap,
    *,
    coeffs: tuple[float, float, float] = SNODGRASS_ROTATION_COEFFS_DEG_PER_DAY,
    radius_km: float = SOLAR_RADIUS_KM,
) -> np.ndarray:
    """Estimate the differential-rotation LOS velocity for every map pixel.

    Positive velocity is defined as motion away from the observer, matching the
    redshift sign convention used by the Doppler conversion.
    """

    yy, xx = np.indices(smap.data.shape, dtype=float)
    coords = smap.wcs.pixel_to_world(xx, yy)
    hgs = coords.transform_to(frames.HeliographicStonyhurst(obstime=smap.date))
    observer = smap.observer_coordinate.transform_to(frames.HeliographicStonyhurst(obstime=smap.date))

    lon = hgs.lon.to_value(u.rad)
    lat = hgs.lat.to_value(u.rad)
    observer_lat = observer.lat.to_value(u.rad)
    omega = differential_rotation_rate(lat, coeffs=coeffs)
    velocity_toward_observer = -omega * radius_km * np.cos(lat) * np.sin(lon) * np.cos(observer_lat)
    velocity_away_observer = -velocity_toward_observer
    return velocity_away_observer


def differential_rotation_rate(
    latitude_rad,
    *,
    coeffs: tuple[float, float, float] = SNODGRASS_ROTATION_COEFFS_DEG_PER_DAY,
) -> np.ndarray:
    """Return the Snodgrass differential rotation rate in rad/s."""

    latitude_rad = np.asarray(latitude_rad, dtype=float)
    sin2 = np.sin(latitude_rad) ** 2
    a, b, c = coeffs
    omega_deg_per_day = a + b * sin2 + c * sin2**2
    return omega_deg_per_day * np.pi / 180.0 / 86400.0


def _spectral_window(wavelength, data, window, reference_line_center: float):
    wavelength = np.asarray(wavelength, dtype=float)
    if window is None:
        half_width = 2.0
        mask = np.abs(wavelength - reference_line_center) <= half_width
    else:
        mask = (wavelength >= min(window)) & (wavelength <= max(window))
    if np.count_nonzero(mask) < 5:
        raise ValueError("Spectral window must contain at least five samples")
    return wavelength[mask], data[mask]


def _absorption_gaussian(x, continuum, depth, center, sigma):
    return continuum - depth * np.exp(-0.5 * ((x - center) / sigma) ** 2)


def _map_pixel_bounds(parent_map, child_map) -> tuple[int, int, int, int]:
    bottom_left = child_map.bottom_left_coord
    top_right = child_map.top_right_coord
    x0, y0 = parent_map.world_to_pixel(bottom_left)
    x1, y1 = parent_map.world_to_pixel(top_right)
    x_min = max(0, int(np.floor(min(x0.value, x1.value))))
    x_max = min(parent_map.data.shape[1], int(np.ceil(max(x0.value, x1.value))) + 1)
    y_min = max(0, int(np.floor(min(y0.value, y1.value))))
    y_max = min(parent_map.data.shape[0], int(np.ceil(max(y0.value, y1.value))) + 1)
    return x_min, x_max, y_min, y_max

"""SunPy map helpers for common SDO/AIA/HMI workflows."""

from __future__ import annotations

from pathlib import Path

import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from scipy import ndimage
from scipy.ndimage import map_coordinates, shift
from sunpy.map import Map
from sunpy.physics.differential_rotation import solar_rotate_coordinate


AIA_EUV_WAVELENGTHS = (94, 131, 171, 193, 211, 304, 335)
DEFAULT_AIA_PSF_DIR = Path("data/calibration/sdo/aia/psf")


def aia_prep(
    aia_map,
    pointing_table=None,
    psf=None,
    *,
    correct_pointing: bool = True,
    deconvolve_psf: bool = False,
    register_map: bool = True,
    respike_spikes=None,
    deconvolve_iterations: int = 25,
):
    """Prepare an AIA level-1 map toward level 1.5.

    The processing order follows the common AIA calibration workflow: update
    pointing, optionally restore despiked pixels, optionally deconvolve the PSF,
    and register the map onto the standard AIA plate scale/orientation.

    Parameters
    ----------
    aia_map
        Input SunPy AIA map.
    pointing_table
        Pointing table returned by `aiapy.calibrate.utils.get_pointing_table`.
        Passing it explicitly is recommended when preparing many maps.
    psf
        Precomputed PSF array. If omitted and `deconvolve_psf=True`, the PSF is
        calculated from the map wavelength.
    correct_pointing, deconvolve_psf, register_map
        Switches for each processing step. They are useful for debugging and
        lightweight tests.
    respike_spikes
        Optional `(positions, values)` tuple to pass to `aiapy.calibrate.respike`.
        Leave as `None` to skip respiking.
    deconvolve_iterations
        Number of Richardson-Lucy iterations used by `aiapy.psf.deconvolve`.
    """

    calibrate, psf_tools = _aiapy_imports()
    prepared_map = aia_map

    if correct_pointing:
        prepared_map = calibrate.update_pointing(prepared_map, pointing_table=pointing_table)

    if respike_spikes is not None:
        prepared_map = calibrate.respike(prepared_map, spikes=respike_spikes)

    if deconvolve_psf:
        if psf is None:
            psf = psf_tools.calculate_psf(prepared_map.wavelength)
        else:
            psf = _as_psf_array(psf)
        prepared_map = psf_tools.deconvolve(prepared_map, psf=psf, iterations=deconvolve_iterations)

    if register_map:
        prepared_map = calibrate.register(prepared_map)

    return prepared_map


def aia_psf_path(wavelength: int, directory: str | Path = DEFAULT_AIA_PSF_DIR) -> Path:
    """Return the local cache path for an AIA PSF."""

    return Path(directory) / f"aia_{int(wavelength):04d}_psf.npz"


def cache_aia_psfs(
    wavelengths=AIA_EUV_WAVELENGTHS,
    directory: str | Path = DEFAULT_AIA_PSF_DIR,
    *,
    overwrite: bool = False,
    dtype=np.float32,
) -> dict[int, Path]:
    """Calculate and cache AIA PSFs for later preprocessing."""

    _, psf_tools = _aiapy_imports()
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    cached = {}
    for wavelength in wavelengths:
        wavelength = int(wavelength)
        path = aia_psf_path(wavelength, directory)
        if path.exists() and not overwrite:
            cached[wavelength] = path
            continue

        psf = psf_tools.calculate_psf(wavelength * u.AA).astype(dtype, copy=False)
        np.savez_compressed(path, psf=psf, wavelength_angstrom=wavelength)
        cached[wavelength] = path
    return cached


def load_aia_psf(wavelength: int, directory: str | Path = DEFAULT_AIA_PSF_DIR):
    """Load a cached AIA PSF array."""

    path = aia_psf_path(wavelength, directory)
    if not path.exists():
        raise FileNotFoundError(f"AIA PSF cache file not found: {path}")
    return _as_psf_array(path)


def deconvolve_submap_with_padding(
    aia_map,
    extent: tuple[float, float, float, float],
    psf=None,
    *,
    padding: float = 200.0,
    iterations: int = 25,
):
    """Locally deconvolve a padded AIA submap and return the target extent.

    This is a faster approximation to full-frame PSF deconvolution followed by
    cropping. The input `extent` and `padding` are in arcsec. `psf` is optional;
    if omitted, it is calculated from the map wavelength.
    """

    _, psf_tools = _aiapy_imports()
    if psf is None:
        psf = psf_tools.calculate_psf(aia_map.wavelength)
    else:
        psf = _as_psf_array(psf)

    padded_extent = pad_extent(extent, padding, aia_map)
    padded_map = ensure_even_shape(get_submap(aia_map, padded_extent))
    local_psf = crop_psf_to_shape(psf, padded_map.data.shape)
    deconvolved = psf_tools.deconvolve(padded_map, psf=local_psf, iterations=iterations)
    return get_submap(deconvolved, extent)


def pad_extent(
    extent: tuple[float, float, float, float],
    padding: float,
    smap,
) -> tuple[float, float, float, float]:
    """Pad an arcsec extent while keeping it inside a map footprint."""

    x_min, x_max, y_min, y_max = extent
    return (
        max(x_min - padding, smap.bottom_left_coord.Tx.to_value("arcsec")),
        min(x_max + padding, smap.top_right_coord.Tx.to_value("arcsec")),
        max(y_min - padding, smap.bottom_left_coord.Ty.to_value("arcsec")),
        min(y_max + padding, smap.top_right_coord.Ty.to_value("arcsec")),
    )


def ensure_even_shape(smap):
    """Trim a map to even pixel dimensions for FFT-based deconvolution."""

    ny = smap.data.shape[0] // 2 * 2
    nx = smap.data.shape[1] // 2 * 2
    if (ny, nx) == smap.data.shape:
        return smap
    return Map(smap.data[:ny, :nx], smap.meta.copy())


def crop_psf_to_shape(psf, shape: tuple[int, int], *, normalize: bool = True):
    """Center-crop a PSF to a target `(ny, nx)` shape.

    This is useful for local/submap deconvolution. The target shape should have
    even dimensions because `aiapy.psf.deconvolve` uses `rfft2/irfft2`.
    """

    psf = _as_psf_array(psf)
    ny, nx = shape
    if ny % 2 or nx % 2:
        raise ValueError("target PSF shape must have even dimensions")
    if ny > psf.shape[0] or nx > psf.shape[1]:
        raise ValueError("target PSF shape cannot be larger than the source PSF")

    y0 = (psf.shape[0] - ny) // 2
    x0 = (psf.shape[1] - nx) // 2
    cropped = psf[y0 : y0 + ny, x0 : x0 + nx].copy()
    if normalize:
        total = cropped.sum()
        if total != 0:
            cropped /= total
    return cropped


def get_submap(smap, extent: tuple[float, float, float, float]):
    """Return a submap from an `(x_min, x_max, y_min, y_max)` arcsec extent."""

    bottom_left = SkyCoord(extent[0] * u.arcsec, extent[2] * u.arcsec, frame=smap.coordinate_frame)
    top_right = SkyCoord(extent[1] * u.arcsec, extent[3] * u.arcsec, frame=smap.coordinate_frame)
    return smap.submap(bottom_left, top_right=top_right)


def sharpen_map(smap, *, sigma: float = 1.0, weight: float = 1.5, dtype=np.float64):
    """Sharpen a SunPy map with an unsharp-mask filter.

    The formula follows the older case notebooks:
    `data + (data - gaussian_filter(data, sigma)) * weight`.
    """

    data = smap.data.astype(dtype, copy=False)
    blurred = ndimage.gaussian_filter(data, sigma=sigma)
    sharpened_data = data + (data - blurred) * weight
    return Map(sharpened_data, smap.meta.copy())


def multiscale_sharpen_map(
    smap,
    *,
    sigmas: tuple[float, ...] = (1.0, 2.0, 4.0),
    weights: float | tuple[float, ...] = (0.8, 0.35, 0.15),
    dtype=np.float64,
):
    """Sharpen a SunPy map by adding unsharp-mask details at multiple scales."""

    if not sigmas:
        raise ValueError("sigmas must contain at least one scale")
    if np.isscalar(weights):
        weights = tuple(float(weights) for _ in sigmas)
    else:
        weights = tuple(weights)
    if len(sigmas) != len(weights):
        raise ValueError("sigmas and weights must have the same length")

    data = smap.data.astype(dtype, copy=False)
    sharpened_data = data.copy()
    for sigma, weight in zip(sigmas, weights, strict=True):
        blurred = ndimage.gaussian_filter(data, sigma=sigma)
        sharpened_data += (data - blurred) * weight
    return Map(sharpened_data, smap.meta.copy())


def derotate_to_reference(
    path: str | Path,
    reference_map,
    aia_prep: bool = False,
    pointing_table=None,
    pointing_source: str = "JSOC",
    psf=None,
    deconvolve_psf: bool = False,
    deconvolve_iterations: int = 25,
    method: str = "pixel",
    order: int = 1,
    fill_value: float = np.nan,
):
    """Load a map and derotate it onto a reference map.

    This is the maintained version of the common `de_rotate` helper. For AIA
    level-1.5 preparation, install the `aia` extra and pass `aia_prep=True`.
    The default `method="pixel"` uses per-pixel differential-rotation remapping.
    Use `method="center"` for the faster center-point approximation.
    """

    if method not in {"center", "pixel"}:
        raise ValueError("method must be 'center' or 'pixel'")

    path = str(path)
    if "hmi" in path.lower():
        smap = Map(path).rotate()
        exposure_time = 1.0
    elif aia_prep:
        from aiapy.calibrate.utils import get_pointing_table

        smap = Map(path)
        if pointing_table is None:
            pointing_table = get_pointing_table(
                pointing_source,
                time_range=(smap.date - 12 * u.h, smap.date + 12 * u.h),
            )
        smap = globals()["aia_prep"](
            smap,
            pointing_table=pointing_table,
            psf=psf,
            deconvolve_psf=deconvolve_psf,
            deconvolve_iterations=deconvolve_iterations,
        )
        exposure_time = smap.exposure_time.value
    else:
        smap = Map(path)
        exposure_time = getattr(smap, "exposure_time", 1.0)
        exposure_time = getattr(exposure_time, "value", exposure_time)

    if method == "pixel":
        return derotate_to_reference_pixel(smap, reference_map, order=order, fill_value=fill_value)

    dy, dx = reference_map.data.shape
    new_center = solar_rotate_coordinate(reference_map.center, observer=smap.observer_coordinate)
    x_center, y_center = smap.wcs.world_to_pixel(new_center)
    x0 = int(x_center - (dx - 1) / 2)
    y0 = int(y_center - (dy - 1) / 2)
    shift_x = (x_center - (dx - 1) / 2) % 1
    shift_y = (y_center - (dy - 1) / 2) % 1

    reference_exposure_time = _map_exposure_value(reference_map, default=None)
    exposure_scale = 1.0 if reference_exposure_time is None else reference_exposure_time
    data = smap.data[y0 : y0 + dy + 1, x0 : x0 + dx + 1] / exposure_time * exposure_scale
    data = shift(data, [-shift_y, -shift_x], mode="nearest")
    data = data[:dy, :dx]
    meta = _derotated_reference_meta(
        reference_map,
        smap,
        exposure_time=reference_exposure_time,
        exposure_normalized=False,
    )
    return Map(data, meta)


def derotate_to_reference_pixel(
    source_map,
    reference_map,
    *,
    order: int = 1,
    fill_value: float = np.nan,
):
    """Derotate a source map onto every pixel of a reference map.

    This is more accurate than the center-point method for large submaps because
    each reference pixel is independently solar-rotated to the source observer
    time before sampling the source image.
    """

    yy, xx = np.indices(reference_map.data.shape, dtype=float)
    reference_coords = reference_map.wcs.pixel_to_world(xx, yy)
    source_coords = solar_rotate_coordinate(reference_coords, observer=source_map.observer_coordinate)
    source_x, source_y = source_map.wcs.world_to_pixel(source_coords)

    source_data = source_map.data.astype(float, copy=False)
    data = map_coordinates(
        source_data,
        [np.asarray(source_y, dtype=float), np.asarray(source_x, dtype=float)],
        order=order,
        mode="constant",
        cval=fill_value,
    )

    source_exposure_time = _map_exposure_value(source_map, default=None)
    reference_exposure_time = _map_exposure_value(reference_map, default=None)
    if source_exposure_time is not None and reference_exposure_time is not None:
        data = data / source_exposure_time * reference_exposure_time

    meta = _derotated_reference_meta(
        reference_map,
        source_map,
        exposure_time=reference_exposure_time,
        exposure_normalized=False,
    )
    return Map(data, meta)


def _derotated_reference_meta(
    reference_map,
    source_map,
    *,
    exposure_time: float = 1.0,
    exposure_normalized: bool = True,
):
    """Return reference-map metadata with source-map timing and exposure."""

    meta = reference_map.meta.copy()
    source_meta = source_map.meta
    meta.pop("blank", None)

    for key in ("date-obs", "date_obs", "t_obs", "date-beg", "date-avg", "date-end"):
        if key in source_meta:
            meta[key] = source_meta[key]
    if "date-obs" not in source_meta and "date_obs" not in source_meta:
        meta["date-obs"] = source_map.date.isot

    for key in ("exptime", "xposure", "texpmean"):
        if exposure_time is None:
            meta.pop(key, None)
        elif key in meta:
            meta[key] = exposure_time
    if exposure_time is not None:
        meta["exptime"] = exposure_time
    if exposure_normalized:
        meta["bunit"] = "DN / s"
    elif "bunit" in meta:
        pass
    elif "bunit" in source_meta:
        meta["bunit"] = source_meta["bunit"]
    else:
        meta["bunit"] = "DN"
    return meta


def _map_exposure_value(smap, default: float | None = 1.0) -> float | None:
    """Return a map exposure time as a plain float."""

    exposure_time = getattr(smap, "exposure_time", default)
    if exposure_time is None:
        return default
    return float(getattr(exposure_time, "value", exposure_time))


def _aiapy_imports():
    from aiapy import calibrate, psf

    return calibrate, psf


def _as_psf_array(psf):
    if isinstance(psf, (str, Path)):
        with np.load(psf) as data:
            return data["psf"]
    return psf

"""CHASE data readers."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from astropy.io import fits
import astropy.units as u


def read_chase(
    path: str | Path,
    hdu_index: int = 1,
    *,
    normalize_header: bool = True,
) -> tuple[np.ndarray, fits.Header, np.ndarray]:
    """Read a CHASE FITS spectral cube and its wavelength axis."""

    with fits.open(path) as hdul:
        data = hdul[hdu_index].data.astype(np.float32)
        header = hdul[hdu_index].header.copy()

    if normalize_header:
        header = normalize_chase_header(header)

    n_wavelength = int(header.get("NAXIS3", data.shape[0]))
    wavelength = header["CRVAL3"] + np.arange(n_wavelength) * header["CDELT3"]
    return data, header, wavelength


def normalize_chase_header(header: fits.Header) -> fits.Header:
    """Return a CHASE header with SunPy-friendly coordinate keywords."""

    header = header.copy()
    header.setdefault("CTYPE1", "HPLN-TAN")
    header.setdefault("CTYPE2", "HPLT-TAN")
    header.setdefault("CUNIT1", "arcsec")
    header.setdefault("CUNIT2", "arcsec")
    if "DATE-OBS" not in header and "DATE_OBS" in header:
        header["DATE-OBS"] = header["DATE_OBS"]
    if "DSUN_OBS" not in header:
        dsun_obs = _dsun_from_solar_radius(header)
        if dsun_obs is not None:
            header["DSUN_OBS"] = dsun_obs
    header.setdefault("HGLN_OBS", 0.0)
    header.setdefault("HGLT_OBS", 0.0)
    return header


def _dsun_from_solar_radius(header: fits.Header) -> float | None:
    rsun_ref = header.get("RSUN_REF") or header.get("R_SUN")
    rsun_obs = header.get("RSUN_OBS")
    if rsun_ref is None or rsun_obs is None:
        return None
    return float(rsun_ref) / np.tan(float(rsun_obs) * u.arcsec.to(u.rad))

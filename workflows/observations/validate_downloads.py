#!/usr/bin/env python3
"""Validate downloaded observation files by opening and summarizing them."""

from __future__ import annotations

import argparse
from pathlib import Path

from astropy.io import fits


FITS_SUFFIXES = {".fits", ".fit", ".fts"}
CDF_SUFFIXES = {".cdf"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("path", type=Path, help="File or directory to validate")
    parser.add_argument("--max-files", type=int, default=20)
    parser.add_argument("--sunpy-map", action="store_true", help="Also try opening FITS files as SunPy maps")
    args = parser.parse_args()

    files = collect_files(args.path, max_files=args.max_files)
    if not files:
        raise SystemExit(f"No FITS/CDF files found under {args.path}")

    for file_path in files:
        print(f"\n{file_path}")
        if file_path.suffix.lower() in FITS_SUFFIXES:
            summarize_fits(file_path, try_sunpy_map=args.sunpy_map)
        elif file_path.suffix.lower() in CDF_SUFFIXES:
            summarize_cdf(file_path)


def collect_files(path: Path, max_files: int) -> list[Path]:
    if path.is_file():
        return [path]

    files: list[Path] = []
    for suffix in sorted(FITS_SUFFIXES | CDF_SUFFIXES):
        files.extend(sorted(path.rglob(f"*{suffix}")))
    return files[:max_files]


def summarize_fits(path: Path, try_sunpy_map: bool = False) -> None:
    with fits.open(path) as hdul:
        print(f"  FITS HDUs: {len(hdul)}")
        for index, hdu in enumerate(hdul[:3]):
            shape = getattr(hdu.data, "shape", None)
            date_obs = hdu.header.get("DATE-OBS") or hdu.header.get("T_OBS")
            instrume = hdu.header.get("INSTRUME")
            detector = hdu.header.get("DETECTOR")
            wavelength = hdu.header.get("WAVELNTH")
            print(
                "  "
                f"HDU {index}: shape={shape}, instrume={instrume}, "
                f"detector={detector}, wavelength={wavelength}, date={date_obs}"
            )

    if try_sunpy_map:
        from sunpy.map import Map

        smap = Map(path)
        print(f"  SunPy map: shape={smap.data.shape}, date={smap.date}, observatory={smap.observatory}")


def summarize_cdf(path: Path) -> None:
    try:
        from cdflib import CDF
    except ImportError as exc:
        raise ImportError("CDF validation requires `cdflib`: pip install cdflib") from exc

    cdf = CDF(str(path))
    info = cdf.cdf_info()
    variables = list(info.zVariables or info.rVariables)
    print(f"  CDF variables: {len(variables)}")
    print(f"  First variables: {variables[:10]}")


if __name__ == "__main__":
    main()


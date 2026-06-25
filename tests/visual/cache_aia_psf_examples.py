#!/usr/bin/env python3
"""Cache AIA PSFs used by visual and preprocessing tests."""

from __future__ import annotations

import argparse
from pathlib import Path

from solarutils.obs.maps import AIA_EUV_WAVELENGTHS, DEFAULT_AIA_PSF_DIR, cache_aia_psfs


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_AIA_PSF_DIR)
    parser.add_argument("--wavelengths", nargs="+", type=int, default=list(AIA_EUV_WAVELENGTHS))
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    cached = cache_aia_psfs(args.wavelengths, args.output_dir, overwrite=args.overwrite)
    for wavelength, path in cached.items():
        print(f"{wavelength}: {path}")


if __name__ == "__main__":
    main()

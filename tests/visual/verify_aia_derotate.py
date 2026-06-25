#!/usr/bin/env python3
"""Derotate one AIA map to a reference submap and save a FITS smoke-test file."""

from __future__ import annotations

import argparse
from pathlib import Path

from sunpy.map import Map

from solarutils.obs.maps import derotate_to_reference, get_submap


DEFAULT_REFERENCE = Path(
    "data/test/observations/20241008_2300_2400/sdo/aia/"
    "aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits"
)
DEFAULT_INPUT = Path(
    "data/test/observations/20241008_2300_2400/sdo/aia/"
    "aia.lev1_euv_12s.2024-10-08T231006Z.193.image_lev1.fits"
)
DEFAULT_OUTPUT = Path("data/test/processed/observations/sdo/aia/193/derotate_20241008T231006_to_230006.fits")
DEFAULT_EXTENT = (100.0, 800.0, 150.0, 500.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--extent", nargs=4, type=float, default=DEFAULT_EXTENT)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    reference = get_submap(Map(args.reference), tuple(args.extent))
    derotated = derotate_to_reference(args.input, reference)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    derotated.save(args.output, overwrite=args.overwrite)

    print(args.output)
    print(f"reference shape: {reference.data.shape}")
    print(f"derotated shape: {derotated.data.shape}")
    print(f"reference date-obs: {reference.meta.get('date-obs')}")
    print(f"source date-obs: {Map(args.input).meta.get('date-obs')}")
    print(f"output date-obs: {derotated.meta.get('date-obs')}")
    print(f"output t_obs: {derotated.meta.get('t_obs')}")
    print(f"output exptime: {derotated.meta.get('exptime')}")
    print(f"output bunit: {derotated.meta.get('bunit')}")


if __name__ == "__main__":
    main()

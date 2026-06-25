#!/usr/bin/env python3
"""Calculate a time-distance array from processed AIA maps."""

from __future__ import annotations

import argparse
from glob import glob
from pathlib import Path

from solarutils.obs import time_distance_series


DEFAULT_INPUT = Path(
    "data/test/processed/observations/20241008_2300_2400/sdo/aia/0193/prep_psf_derotate_submap"
)
DEFAULT_OUTPUT = Path(
    "data/test/processed/observations/20241008_2300_2400/sdo/aia/0193/time_distance/aia193_time_distance.npz"
)
DEFAULT_START = (200.0, 260.0)
DEFAULT_END = (720.0, 430.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="Processed AIA FITS file, directory, or glob pattern",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--start", nargs=2, type=float, default=DEFAULT_START, metavar=("X", "Y"))
    parser.add_argument("--end", nargs=2, type=float, default=DEFAULT_END, metavar=("X", "Y"))
    parser.add_argument("--samples", type=int, default=200)
    parser.add_argument("--order", type=int, default=1)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    if args.output.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists, use --overwrite: {args.output}")

    paths = resolve_input_paths(args.input)
    result = time_distance_series(
        paths,
        tuple(args.start),
        tuple(args.end),
        args.samples,
        order=args.order,
        output=args.output,
    )

    print(args.output, flush=True)
    print(f"frames={len(result['time_iso'])}", flush=True)
    print(f"samples={len(result['distance_km'])}", flush=True)
    print(f"total_distance_km={result['total_distance_km']:.6g}", flush=True)


def resolve_input_paths(input_path: Path) -> list[Path]:
    input_text = str(input_path)
    if any(char in input_text for char in "*?[]"):
        paths = [Path(path) for path in sorted(glob(input_text))]
    elif input_path.is_dir():
        paths = sorted(input_path.glob("*.fits"))
    else:
        paths = [input_path]

    paths = [path for path in paths if path.exists()]
    if not paths:
        raise FileNotFoundError(f"No processed AIA FITS files found for {input_path}")
    return paths


if __name__ == "__main__":
    main()

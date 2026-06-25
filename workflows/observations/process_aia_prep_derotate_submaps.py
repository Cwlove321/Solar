#!/usr/bin/env python3
"""Prepare AIA maps with PSF deconvolution, derotate submaps, and save FITS files."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import astropy.units as u
from aiapy.calibrate.utils import get_pointing_table
from sunpy.map import Map

from solarutils.obs.maps import aia_prep, aia_psf_path, derotate_to_reference, get_submap


DEFAULT_INPUT_DIR = Path("data/test/observations/20241008_2300_2400/sdo/aia")
DEFAULT_OUTPUT_DIR = Path("data/test/processed/observations/20241008_2300_2400/sdo/aia")
DEFAULT_EXTENT = (100.0, 800.0, 150.0, 500.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--wavelengths", nargs="*", type=int, default=None)
    parser.add_argument("--extent", nargs=4, type=float, default=DEFAULT_EXTENT)
    parser.add_argument("--iterations", type=int, default=25)
    parser.add_argument("--pointing-source", default="JSOC")
    parser.add_argument(
        "--method",
        choices=("center", "pixel"),
        default="pixel",
        help="Use center-point derotation or per-pixel derotation",
    )
    parser.add_argument("--order", type=int, default=1, help="Interpolation order for --method pixel")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    files = sorted(args.input_dir.glob("*.fits"))
    if not files:
        raise FileNotFoundError(f"No FITS files found in {args.input_dir}")

    grouped = group_by_wavelength(files)
    if args.wavelengths is not None:
        grouped = {wave: grouped[wave] for wave in args.wavelengths if wave in grouped}
    if not grouped:
        raise ValueError("No AIA files matched the requested wavelengths")

    first_map = Map(files[0])
    pointing_table = get_pointing_table(
        args.pointing_source,
        time_range=(first_map.date - 12 * u.h, first_map.date + 12 * u.h),
    )

    written = []
    for wave, wave_files in sorted(grouped.items()):
        print(f"[{wave:04d}] {len(wave_files)} files", flush=True)
        psf = aia_psf_path(wave)
        reference_full = aia_prep(
            Map(wave_files[0]),
            pointing_table=pointing_table,
            psf=psf,
            deconvolve_psf=True,
            deconvolve_iterations=args.iterations,
        )
        reference_submap = get_submap(reference_full, tuple(args.extent))
        output_dir = args.output_dir / f"{wave:04d}" / output_dir_name(args.method)
        output_dir.mkdir(parents=True, exist_ok=True)

        for index, path in enumerate(wave_files):
            output = output_dir / output_name(path, wave_files[0], args.method)
            if output.exists() and not args.overwrite:
                print(f"  skip {output}", flush=True)
                continue

            if index == 0:
                prepared = clean_reference_output(reference_submap)
            else:
                prepared = derotate_to_reference(
                    path,
                    reference_submap,
                    aia_prep=True,
                    pointing_table=pointing_table,
                    psf=psf,
                    deconvolve_psf=True,
                    deconvolve_iterations=args.iterations,
                    method=args.method,
                    order=args.order,
                )
            prepared.save(output, overwrite=args.overwrite)
            written.append(output)
            print(f"  wrote {output}", flush=True)

    print(f"written={len(written)}", flush=True)


def group_by_wavelength(files: list[Path]) -> dict[int, list[Path]]:
    grouped = defaultdict(list)
    for path in files:
        wave = int(Map(path).wavelength.to_value("Angstrom"))
        grouped[wave].append(path)
    return dict(grouped)


def clean_reference_output(smap):
    meta = smap.meta.copy()
    meta.pop("blank", None)
    if "bunit" not in meta:
        meta["bunit"] = "DN"
    return Map(smap.data, meta)


def output_dir_name(method: str) -> str:
    if method == "center":
        return "prep_psf_derotate_submap"
    return f"prep_psf_derotate_submap_{method}"


def output_name(path: Path, reference_path: Path, method: str) -> str:
    date_token = _date_token(path)
    reference_token = _date_token(reference_path)
    method_suffix = "" if method == "center" else f"_{method}"
    return f"{date_token}_to_{reference_token}.prep_psf_derotate_submap{method_suffix}.fits"


def _date_token(path: Path) -> str:
    for part in path.name.split("."):
        if part.startswith("2024-"):
            return part.replace("-", "").replace(":", "")
    return path.stem


if __name__ == "__main__":
    main()

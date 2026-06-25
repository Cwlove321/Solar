#!/usr/bin/env python3
"""Command-line entry point for observational data search and download."""

from __future__ import annotations

import argparse
from pathlib import Path

from solarutils.obs.download import (
    carrington_rotation_for_time,
    fetch_results,
    search_hmi,
    search_hmi_mrdailysynframe_polfil,
    search_hmi_synoptic_mr_polfil,
    search_hmi_vector,
    search_sdo_aia,
    search_sharp,
    search_solar_orbiter,
    search_stereo_a,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    add_sdo_aia_parser(subparsers)
    add_hmi_parser(subparsers)
    add_hmi_vector_parser(subparsers)
    add_hmi_synoptic_parser(subparsers)
    add_hmi_synoptic_cr_parser(subparsers)
    add_sharp_parser(subparsers)
    add_stereo_parser(subparsers)
    add_solar_orbiter_parser(subparsers)

    args = parser.parse_args()
    plan = args.func(args)
    print(plan.description)
    print(plan.query_result)

    if args.fetch:
        files = fetch_results(plan, overwrite=args.overwrite)
        print(files)
    else:
        print("Search only. Add --fetch to download files.")


def add_common_time(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--start", required=True, help="Start time, e.g. 2024-10-08T22:00:00")
    parser.add_argument("--end", required=True, help="End time, e.g. 2024-10-08T23:00:00")
    parser.add_argument("--out", type=Path, default=None, help="Output directory")
    parser.add_argument("--sample-minutes", type=float, default=None, help="Sample cadence in minutes")
    parser.add_argument("--fetch", action="store_true", help="Download after search")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files when fetching")


def add_jsoc_email(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("--email", required=True, help="JSOC registered email address")


def add_sdo_aia_parser(subparsers) -> None:
    parser = subparsers.add_parser("sdo-aia", help="Search/download SDO AIA data from JSOC")
    add_common_time(parser)
    add_jsoc_email(parser)
    parser.add_argument("--wavelengths", nargs="+", type=int, default=[171, 193, 211, 304])
    parser.add_argument("--segments", nargs="+", default=["image"])
    parser.set_defaults(func=lambda args: search_sdo_aia(
        args.start,
        args.end,
        args.email,
        wavelengths=args.wavelengths,
        segments=args.segments,
        sample_minutes=args.sample_minutes,
        output_dir=args.out or "data/raw/sdo/aia",
    ))


def add_hmi_parser(subparsers) -> None:
    parser = subparsers.add_parser("hmi", help="Search/download HMI data from JSOC")
    add_common_time(parser)
    add_jsoc_email(parser)
    parser.add_argument("--series", default="hmi_m_45s", help="Alias or JSOC series, e.g. hmi_m_45s")
    parser.add_argument("--segments", nargs="+", default=["magnetogram"])
    parser.set_defaults(func=lambda args: search_hmi(
        args.start,
        args.end,
        args.email,
        series=args.series,
        segments=args.segments,
        sample_minutes=args.sample_minutes,
        output_dir=args.out or "data/raw/sdo/hmi",
    ))


def add_hmi_vector_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "hmi-vector",
        help="Search/download full-disk HMI 720s vector data from JSOC",
    )
    add_common_time(parser)
    add_jsoc_email(parser)
    parser.add_argument(
        "--segments",
        nargs="+",
        default=["field", "inclination", "azimuth", "disambig"],
        help="HMI vector segments to download",
    )
    parser.set_defaults(func=lambda args: search_hmi_vector(
        args.start,
        args.end,
        args.email,
        segments=args.segments,
        sample_minutes=args.sample_minutes,
        output_dir=args.out or "data/raw/sdo/hmi_b_720s",
    ))


def add_hmi_synoptic_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "hmi-synoptic",
        help="Search/download HMI mrdailysynframe_polfil 720s data from JSOC",
    )
    add_common_time(parser)
    add_jsoc_email(parser)
    parser.add_argument(
        "--segments",
        nargs="+",
        default=["Mr_polfil"],
        help="Synoptic-frame segments to download",
    )
    parser.set_defaults(func=lambda args: search_hmi_mrdailysynframe_polfil(
        args.start,
        args.end,
        args.email,
        segments=args.segments,
        sample_minutes=args.sample_minutes,
        output_dir=args.out or "data/raw/sdo/hmi_mrdailysynframe_polfil_720s",
    ))


def add_hmi_synoptic_cr_parser(subparsers) -> None:
    parser = subparsers.add_parser(
        "hmi-synoptic-cr",
        help="Search/download HMI Carrington-rotation synoptic Mr_polfil data from JSOC",
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--car-rotation", type=int, help="Carrington rotation number, e.g. 2289")
    group.add_argument("--time", help="Observation time used to infer the Carrington rotation")
    parser.add_argument("--out", type=Path, default=None, help="Output directory")
    parser.add_argument("--segments", nargs="+", default=["Mr_polfil"])
    parser.add_argument("--fetch", action="store_true", help="Download after search")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files when fetching")
    add_jsoc_email(parser)
    parser.set_defaults(func=lambda args: search_hmi_synoptic_mr_polfil(
        args.car_rotation if args.car_rotation is not None else carrington_rotation_for_time(args.time),
        args.email,
        segments=args.segments,
        output_dir=args.out or "data/raw/sdo/hmi_synoptic_mr_polfil_720s",
    ))


def add_sharp_parser(subparsers) -> None:
    parser = subparsers.add_parser("sharp", help="Search/download SHARP CEA data from JSOC")
    add_common_time(parser)
    add_jsoc_email(parser)
    parser.add_argument("--harpnum", type=int, default=None)
    parser.add_argument("--series", default="hmi.sharp_cea_720s")
    parser.add_argument("--segments", nargs="+", default=["Br", "Bp", "Bt", "bitmap", "conf_disambig"])
    parser.set_defaults(func=lambda args: search_sharp(
        args.start,
        args.end,
        args.email,
        harpnum=args.harpnum,
        series=args.series,
        segments=args.segments,
        sample_minutes=args.sample_minutes,
        output_dir=args.out or "data/raw/sdo/sharp",
    ))


def add_stereo_parser(subparsers) -> None:
    parser = subparsers.add_parser("stereo-a", help="Search/download STEREO-A data through VSO")
    add_common_time(parser)
    parser.add_argument("--instrument", default="SECCHI")
    parser.add_argument("--detector", default="EUVI")
    parser.add_argument("--wavelength", type=int, default=None)
    parser.set_defaults(func=lambda args: search_stereo_a(
        args.start,
        args.end,
        instrument=args.instrument,
        detector=args.detector,
        wavelength=args.wavelength,
        sample_minutes=args.sample_minutes,
        output_dir=args.out or "data/raw/stereo_a",
    ))


def add_solar_orbiter_parser(subparsers) -> None:
    parser = subparsers.add_parser("solar-orbiter", help="Search/download Solar Orbiter data through SOAR")
    add_common_time(parser)
    parser.add_argument("--instrument", default="EUI")
    parser.add_argument("--product", default=None)
    parser.add_argument("--level", default="L1")
    parser.set_defaults(func=lambda args: search_solar_orbiter(
        args.start,
        args.end,
        instrument=args.instrument,
        product=args.product,
        level=args.level,
        sample_minutes=args.sample_minutes,
        output_dir=args.out or "data/raw/solar_orbiter",
    ))


if __name__ == "__main__":
    main()

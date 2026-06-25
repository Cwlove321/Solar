#!/usr/bin/env python3
"""Derotate HMI 45-second LOS magnetogram submaps and save FITS files."""

from __future__ import annotations

import argparse
from pathlib import Path

from sunpy.map import Map
from tqdm import tqdm

from solarutils.obs.maps import derotate_to_reference, get_submap


DEFAULT_INPUT_DIR = Path("data/test/observations/20241008_2300_2400/sdo/hmi_45s")
DEFAULT_EXTENT = (100.0, 800.0, 150.0, 500.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to a sibling `sub_m45` directory next to --input-dir.",
    )
    parser.add_argument("--extent", nargs=4, type=float, default=DEFAULT_EXTENT)
    parser.add_argument(
        "--reference-time",
        default=None,
        help="Reference time token, for example 20241008_230000. Defaults to the first FITS file.",
    )
    parser.add_argument(
        "--method",
        choices=("center", "pixel"),
        default="pixel",
        help="Use center-point derotation or per-pixel derotation",
    )
    parser.add_argument("--order", type=int, default=1, help="Interpolation order for --method pixel")
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    files = sorted(path for path in args.input_dir.glob("*.fits") if not path.name.startswith("._"))
    if not files:
        raise FileNotFoundError(f"No HMI 45s FITS files found in {args.input_dir}")

    reference_path = select_reference_file(files, args.reference_time)
    reference_full = Map(reference_path).rotate()
    reference_submap = get_submap(reference_full, tuple(args.extent))
    reference_submap = clean_hmi_output(reference_submap)
    output_dir = args.output_dir or args.input_dir.parent / "sub_m45"
    output_dir.mkdir(parents=True, exist_ok=True)

    written = []
    skipped = []
    for path in tqdm(files, desc="Derotate HMI45", unit="frame"):
        output = output_dir / output_name(path, reference_path, args.method)
        if output.exists() and not args.overwrite:
            skipped.append(output)
            continue

        if path == reference_path:
            prepared = reference_submap
        else:
            prepared = derotate_to_reference(
                path,
                reference_submap,
                method=args.method,
                order=args.order,
            )
            prepared = clean_hmi_output(prepared)
        prepared.save(output, overwrite=args.overwrite)
        written.append(output)

    print(f"written={len(written)}", flush=True)
    print(f"skipped={len(skipped)}", flush=True)
    print(f"reference={reference_path}", flush=True)
    print(f"output_dir={output_dir}", flush=True)


def clean_hmi_output(smap):
    meta = smap.meta.copy()
    meta.pop("blank", None)
    meta.pop("exptime", None)
    meta.pop("xposure", None)
    meta.pop("texpmean", None)
    if not meta.get("bunit"):
        meta["bunit"] = "Gauss"
    return Map(smap.data, meta)


def select_reference_file(files: list[Path], reference_time: str | None = None) -> Path:
    if reference_time is None:
        return files[0]
    normalized = reference_time.replace("_TAI", "")
    for path in files:
        if _date_token(path) == normalized:
            return path
    raise FileNotFoundError(f"No reference FITS matched --reference-time {reference_time}")


def output_name(path: Path, reference_path: Path, method: str) -> str:
    date_token = _date_token(path)
    reference_token = _date_token(reference_path)
    return f"{date_token}_to_{reference_token}.hmi45_los_derotate_submap_{method}.fits"


def _date_token(path: Path) -> str:
    parts = path.name.split(".")
    for part in parts:
        if "_TAI" in part:
            return part.replace("_TAI", "")
    return path.stem


if __name__ == "__main__":
    main()

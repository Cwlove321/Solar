#!/usr/bin/env python
"""Convert a chronological HMI B_720s sequence to flat SHARP-like CEA Br/Bt/Bp files."""

from __future__ import annotations

import argparse
from pathlib import Path
import re
import shutil
import subprocess
import sys


SEGMENTS = ("field", "inclination", "azimuth", "disambig")
COMPONENTS = ("Br", "Bt", "Bp")
NAME_PATTERN = re.compile(
    r"hmi\.b_720s\.(?P<time>\d{8}_\d{6}_TAI)\.(?P<segment>field|inclination|azimuth|disambig)\.fits$"
)


def main() -> None:
    args = parse_args()
    records = discover_records(args.input_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    metadata_dir = args.output_dir / "metadata"
    work_root = args.output_dir / ".work"
    metadata_dir.mkdir(exist_ok=True)
    work_root.mkdir(exist_ok=True)

    script = Path(__file__).resolve().with_name("make_cea_patch.py")
    for index, (timestamp, segments) in enumerate(records.items(), start=1):
        prefix = f"hmi.sharp_cea_720s.{args.harpnum}.{timestamp}"
        outputs = {component: args.output_dir / f"{prefix}.{component}.fits" for component in COMPONENTS}
        if all(path.exists() and path.stat().st_size > 0 for path in outputs.values()):
            print(f"[{index}/{len(records)}] skip complete {timestamp}", flush=True)
            continue

        work_dir = work_root / timestamp
        if work_dir.exists():
            shutil.rmtree(work_dir)
        work_dir.mkdir(parents=True)
        command = [
            sys.executable,
            str(script),
            "--field", str(segments["field"]),
            "--inclination", str(segments["inclination"]),
            "--azimuth", str(segments["azimuth"]),
            "--disambig", str(segments["disambig"]),
            "--output-dir", str(work_dir),
            "--center-lon", str(args.center_lon),
            "--center-lat", str(args.center_lat),
            "--width-degree", str(args.width_degree),
            "--height-degree", str(args.height_degree),
            "--resolution-degree", str(args.resolution_degree),
            "--sampling-mode", args.sampling_mode,
            "--oversample-resolution-degree", str(args.oversample_resolution_degree),
            "--smooth-sigma-degree", str(args.smooth_sigma_degree),
            "--smooth-truncate", str(args.smooth_truncate),
            "--components-only",
            "--no-quicklook",
        ]
        print(f"[{index}/{len(records)}] convert {timestamp}", flush=True)
        subprocess.run(command, check=True)
        for component, destination in outputs.items():
            source = work_dir / f"{component}.fits"
            if not source.exists() or source.stat().st_size == 0:
                raise RuntimeError(f"missing generated component: {source}")
            source.replace(destination)
        (work_dir / "parameters.json").replace(metadata_dir / f"{timestamp}.json")
        shutil.rmtree(work_dir)

    if work_root.exists() and not any(work_root.iterdir()):
        work_root.rmdir()
    print(f"Completed {len(records)} timestamps in {args.output_dir}", flush=True)


def discover_records(input_dir: Path) -> dict[str, dict[str, Path]]:
    records: dict[str, dict[str, Path]] = {}
    for path in sorted(input_dir.glob("hmi.b_720s.*.fits")):
        if path.name.startswith("._"):
            continue
        match = NAME_PATTERN.fullmatch(path.name)
        if match:
            records.setdefault(match.group("time"), {})[match.group("segment")] = path
    if not records:
        raise FileNotFoundError(f"no HMI B_720s files found in {input_dir}")
    incomplete = {time: sorted(set(SEGMENTS) - set(parts)) for time, parts in records.items() if set(parts) != set(SEGMENTS)}
    if incomplete:
        raise ValueError(f"incomplete timestamps: {incomplete}")
    return dict(sorted(records.items()))


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--harpnum", required=True, type=int)
    parser.add_argument("--center-lon", required=True, type=float)
    parser.add_argument("--center-lat", required=True, type=float)
    parser.add_argument("--width-degree", required=True, type=float)
    parser.add_argument("--height-degree", required=True, type=float)
    parser.add_argument("--resolution-degree", type=float, default=0.03)
    parser.add_argument("--sampling-mode", choices=("direct", "sharp"), default="sharp")
    parser.add_argument("--oversample-resolution-degree", type=float, default=0.01)
    parser.add_argument("--smooth-sigma-degree", type=float, default=0.01)
    parser.add_argument("--smooth-truncate", type=float, default=2.0)
    return parser.parse_args()


if __name__ == "__main__":
    main()

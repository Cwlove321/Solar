#!/usr/bin/env python3
"""Create a differentially tracked SHARP-like CEA series from HMI B_720s files.

The user supplies the Stonyhurst longitude and latitude of the patch center at
one reference record.  For every other record, SunPy differentially rotates
that center to the target observation time before ``make_cea_patch.py`` remaps
the four HMI vector segments.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys

import astropy.units as u
from astropy.coordinates import SkyCoord
import numpy as np
from sunpy.coordinates import frames
from sunpy.map import Map
from sunpy.physics.differential_rotation import solar_rotate_coordinate


SEGMENTS = ("field", "inclination", "azimuth", "disambig")
COMPONENTS = ("Br", "Bt", "Bp")
ROTATION_MODELS = ("howard", "snodgrass", "allen", "rigid")
NAME_PATTERN = re.compile(
    r"hmi\.[bB]_720s\.(?P<time>\d{8}_\d{6}_TAI)\."
    r"(?P<segment>field|inclination|azimuth|disambig)\.fits(?:\.gz)?$",
    re.IGNORECASE,
)
TRACKING_SCHEMA_VERSION = 2


def map_reference_time(smap):
    """Return the WCS reference time (HMI T_OBS), not exposure start DATE-OBS."""

    return getattr(smap, "reference_date", smap.date)


def discover_records(input_dir: Path) -> dict[str, dict[str, Path]]:
    """Return complete B_720s records ordered by the timestamp in the filename."""

    if not input_dir.is_dir():
        raise FileNotFoundError(f"input directory does not exist: {input_dir}")

    records: dict[str, dict[str, Path]] = {}
    for path in sorted(input_dir.iterdir()):
        if not path.is_file() or path.name.startswith("._"):
            continue
        match = NAME_PATTERN.fullmatch(path.name)
        if match is None:
            continue
        timestamp = match.group("time").upper()
        segment = match.group("segment").lower()
        previous = records.setdefault(timestamp, {}).get(segment)
        if previous is not None:
            raise ValueError(
                f"duplicate {segment!r} segment for {timestamp}: {previous} and {path}"
            )
        records[timestamp][segment] = path

    if not records:
        raise FileNotFoundError(f"no HMI B_720s vector-segment files found in {input_dir}")

    incomplete = {
        timestamp: sorted(set(SEGMENTS) - set(parts))
        for timestamp, parts in records.items()
        if set(parts) != set(SEGMENTS)
    }
    if incomplete:
        raise ValueError(f"incomplete timestamps: {incomplete}")
    return dict(sorted(records.items()))


def select_reference_timestamp(records: dict[str, dict[str, Path]], requested: str | None) -> str:
    """Resolve the reference record, defaulting to the first input timestamp."""

    if requested is None:
        return next(iter(records))
    requested = requested.upper()
    if requested not in records:
        available = ", ".join(records)
        raise ValueError(
            f"reference timestamp {requested!r} is not a complete input record; "
            f"available timestamps: {available}"
        )
    return requested


def reference_coordinate(reference_map, center_lon: float, center_lat: float) -> SkyCoord:
    """Build the physical surface coordinate that anchors the tracked patch."""

    return SkyCoord(
        center_lon * u.deg,
        center_lat * u.deg,
        reference_map.rsun_meters,
        frame=frames.HeliographicStonyhurst,
        obstime=map_reference_time(reference_map),
    )


def rotate_center(reference_center: SkyCoord, target_map, model: str = "howard") -> tuple[float, float]:
    """Rotate a reference center to a target map and return Stonyhurst lon/lat."""

    rotated = solar_rotate_coordinate(
        reference_center,
        observer=target_map.observer_coordinate,
        model=model,
    )
    target_hgs = rotated.transform_to(
        frames.HeliographicStonyhurst(obstime=map_reference_time(target_map))
    )
    longitude = float(target_hgs.lon.wrap_at(180 * u.deg).to_value(u.deg))
    latitude = float(target_hgs.lat.to_value(u.deg))
    if not np.isfinite([longitude, latitude]).all():
        raise ValueError(f"non-finite tracked center at {map_reference_time(target_map).isot}")
    return longitude, latitude


def build_tracking_plan(
    records: dict[str, dict[str, Path]],
    reference_timestamp: str,
    center_lon: float,
    center_lat: float,
    model: str,
) -> tuple[dict[str, object], list[dict[str, object]]]:
    """Read map geometry and calculate the tracked center at every timestamp."""

    reference_map = Map(records[reference_timestamp]["field"])
    reference_time = map_reference_time(reference_map)
    center = reference_coordinate(reference_map, center_lon, center_lat)
    plan: list[dict[str, object]] = []
    for timestamp, segments in records.items():
        target_map = Map(segments["field"])
        target_time = map_reference_time(target_map)
        tracked_lon, tracked_lat = rotate_center(center, target_map, model=model)
        plan.append(
            {
                "timestamp": timestamp,
                "observation_time": target_time.utc.isot,
                "observation_time_scale": "UTC",
                "delta_days": float((target_time - reference_time).to_value(u.day)),
                "center_lon": tracked_lon,
                "center_lat": tracked_lat,
                "segments": {name: str(path) for name, path in segments.items()},
            }
        )

    reference = {
        "timestamp": reference_timestamp,
        "observation_time": reference_time.utc.isot,
        "observation_time_scale": "UTC",
        "center_lon": float(center_lon),
        "center_lat": float(center_lat),
        "longitude_frame": "HeliographicStonyhurst",
    }
    return reference, plan


def write_center_table(path: Path, plan: list[dict[str, object]]) -> None:
    """Write a compact, human-readable record of the center trajectory."""

    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=(
                "timestamp",
                "observation_time",
                "observation_time_scale",
                "delta_days",
                "center_lon",
                "center_lat",
            ),
        )
        writer.writeheader()
        for item in plan:
            writer.writerow({key: item[key] for key in writer.fieldnames})


def tracking_metadata(
    *,
    reference: dict[str, object],
    record: dict[str, object],
    rotation_model: str,
    processing: dict[str, object],
) -> dict[str, object]:
    return {
        "schema_version": TRACKING_SCHEMA_VERSION,
        "method": "sunpy.physics.differential_rotation.solar_rotate_coordinate",
        "rotation_model": rotation_model,
        "reference": reference,
        "record": {
            key: record[key]
            for key in (
                "timestamp",
                "observation_time",
                "observation_time_scale",
                "delta_days",
                "center_lon",
                "center_lat",
            )
        },
        "processing": processing,
        "scope": "center-point tracking; the CEA patch is not warped pixel by pixel",
    }


def processing_settings(args: argparse.Namespace) -> dict[str, object]:
    """Return output-affecting options used for safe resume checks."""

    fill_value = float(args.fill_value)
    return {
        "grid": {
            "width_degree": float(args.width_degree),
            "height_degree": float(args.height_degree),
            "resolution_degree": float(args.resolution_degree),
        },
        "sampling": {
            "backend": "jsoc-wiener6-oversample3-gaussian1-v1",
            "mode": args.sampling_mode,
            "interpolation_order": int(args.interpolation_order),
            "oversample_resolution_degree": float(args.oversample_resolution_degree),
            "smooth_sigma_degree": float(args.smooth_sigma_degree),
            "smooth_truncate": float(args.smooth_truncate),
            "fill_value": fill_value if np.isfinite(fill_value) else None,
            "disambig_method": int(args.disambig_method),
        },
    }


def outputs_are_current(
    outputs: dict[str, Path],
    metadata_path: Path,
    expected_tracking: dict[str, object],
) -> bool:
    """Return whether complete outputs were made with this tracking trajectory."""

    if not all(path.is_file() and path.stat().st_size > 0 for path in outputs.values()):
        return False
    try:
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return False
    actual = metadata.get("tracking")
    if not isinstance(actual, dict):
        return False
    return tracking_values_match(actual, expected_tracking)


def tracking_values_match(actual: dict[str, object], expected: dict[str, object]) -> bool:
    """Compare tracking metadata, allowing harmless floating-point roundoff."""

    return metadata_values_match(actual, expected)


def metadata_values_match(actual: object, expected: object) -> bool:
    """Recursively compare JSON-like values with a float tolerance."""

    if isinstance(expected, dict):
        if not isinstance(actual, dict) or set(actual) != set(expected):
            return False
        return all(metadata_values_match(actual[key], value) for key, value in expected.items())
    if isinstance(expected, list):
        if not isinstance(actual, list) or len(actual) != len(expected):
            return False
        return all(metadata_values_match(a, e) for a, e in zip(actual, expected, strict=True))
    if isinstance(expected, float):
        try:
            return bool(np.isclose(float(actual), expected, rtol=0.0, atol=1e-10))
        except (TypeError, ValueError):
            return False
    return actual == expected


def conversion_command(
    args: argparse.Namespace,
    script: Path,
    segments: dict[str, Path],
    work_dir: Path,
    center_lon: float,
    center_lat: float,
) -> list[str]:
    command = [
        sys.executable,
        str(script),
        "--field", str(segments["field"]),
        "--inclination", str(segments["inclination"]),
        "--azimuth", str(segments["azimuth"]),
        "--disambig", str(segments["disambig"]),
        "--output-dir", str(work_dir),
        "--center-lon", f"{center_lon:.12f}",
        "--center-lat", f"{center_lat:.12f}",
        "--width-degree", str(args.width_degree),
        "--height-degree", str(args.height_degree),
        "--resolution-degree", str(args.resolution_degree),
        "--disambig-method", str(args.disambig_method),
        "--interpolation-order", str(args.interpolation_order),
        "--sampling-mode", args.sampling_mode,
        "--oversample-resolution-degree", str(args.oversample_resolution_degree),
        "--smooth-sigma-degree", str(args.smooth_sigma_degree),
        "--smooth-truncate", str(args.smooth_truncate),
        "--fill-value", str(args.fill_value),
        "--components-only",
    ]
    command.append("--quicklook" if args.quicklook else "--no-quicklook")
    return command


def write_json(path: Path, value: dict[str, object]) -> None:
    path.write_text(json.dumps(value, indent=2, allow_nan=False), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path, help="Directory containing HMI B_720s segments.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Directory for tracked CEA component files.")
    parser.add_argument("--harpnum", required=True, type=int, help="Integer used in SHARP-like output filenames.")
    parser.add_argument(
        "--reference-time",
        help="Reference filename timestamp (YYYYMMDD_HHMMSS_TAI); default is the first complete record.",
    )
    parser.add_argument("--center-lon", required=True, type=float, help="Reference-center Stonyhurst longitude in deg.")
    parser.add_argument("--center-lat", required=True, type=float, help="Reference-center Stonyhurst latitude in deg.")
    parser.add_argument("--width-degree", required=True, type=float, help="CEA patch width in degrees.")
    parser.add_argument("--height-degree", required=True, type=float, help="CEA patch height in degrees.")
    parser.add_argument("--resolution-degree", type=float, default=0.03, help="CEA pixel size in degrees.")
    parser.add_argument(
        "--rotation-model",
        choices=ROTATION_MODELS,
        default="howard",
        help="SunPy differential-rotation profile; default: howard.",
    )
    parser.add_argument("--disambig-method", type=int, default=2, help="HMI disambiguation bit method.")
    parser.add_argument("--interpolation-order", type=int, default=3, choices=range(0, 6), help="Spline order.")
    parser.add_argument("--sampling-mode", choices=("sharp",), default="sharp")
    parser.add_argument("--oversample-resolution-degree", type=float)
    parser.add_argument("--smooth-sigma-degree", type=float)
    parser.add_argument("--smooth-truncate", type=float, default=2.0)
    parser.add_argument("--fill-value", type=float, default=np.nan)
    parser.add_argument("--quicklook", action="store_true", help="Write one quicklook PNG per timestamp.")
    parser.add_argument("--overwrite", action="store_true", help="Rebuild outputs even if their tracking metadata matches.")
    parser.add_argument(
        "--centers-only",
        action="store_true",
        help="Only calculate tracked_centers.csv and tracking_manifest.json; do not remap FITS files.",
    )
    args = parser.parse_args(argv)
    for name in ("oversample_resolution_degree", "smooth_sigma_degree"):
        if getattr(args, name) is None:
            setattr(args, name, args.resolution_degree / 3)
    if args.width_degree <= 0 or args.height_degree <= 0 or args.resolution_degree <= 0:
        parser.error("patch dimensions and --resolution-degree must be positive")
    if not -90.0 <= args.center_lat <= 90.0:
        parser.error("--center-lat must be between -90 and 90 degrees")
    return args


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    records = discover_records(args.input_dir)
    reference_timestamp = select_reference_timestamp(records, args.reference_time)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Calculating {len(records)} tracked centers using {args.rotation_model!r} rotation...", flush=True)
    reference, plan = build_tracking_plan(
        records,
        reference_timestamp,
        args.center_lon,
        args.center_lat,
        args.rotation_model,
    )
    center_table = args.output_dir / "tracked_centers.csv"
    manifest_path = args.output_dir / "tracking_manifest.json"
    processing = processing_settings(args)
    write_center_table(center_table, plan)

    manifest: dict[str, object] = {
        "schema_version": TRACKING_SCHEMA_VERSION,
        "input_dir": str(args.input_dir),
        "output_dir": str(args.output_dir),
        "reference": reference,
        "tracking": {
            "method": "sunpy.physics.differential_rotation.solar_rotate_coordinate",
            "rotation_model": args.rotation_model,
            "scope": "center-point tracking; the CEA patch is not warped pixel by pixel",
        },
        "processing": processing,
        "center_table": str(center_table),
        "records": plan,
    }
    write_json(manifest_path, manifest)

    if args.centers_only:
        for item in plan:
            print(
                f"{item['timestamp']}  lon={item['center_lon']:.6f} deg  "
                f"lat={item['center_lat']:.6f} deg",
                flush=True,
            )
        print(f"Wrote {center_table} and {manifest_path}", flush=True)
        return 0

    metadata_dir = args.output_dir / "metadata"
    work_root = args.output_dir / ".work"
    metadata_dir.mkdir(exist_ok=True)
    work_root.mkdir(exist_ok=True)
    script = Path(__file__).resolve().with_name("make_cea_patch.py")

    for index, record in enumerate(plan, start=1):
        timestamp = str(record["timestamp"])
        segments = records[timestamp]
        prefix = f"hmi.sharp_cea_720s.{args.harpnum}.{timestamp}"
        outputs = {
            component: args.output_dir / f"{prefix}.{component}.fits"
            for component in COMPONENTS
        }
        output_metadata = metadata_dir / f"{timestamp}.json"
        expected_tracking = tracking_metadata(
            reference=reference,
            record=record,
            rotation_model=args.rotation_model,
            processing=processing,
        )
        if not args.overwrite and outputs_are_current(outputs, output_metadata, expected_tracking):
            record["status"] = "skipped_current"
            print(f"[{index}/{len(plan)}] skip current {timestamp}", flush=True)
            write_json(manifest_path, manifest)
            continue

        work_dir = work_root / timestamp
        if work_dir.exists():
            shutil.rmtree(work_dir)
        work_dir.mkdir(parents=True)
        command = conversion_command(
            args,
            script,
            segments,
            work_dir,
            float(record["center_lon"]),
            float(record["center_lat"]),
        )
        print(
            f"[{index}/{len(plan)}] convert {timestamp}: "
            f"center=({record['center_lon']:.6f}, {record['center_lat']:.6f}) deg",
            flush=True,
        )
        subprocess.run(command, check=True)

        for component, destination in outputs.items():
            source = work_dir / f"{component}.fits"
            if not source.is_file() or source.stat().st_size == 0:
                raise RuntimeError(f"missing generated component: {source}")

        generated_metadata_path = work_dir / "parameters.json"
        generated_metadata = json.loads(generated_metadata_path.read_text(encoding="utf-8"))
        generated_metadata["tracking"] = expected_tracking
        write_json(generated_metadata_path, generated_metadata)
        for component, destination in outputs.items():
            (work_dir / f"{component}.fits").replace(destination)
        generated_metadata_path.replace(output_metadata)
        if args.quicklook:
            (work_dir / "quicklook.png").replace(args.output_dir / f"{prefix}.quicklook.png")
        shutil.rmtree(work_dir)
        record["status"] = "converted"
        record["outputs"] = {name: str(path) for name, path in outputs.items()}
        record["metadata"] = str(output_metadata)
        write_json(manifest_path, manifest)

    if work_root.exists() and not any(work_root.iterdir()):
        work_root.rmdir()
    print(f"Completed {len(plan)} timestamps in {args.output_dir}", flush=True)
    print(f"Center trajectory: {center_table}", flush=True)
    print(f"Tracking manifest: {manifest_path}", flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

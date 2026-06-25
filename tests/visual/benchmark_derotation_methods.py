#!/usr/bin/env python3
"""Benchmark center-point and per-pixel derotation on local AIA/HMI test data."""

from __future__ import annotations

import argparse
from pathlib import Path
import time

from sunpy.map import Map

from solarutils.obs.maps import derotate_to_reference, get_submap


DEFAULT_AIA_DIR = Path("data/test/observations/20241008_2300_2400/sdo/aia")
DEFAULT_HMI_DIR = Path("data/test/observations/20241008_2300_2400/sdo/hmi_45s")
DEFAULT_OUTPUT = Path("data/test/processed/benchmarks/derotation_methods.txt")
DEFAULT_EXTENT = (100.0, 800.0, 150.0, 500.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aia-dir", type=Path, default=DEFAULT_AIA_DIR)
    parser.add_argument("--hmi-dir", type=Path, default=DEFAULT_HMI_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--extent", nargs=4, type=float, default=DEFAULT_EXTENT)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--aia-wavelength", type=int, default=193)
    parser.add_argument("--order", type=int, default=1)
    args = parser.parse_args()

    rows = []
    hmi_files = sorted(args.hmi_dir.glob("*.fits"))
    if len(hmi_files) >= 2:
        reference = get_submap(Map(hmi_files[0]).rotate(), tuple(args.extent))
        rows.extend(benchmark_pair("hmi45", hmi_files[-1], reference, args.repeats, args.order))

    aia_files = sorted(args.aia_dir.glob(f"*{args.aia_wavelength}*.fits"))
    if len(aia_files) >= 2:
        reference = get_submap(Map(aia_files[0]), tuple(args.extent))
        rows.extend(benchmark_pair(f"aia{args.aia_wavelength}", aia_files[-1], reference, args.repeats, args.order))

    if not rows:
        raise FileNotFoundError("Need at least two AIA or HMI FITS files for benchmarking")

    lines = [
        "# Derotation method benchmark",
        f"extent_arcsec={tuple(args.extent)}",
        f"repeats={args.repeats}",
        f"pixel_order={args.order}",
        "",
        "dataset method seconds_best seconds_mean",
    ]
    for row in rows:
        lines.append(
            f"{row['dataset']} {row['method']} {row['seconds_best']:.6f} {row['seconds_mean']:.6f}"
        )

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print("\n".join(lines))
    print(args.output)


def benchmark_pair(dataset: str, source_path: Path, reference_map, repeats: int, order: int) -> list[dict[str, float | str]]:
    rows = []
    for method in ("center", "pixel"):
        seconds = []
        for _ in range(max(repeats, 1)):
            start = time.perf_counter()
            derotate_to_reference(source_path, reference_map, method=method, order=order)
            seconds.append(time.perf_counter() - start)
        rows.append(
            {
                "dataset": dataset,
                "method": method,
                "seconds_best": min(seconds),
                "seconds_mean": sum(seconds) / len(seconds),
            }
        )
    return rows


if __name__ == "__main__":
    main()

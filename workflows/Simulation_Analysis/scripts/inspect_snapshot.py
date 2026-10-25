#!/usr/bin/env python3
"""Inspect AMRVAC DAT metadata, ranges, and one physical extremum."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import data_paths, load_config
from simulation_analysis.snapshot import inspect_snapshot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dat", type=Path)
    parser.add_argument("--variable", default="temperature_K")
    parser.add_argument("--extreme", choices=("min", "max", "absmax"), default="max")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    paths = data_paths(config)
    source = args.dat or (paths[-1] if paths else None)
    if source is None:
        raise SystemExit("No DAT selected and data_glob matched nothing")
    report = inspect_snapshot(source, config, variable=args.variable, extreme=args.extreme)
    output = args.output or Path(config["output_dir"]) / "snapshot_inspection.json"
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2), flush=True)
    print(output.resolve(), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Plot derived fields on one fixed-r, fixed-theta, or fixed-phi slice."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import data_paths, load_config
from simulation_analysis.plotting import plot_plane_panels
from simulation_analysis.snapshot import extract_axis_slice, to_plot_plane


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dat", type=Path)
    parser.add_argument("--axis", choices=("r", "theta", "phi"))
    parser.add_argument("--coordinate", type=float)
    parser.add_argument("--fields", nargs="+", default=["temperature_K", "number_density_cm3", "amr_level"])
    parser.add_argument("--fieldlines", action="store_true")
    parser.add_argument("--save-npz", action="store_true")
    parser.add_argument("--columns", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    paths = data_paths(config)
    source = args.dat or (paths[-1] if paths else None)
    if source is None:
        raise SystemExit("No DAT selected and data_glob matched nothing")
    section = extract_axis_slice(
        source,
        config,
        axis=args.axis,
        coordinate=args.coordinate,
    )
    plane = to_plot_plane(section, config)
    output = args.output or Path(config["output_dir"]) / (
        f"slice_{section.fixed_axis}_{section.coordinate:.3f}_{Path(source).stem}.png"
    )
    result = plot_plane_panels(
        plane,
        args.fields,
        config,
        output,
        fieldlines=args.fieldlines,
        columns=args.columns,
    )
    metrics = {
        "source": str(Path(source).resolve()),
        "output": str(result),
        "axis": section.fixed_axis,
        "requested_coordinate": section.requested_coordinate,
        "actual_coordinate": section.coordinate,
        "time": section.time,
        "target_level": section.target_level,
        "fields": {
            field: {
                "min": float(np.nanmin(section.level if field == "amr_level" else plane.values[field])),
                "max": float(np.nanmax(section.level if field == "amr_level" else plane.values[field])),
            }
            for field in args.fields
        },
    }
    output.with_suffix(".json").write_text(json.dumps(metrics, indent=2) + "\n")
    if args.save_npz:
        np.savez_compressed(
            output.with_suffix(".npz"),
            x=plane.x,
            y=plane.y,
            level=plane.level,
            time=plane.time,
            **{
                field: plane.level if field == "amr_level" else plane.values[field]
                for field in args.fields
            },
        )
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Locate a 3D extremum and mark its angular position on the first-layer Br map."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import data_paths, load_config
from simulation_analysis.snapshot import extract_axis_slice, inspect_snapshot


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dat", type=Path)
    parser.add_argument("--variable", default="temperature_K")
    parser.add_argument("--extreme", choices=("min", "max", "absmax"), default="max")
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    config = load_config(args.config)
    paths = data_paths(config)
    source = args.dat or (paths[-1] if paths else None)
    if source is None:
        raise SystemExit("No DAT selected")
    report = inspect_snapshot(source, config, variable=args.variable, extreme=args.extreme)
    location = report["extremum"]
    bottom = extract_axis_slice(
        source,
        config,
        axis="r",
        coordinate=float(config["grid"]["surface_radius_code"]),
    )
    theta, phi = bottom.coordinates
    latitude = (90.0 - np.degrees(theta))[::-1]
    longitude = np.degrees(phi)
    br = bottom.values["br_G"][::-1, :]
    limit = args.vmax or float(np.nanpercentile(np.abs(br), 99.5))
    figure, axis = plt.subplots(figsize=(11, 6.5), constrained_layout=True)
    image = axis.imshow(
        br,
        origin="lower",
        extent=[longitude[0], longitude[-1], latitude[0], latitude[-1]],
        aspect="auto",
        cmap="coolwarm",
        norm=TwoSlopeNorm(0.0, -limit, limit),
    )
    axis.scatter(
        location["phi_degree"],
        location["latitude_degree"],
        marker="x",
        s=100,
        linewidths=2.0,
        color="black",
        label=(
            f"{args.variable} {args.extreme}; h={location['height_Mm']:.3f} Mm, "
            f"L{location['level']}"
        ),
    )
    axis.set_xlabel("Longitude [deg]")
    axis.set_ylabel("Latitude [deg]")
    axis.legend(loc="best")
    figure.colorbar(image, ax=axis, label=r"First-layer $B_r$ [G]")
    axis.set_title(f"{Path(source).name}, t={report['time_code']:.6g}")
    output = args.output or Path(config["output_dir"]) / f"{args.variable}_{args.extreme}_on_bottom_br.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=int(config["plot"]["dpi"]), bbox_inches="tight")
    plt.close(figure)
    output.with_suffix(".json").write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(location, indent=2), flush=True)
    print(output.resolve(), flush=True)


if __name__ == "__main__":
    main()

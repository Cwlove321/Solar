#!/usr/bin/env python3
"""Plot dt, momentum extrema, density, and internal-energy diagnostics."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import load_config
from simulation_analysis.logs import read_rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    config = load_config(args.config)
    rows = read_rows(config)
    if not rows:
        raise SystemExit("No per-step rows were parsed")

    time = np.asarray([float(row["time"]) for row in rows])
    dt = np.asarray([float(row["dt"]) for row in rows])
    figure, axes = plt.subplots(3, 1, figsize=(11, 11), constrained_layout=True, sharex=True)
    axes[0].plot(time, dt, color="#1d5f8a", linewidth=1.0)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("dt [code units]")

    for name, label, color in (
        ("m1_extreme", "m1", "#c23b32"),
        ("m2_extreme", "m2", "#2769a3"),
        ("m3_extreme", "m3", "#23845d"),
        ("momentum_extreme", "max |m|", "black"),
    ):
        if name in rows[0]:
            axes[1].plot(time, [float(row[name]) for row in rows], label=label, color=color, linewidth=1.0)
    axes[1].set_ylabel("Momentum extrema [code units]")
    axes[1].legend(ncol=4)

    if "rho_max" in rows[0]:
        axes[2].plot(time, [float(row["rho_max"]) for row in rows], color="#23845d", label="rho max")
    if "e_max" in rows[0]:
        twin = axes[2].twinx()
        twin.plot(time, [float(row["e_max"]) for row in rows], color="#713c97", label="e max")
        twin.set_ylabel("Internal energy max [code units]", color="#713c97")
    axes[2].set_ylabel("Density max [code units]", color="#23845d")
    axes[2].set_xlabel("Simulation time [code units]")

    for axis in axes:
        axis.grid(alpha=0.2, linewidth=0.5)
    top = axes[0].secondary_xaxis(
        "top",
        functions=(
            lambda value: value * float(config["units"]["time_s"]) / 60.0,
            lambda value: value * 60.0 / float(config["units"]["time_s"]),
        ),
    )
    top.set_xlabel("Physical time [min]")
    figure.suptitle(f"Run history through it={int(rows[-1]['it'])}, t={time[-1]:.6g}")
    output = args.output or Path(config["output_dir"]) / "run_history.png"
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=int(config["plot"]["dpi"]), bbox_inches="tight")
    plt.close(figure)
    print(output.resolve(), flush=True)


if __name__ == "__main__":
    main()

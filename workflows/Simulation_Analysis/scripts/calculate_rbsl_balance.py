#!/usr/bin/env python3
"""Calculate the RBSL balance-current profile in an external AMR field."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import load_config
from simulation_analysis.rbsl_balance import estimate_rbsl_balance


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--dat", type=Path)
    source.add_argument("--vtu", type=Path)
    parser.add_argument("--path-dat", type=Path, required=True)
    parser.add_argument("--path-summary", type=Path, required=True)
    parser.add_argument("--minor-radius-mm", type=float, required=True)
    parser.add_argument("--current-flux-mx", type=float)
    parser.add_argument("--helicity", choices=("negative", "positive"), default="negative")
    parser.add_argument("--central-fraction", nargs=2, type=float, default=(0.1, 0.9))
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--prefix", default="rbsl_balance")
    args = parser.parse_args()

    config = load_config(args.config)
    output_dir = args.output_dir or Path(config["output_dir"]) / "rbsl_balance"
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    path_info = json.loads(args.path_summary.expanduser().resolve().read_text())
    helicity_sign = -1 if args.helicity == "negative" else 1

    summary, profile = estimate_rbsl_balance(
        args.dat if args.dat is not None else args.vtu,
        args.path_dat,
        path_info,
        config,
        minor_radius_mm=args.minor_radius_mm,
        current_flux_mx=args.current_flux_mx,
        helicity_sign=helicity_sign,
        central_fraction=tuple(args.central_fraction),
    )

    summary_path = output_dir / f"{args.prefix}_summary.json"
    csv_path = output_dir / f"{args.prefix}_profile.csv"
    plot_path = output_dir / f"{args.prefix}_profile.png"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    columns = list(profile)
    np.savetxt(
        csv_path,
        np.column_stack([profile[name] for name in columns]),
        delimiter=",",
        header=",".join(columns),
        comments="",
    )

    import matplotlib.pyplot as plt

    fraction = profile["fraction"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True, constrained_layout=True)
    axes[0].plot(fraction, profile["b_strap_gauss"], lw=1.2)
    axes[0].axhline(summary["b_strap_gauss_global_mean"], color="black", ls="--", lw=1.0)
    axes[0].set_ylabel(r"$B_{strap}$ [G]")

    axes[1].plot(fraction, profile["i_balance_code"], lw=1.2, label=r"$I_{balance}(s)$")
    axes[1].axhline(
        summary["i_balance_code_global_mean"],
        color="black",
        ls="--",
        lw=1.0,
        label="central weighted mean",
    )
    if "i_current_code" in summary:
        axes[1].axhline(
            summary["i_current_code"], color="tab:red", ls="--", lw=1.2, label="input flux"
        )
    axes[1].set_ylabel("normalized current")
    axes[1].legend()

    axes[2].plot(fraction, profile["flux_balance_local_mx"] / 1.0e20, lw=1.2)
    axes[2].axhline(
        summary["flux_balance_global_mx"] / 1.0e20, color="black", ls="--", lw=1.0
    )
    if args.current_flux_mx is not None:
        axes[2].axhline(args.current_flux_mx / 1.0e20, color="tab:red", ls="--", lw=1.2)
    axes[2].set_xlabel("axis arc-length fraction")
    axes[2].set_ylabel(r"local balance flux [$10^{20}$ Mx]")

    lower, upper = args.central_fraction
    for axis in axes:
        axis.axvspan(lower, upper, color="0.9", zorder=-2)
        axis.grid(alpha=0.25)
    fig.savefig(plot_path, dpi=int(config["plot"]["dpi"]))
    plt.close(fig)

    print(json.dumps(summary, indent=2), flush=True)
    print(output_dir, flush=True)


if __name__ == "__main__":
    main()

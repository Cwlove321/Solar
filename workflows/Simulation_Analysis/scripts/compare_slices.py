#!/usr/bin/env python3
"""Compare the same coordinate slice in two AMRVAC DAT snapshots."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import load_config
from simulation_analysis.plotting import draw_scalar, field_style
from simulation_analysis.snapshot import extract_axis_slice, to_plot_plane


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--dat-a", required=True, type=Path)
    parser.add_argument("--dat-b", required=True, type=Path)
    parser.add_argument("--axis", choices=("r", "theta", "phi"))
    parser.add_argument("--coordinate", type=float)
    parser.add_argument("--field", default="br_G")
    parser.add_argument("--relative", action="store_true")
    parser.add_argument("--difference-limit", type=float)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    config = load_config(args.config)

    planes = []
    for path in (args.dat_a, args.dat_b):
        section = extract_axis_slice(path, config, axis=args.axis, coordinate=args.coordinate)
        planes.append(to_plot_plane(section, config))
    first, second = planes
    if not np.allclose(first.x, second.x) or not np.allclose(first.y, second.y):
        raise RuntimeError("The two slices do not share the same plotting grid")
    a = first.level.astype(float) if args.field == "amr_level" else first.values[args.field]
    b = second.level.astype(float) if args.field == "amr_level" else second.values[args.field]
    if args.relative:
        difference = np.divide(
            b - a,
            np.abs(a),
            out=np.full_like(a, np.nan),
            where=np.abs(a) > 1.0e-30,
        )
        difference_label = f"Relative change in {args.field}"
    else:
        difference = b - a
        difference_label = f"Delta {args.field} (B-A)"

    combined = np.concatenate((a.ravel(), b.ravel()))
    label, cmap, norm = field_style(args.field, combined, config)
    finite = np.abs(difference[np.isfinite(difference)])
    limit = args.difference_limit or (
        float(np.percentile(finite, 99.0)) if finite.size else 1.0
    )
    limit = max(limit, 1.0e-30)
    difference_norm = TwoSlopeNorm(0.0, -limit, limit)

    figure, axes = plt.subplots(
        1, 3, figsize=(19, 5.5), constrained_layout=True, sharex=True, sharey=True
    )
    images = (
        draw_scalar(axes[0], first, a, cmap=cmap, norm=norm),
        draw_scalar(axes[1], second, b, cmap=cmap, norm=norm),
        draw_scalar(axes[2], second, difference, cmap="coolwarm", norm=difference_norm),
    )
    titles = (f"A: t={first.time:.6g}", f"B: t={second.time:.6g}", difference_label)
    for axis, title in zip(axes, titles):
        axis.set_title(title)
        axis.set_xlabel(first.x_name)
        axis.set_ylabel(first.y_name)
        axis.set_aspect(
            "equal" if config["slice"]["equal_physical_scale"] else "auto",
            adjustable="box",
        )
    figure.colorbar(images[0], ax=axes[:2], label=label, shrink=0.88)
    figure.colorbar(images[2], ax=axes[2], label=difference_label, shrink=0.88)
    figure.suptitle(first.fixed_label)
    output = args.output or Path(config["output_dir"]) / (
        f"compare_{args.field}_{args.dat_a.stem.rsplit('_', 1)[-1]}_"
        f"{args.dat_b.stem.rsplit('_', 1)[-1]}.png"
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output, dpi=int(config["plot"]["dpi"]), bbox_inches="tight")
    plt.close(figure)
    metrics = {
        "dat_a": str(args.dat_a.resolve()),
        "dat_b": str(args.dat_b.resolve()),
        "field": args.field,
        "relative": args.relative,
        "difference_min": float(np.nanmin(difference)),
        "difference_max": float(np.nanmax(difference)),
        "difference_rms": float(np.sqrt(np.nanmean(difference**2))),
        "difference_limit": limit,
        "output": str(output.resolve()),
    }
    output.with_suffix(".json").write_text(json.dumps(metrics, indent=2) + "\n")
    print(json.dumps(metrics, indent=2), flush=True)


if __name__ == "__main__":
    main()

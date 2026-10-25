#!/usr/bin/env python3
"""Plot several DAT times with shared limits and true local physical aspect."""

from __future__ import annotations

import argparse
import gc
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import data_paths, load_config
from simulation_analysis.plotting import add_fieldlines, draw_scalar, field_style
from simulation_analysis.snapshot import extract_axis_slice, to_plot_plane


def parse_layout(value: str) -> tuple[int, int]:
    try:
        rows, columns = (int(part) for part in value.lower().split("x", 1))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("layout must look like 6x1 or 2x3") from exc
    if rows <= 0 or columns <= 0:
        raise argparse.ArgumentTypeError("layout dimensions must be positive")
    return rows, columns


def select_paths(paths: list[Path], first: int, count: int, stride: int) -> list[Path]:
    selected = paths[first::stride]
    return selected[:count] if count > 0 else selected


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--axis", choices=("r", "theta", "phi"))
    parser.add_argument("--coordinate", type=float)
    parser.add_argument("--fields", nargs="+", default=["temperature_K", "number_density_cm3"])
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--count", type=int, default=6)
    parser.add_argument("--stride", type=int, default=1)
    parser.add_argument("--layout", type=parse_layout, default=(6, 1))
    parser.add_argument("--fieldlines", action="store_true")
    parser.add_argument("--output-prefix", type=Path)
    args = parser.parse_args()

    config = load_config(args.config)
    selected = select_paths(data_paths(config), args.first, args.count, args.stride)
    if not selected:
        raise SystemExit("No DAT files selected")
    rows, columns = args.layout
    if rows * columns < len(selected):
        raise SystemExit(f"layout {rows}x{columns} has fewer panels than {len(selected)} files")

    planes = []
    for position, path in enumerate(selected, start=1):
        section = extract_axis_slice(path, config, axis=args.axis, coordinate=args.coordinate)
        plane = to_plot_plane(section, config)
        keep = {name: plane.values[name] for name in args.fields}
        if args.fieldlines:
            keep.update({name: plane.values[name] for name in ("b1", "b2", "b3")})
        plane.values = keep
        planes.append(plane)
        print(f"[{position}/{len(selected)}] {path.name}", flush=True)
        del section
        gc.collect()

    prefix = args.output_prefix or Path(config["output_dir"]) / (
        f"series_{config['slice']['axis']}_{config['slice']['coordinate']}"
    )
    prefix = prefix.expanduser().resolve()
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for field in args.fields:
        all_values = np.concatenate([plane.values[field].ravel() for plane in planes])
        label, cmap, norm = field_style(field, all_values, config)
        width = max(7.0, 6.2 * columns)
        height = max(5.0, 4.5 * rows)
        figure, axes = plt.subplots(
            rows,
            columns,
            figsize=(width, height),
            constrained_layout=True,
            squeeze=False,
            sharex=True,
            sharey=True,
        )
        image = None
        for panel, plane, path in zip(axes.flat, planes, selected):
            image = draw_scalar(panel, plane, plane.values[field], cmap=cmap, norm=norm)
            panel.set_aspect(
                "equal" if config["slice"]["equal_physical_scale"] else "auto",
                adjustable="box",
            )
            if args.fieldlines:
                add_fieldlines(panel, plane)
            panel.set_title(f"{path.stem.rsplit('_', 1)[-1]}, t={plane.time:.5g}", loc="left")
            panel.set_xlabel(plane.x_name)
            panel.set_ylabel(plane.y_name)
            panel.grid(alpha=0.08, linewidth=0.3)
        for panel in axes.flat[len(planes) :]:
            panel.remove()
        if image is None:
            raise RuntimeError("No panels drawn")
        figure.colorbar(image, ax=list(axes.flat[: len(planes)]), label=label, shrink=0.88)
        figure.suptitle(f"{planes[0].fixed_label}: {label} evolution")
        output = Path(f"{prefix}_{field}.png")
        figure.savefig(output, dpi=int(config["plot"]["dpi"]), bbox_inches="tight")
        plt.close(figure)
        print(output, flush=True)


if __name__ == "__main__":
    main()

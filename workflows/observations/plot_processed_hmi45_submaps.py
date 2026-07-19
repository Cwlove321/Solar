#!/usr/bin/env python3
"""Plot processed HMI 45-second derotated LOS magnetogram submaps."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from sunpy.map import Map
from tqdm import tqdm

from solarutils.plotting import draw_submap_box, plot_hmi_los


DEFAULT_INPUT_DIR = Path("cases/20241008/products/processed/observations/2300_2400/sdo/hmi_45s/los_derotate_submap_pixel")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=None,
        help="Output directory. Defaults to a sibling directory named `<input>_figures`.",
    )
    parser.add_argument(
        "--box",
        nargs=4,
        type=float,
        action="append",
        default=None,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Draw a box in arcsec. Repeat for multiple boxes.",
    )
    parser.add_argument("--box-labels", nargs="*", default=None)
    parser.add_argument("--vmax", type=float, default=800.0)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    input_dir = args.input_dir
    output_dir = args.output_dir or input_dir.parent / f"{input_dir.name}_figures"
    files = sorted(path for path in input_dir.glob("*_submap.fits") if not path.name.startswith("._"))
    if not files:
        files = sorted(path for path in input_dir.glob("*.fits") if not path.name.startswith("._"))
    if not files:
        raise FileNotFoundError(f"No processed HMI FITS files found in {input_dir}")

    maps = [Map(path) for path in files]
    boxes = [tuple(box) for box in (args.box or [])]
    labels = args.box_labels
    if labels is not None and len(labels) != len(boxes):
        raise ValueError("--box-labels must have the same length as --box")
    output_dir.mkdir(parents=True, exist_ok=True)
    single_count = 0
    for path, smap in tqdm(list(zip(files, maps, strict=True)), desc="Plot HMI", unit="png"):
        output = output_dir / f"{path.stem}.png"
        if output.exists() and not args.overwrite:
            continue
        plot_single(smap, output, vmax=args.vmax, dpi=args.dpi, boxes=boxes, labels=labels)
        single_count += 1
        print(output, flush=True)

    sequence_output = output_dir / "hmi45_submap_sequence.png"
    if args.overwrite or not sequence_output.exists():
        plot_sequence(maps, sequence_output, vmax=args.vmax, dpi=args.dpi, boxes=boxes, labels=labels)
        print(sequence_output, flush=True)

    print(f"single_png={single_count}", flush=True)


def plot_single(smap, output: Path, *, vmax: float, dpi: int, boxes=None, labels=None) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.4), dpi=dpi)
    plot_hmi_los(smap, ax=ax, vmax=vmax, major_tick=100)
    draw_boxes(ax, boxes or [], labels)
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.18, top=0.86)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_sequence(maps: list, output: Path, *, vmax: float, dpi: int, boxes=None, labels=None) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(12.2, 6.1), dpi=dpi)
    axes = axes.ravel()
    for ax, smap in zip(axes, maps, strict=False):
        plot_hmi_los(smap, ax=ax, vmax=vmax, major_tick=200)
        draw_boxes(ax, boxes or [], labels)
        ax.set_xlabel("")
        ax.set_ylabel("")
    for ax in axes[len(maps) :]:
        ax.axis("off")
    fig.suptitle("HMI 45s LOS derotate submap", y=0.98)
    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.07, top=0.91, wspace=0.18, hspace=0.28)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def draw_boxes(ax, boxes: list[tuple[float, float, float, float]], labels=None) -> None:
    colors = ["yellow", "lime", "cyan", "magenta"]
    for index, box in enumerate(boxes):
        color = colors[index % len(colors)]
        draw_submap_box(ax, box, edgecolor=color, linewidth=1.8)
        if labels:
            x_min, _, _, y_max = box
            ax.text(x_min, y_max, f" {labels[index]}", color=color, fontsize=9, va="bottom", ha="left")


if __name__ == "__main__":
    main()

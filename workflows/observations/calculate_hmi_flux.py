#!/usr/bin/env python3
"""Calculate magnetic flux in one or more HMI submap regions."""

from __future__ import annotations

import argparse
from datetime import datetime
from glob import glob
from pathlib import Path

import matplotlib.dates as mdates
import matplotlib.pyplot as plt
import numpy as np
from sunpy.map import Map

from solarutils.fields import VALID_FLUX_MODES, magnetic_flux_series
from solarutils.plotting import draw_submap_box, plot_hmi_los


DEFAULT_INPUT = Path(
    "cases/20241008/products/processed/observations/2300_2400/sdo/hmi_45s/sub_m45"
)
DEFAULT_OUTPUT = Path("cases/20241008/products/processed/observations/2300_2400/sdo/hmi_45s/flux/hmi45_flux.npz")
DEFAULT_FIGURE_DIR = Path("cases/20241008/products/figures/observations/2300_2400/sdo/hmi_45s/flux")
DEFAULT_EXTENT = (200.0, 500.0, 220.0, 420.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--input",
        type=Path,
        default=DEFAULT_INPUT,
        help="FITS file, directory, or glob pattern for HMI submaps",
    )
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--figure-dir", type=Path, default=DEFAULT_FIGURE_DIR)
    parser.add_argument(
        "--extent",
        nargs=4,
        type=float,
        action="append",
        default=None,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Region in arcsec. Repeat for multiple regions.",
    )
    parser.add_argument("--region-labels", nargs="*", default=None)
    parser.add_argument(
        "--modes",
        nargs="+",
        choices=VALID_FLUX_MODES,
        default=["unsigned", "positive", "negative"],
    )
    parser.add_argument("--vmax", type=float, default=800.0)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    paths = resolve_input_paths(args.input)
    extents = [tuple(extent) for extent in (args.extent or [DEFAULT_EXTENT])]
    if args.output.exists() and not args.overwrite:
        raise FileExistsError(f"Output exists, use --overwrite: {args.output}")

    result = magnetic_flux_series(
        paths,
        extents,
        modes=tuple(args.modes),
        region_labels=args.region_labels,
        output=args.output,
    )
    args.figure_dir.mkdir(parents=True, exist_ok=True)

    region_figure = args.figure_dir / f"{args.output.stem}_regions.png"
    curve_figure = args.figure_dir / f"{args.output.stem}_curves.png"
    plot_flux_regions(Map(paths[0]), extents, result["region_labels"], region_figure, vmax=args.vmax, dpi=args.dpi)
    plot_flux_curves(result, curve_figure, dpi=args.dpi)

    print(args.output, flush=True)
    print(region_figure, flush=True)
    print(curve_figure, flush=True)


def resolve_input_paths(input_path: Path) -> list[Path]:
    input_text = str(input_path)
    if any(char in input_text for char in "*?[]"):
        paths = [Path(path) for path in sorted(glob(input_text))]
    elif input_path.is_dir():
        paths = sorted(input_path.glob("*_submap.fits"))
        if not paths:
            paths = sorted(input_path.glob("*.fits"))
    else:
        paths = [input_path]

    paths = [path for path in paths if path.exists() and not path.name.startswith("._")]
    if not paths:
        raise FileNotFoundError(f"No HMI FITS files found for {input_path}")
    return paths


def plot_flux_regions(smap, extents, labels, output: Path, *, vmax: float, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(6.2, 3.8), dpi=dpi)
    plot_hmi_los(smap, ax=ax, vmax=vmax, major_tick=100)
    for extent, label in zip(extents, labels, strict=True):
        box = draw_submap_box(ax, extent, edgecolor="yellow", linewidth=1.8)
        x_min, x_max, y_min, y_max = extent
        ax.text(
            x_min,
            y_max,
            f" {label}",
            color="yellow",
            fontsize=10,
            va="bottom",
            ha="left",
            path_effects=[],
        )
        box.set_linestyle("--")
    fig.subplots_adjust(left=0.12, right=0.98, bottom=0.16, top=0.86)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_flux_curves(result: dict[str, np.ndarray], output: Path, *, dpi: int) -> None:
    times = [datetime.fromisoformat(str(value)) for value in result["time_iso"]]
    labels = result["region_labels"]
    modes = result["modes"]
    flux = result["flux_mx"]

    fig, axes = plt.subplots(len(modes), 1, figsize=(7.2, 2.4 * len(modes)), dpi=dpi, sharex=True)
    axes = np.atleast_1d(axes)
    for mode_index, (ax, mode) in enumerate(zip(axes, modes, strict=True)):
        for region_index, label in enumerate(labels):
            values = flux[:, region_index, mode_index]
            ax.plot(times, values, marker="o", linewidth=1.4, label=str(label))
        ax.set_ylabel(f"{mode} flux (Mx)")
        set_axis_limits_from_data(ax, flux[:, :, mode_index])
        ax.grid(alpha=0.25)
        ax.legend(loc="best", fontsize=8)
    axes[-1].set_xlabel("Time")
    axes[-1].xaxis.set_major_formatter(mdates.DateFormatter("%H:%M"))
    fig.suptitle("HMI magnetic flux", y=0.995)
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.12, top=0.92, hspace=0.25)
    fig.autofmt_xdate()
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def set_axis_limits_from_data(ax, values, padding: float = 0.08) -> None:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return
    ymin = float(np.nanmin(finite))
    ymax = float(np.nanmax(finite))
    if ymin == ymax:
        delta = abs(ymin) * padding or 1.0
    else:
        delta = (ymax - ymin) * padding
    ax.set_ylim(ymin - delta, ymax + delta)


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Generate a visual-test figure for AIA sharpening methods."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sunpy.map import Map

from solarutils.obs.maps import get_submap, multiscale_sharpen_map, sharpen_map
from solarutils.plotting import plot_aia_map


DEFAULT_INPUT = Path(
    "data/test/observations/20241008_2300_2400/sdo/aia/"
    "aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits"
)
DEFAULT_OUTPUT = Path("data/test/figures/observations/sdo/aia/193/sharpen_compare.png")
DEFAULT_EXTENT = (100.0, 800.0, 150.0, 500.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--sigma", type=float, default=1.0)
    parser.add_argument("--weight", type=float, default=1.5)
    parser.add_argument("--multiscale", action="store_true", help="Use multi-scale unsharp-mask sharpening")
    parser.add_argument(
        "--sigmas",
        nargs="+",
        type=float,
        default=(1.0, 2.0, 4.0),
        help="Gaussian scales for --multiscale",
    )
    parser.add_argument(
        "--weights",
        nargs="+",
        type=float,
        default=(0.8, 0.35, 0.15),
        help="Weights for --multiscale",
    )
    parser.add_argument(
        "--extent",
        nargs=4,
        type=float,
        default=DEFAULT_EXTENT,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Submap extent in arcsec",
    )
    args = parser.parse_args()

    original = get_submap(Map(args.input), tuple(args.extent))
    if args.multiscale:
        sharpened = multiscale_sharpen_map(original, sigmas=tuple(args.sigmas), weights=tuple(args.weights))
        method_label = "Multiscale sharpen"
        detail_label = f"sigmas={tuple(args.sigmas)}, weights={tuple(args.weights)}"
    else:
        sharpened = sharpen_map(original, sigma=args.sigma, weight=args.weight)
        method_label = "Sharpen"
        detail_label = f"sigma={args.sigma:g}, weight={args.weight:g}"
    difference = sharpened.data - original.data

    finite = np.isfinite(original.data)
    norm_max = float(np.nanpercentile(original.data[finite], 99.7))
    diff_limit = float(np.nanpercentile(np.abs(difference[np.isfinite(difference)]), 99.5))

    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.8), dpi=180)
    plot_aia_map(original, ax=axes[0], norm_max=norm_max, major_tick=100)
    axes[0].set_title(f"Original\n{axes[0].get_title()}")
    plot_aia_map(sharpened, ax=axes[1], norm_max=norm_max, major_tick=100)
    axes[1].set_title(f"{method_label}\n{axes[1].get_title()}")

    extent = (
        original.bottom_left_coord.Tx.to_value("arcsec"),
        original.top_right_coord.Tx.to_value("arcsec"),
        original.bottom_left_coord.Ty.to_value("arcsec"),
        original.top_right_coord.Ty.to_value("arcsec"),
    )
    image = axes[2].imshow(
        difference,
        origin="lower",
        extent=extent,
        cmap="RdBu_r",
        vmin=-diff_limit,
        vmax=diff_limit,
    )
    axes[2].set_title(f"{method_label} - Original")
    axes[2].set_xlabel("Solar-X (arcsec)")
    axes[2].set_ylabel("Solar-Y (arcsec)")
    fig.colorbar(image, ax=axes[2], fraction=0.046, pad=0.035, label="DN")

    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.16, top=0.82, wspace=0.28)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)

    print(args.output)
    print(detail_label)
    print(f"original range: {np.nanmin(original.data):.3f} to {np.nanmax(original.data):.3f}")
    print(f"sharpened range: {np.nanmin(sharpened.data):.3f} to {np.nanmax(sharpened.data):.3f}")
    print(f"difference mean/std: {np.nanmean(difference):.3f} / {np.nanstd(difference):.3f}")


if __name__ == "__main__":
    main()

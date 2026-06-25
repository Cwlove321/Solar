#!/usr/bin/env python3
"""Compare center-point and per-pixel HMI 45s derotation methods."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sunpy.map import Map

from solarutils.plotting import plot_hmi_los


DEFAULT_CENTER = Path(
    "data/test/processed/observations/20241008_2300_2400/sdo/hmi_45s/"
    "los_derotate_submap_center/20241008_235915_to_20241008_230045.hmi45_los_derotate_submap_center.fits"
)
DEFAULT_PIXEL = Path(
    "data/test/processed/observations/20241008_2300_2400/sdo/hmi_45s/"
    "los_derotate_submap_pixel/20241008_235915_to_20241008_230045.hmi45_los_derotate_submap_pixel.fits"
)
DEFAULT_OUTPUT = Path(
    "data/test/figures/observations/20241008_2300_2400/sdo/hmi_45s/"
    "compare_derotation_center_pixel.png"
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--center", type=Path, default=DEFAULT_CENTER)
    parser.add_argument("--pixel", type=Path, default=DEFAULT_PIXEL)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--vmax", type=float, default=800.0)
    parser.add_argument("--diff-vmax", type=float, default=None)
    args = parser.parse_args()

    center = Map(args.center)
    pixel = Map(args.pixel)
    difference = pixel.data - center.data
    finite = np.isfinite(difference)
    diff_vmax = args.diff_vmax
    if diff_vmax is None:
        diff_vmax = float(np.nanpercentile(np.abs(difference[finite]), 99.0))

    extent = (
        center.bottom_left_coord.Tx.to_value("arcsec"),
        center.top_right_coord.Tx.to_value("arcsec"),
        center.bottom_left_coord.Ty.to_value("arcsec"),
        center.top_right_coord.Ty.to_value("arcsec"),
    )

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.7), dpi=180)
    plot_hmi_los(center, ax=axes[0], vmax=args.vmax, major_tick=100)
    axes[0].set_title(f"Center\n{axes[0].get_title()}")
    plot_hmi_los(pixel, ax=axes[1], vmax=args.vmax, major_tick=100)
    axes[1].set_title(f"Pixel\n{axes[1].get_title()}")
    image = axes[2].imshow(
        difference,
        origin="lower",
        extent=extent,
        cmap="RdBu_r",
        vmin=-diff_vmax,
        vmax=diff_vmax,
    )
    axes[2].set_title("Pixel - Center")
    axes[2].set_xlabel("Solar-X (arcsec)")
    axes[2].set_ylabel("Solar-Y (arcsec)")
    fig.colorbar(image, ax=axes[2], fraction=0.046, pad=0.035, label="Gauss")

    fig.subplots_adjust(left=0.06, right=0.985, bottom=0.17, top=0.8, wspace=0.28)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)

    print(args.output)
    print(f"mean_abs_diff={np.nanmean(np.abs(difference)):.6g}")
    print(f"p95_abs_diff={np.nanpercentile(np.abs(difference[finite]), 95):.6g}")
    print(f"p99_abs_diff={np.nanpercentile(np.abs(difference[finite]), 99):.6g}")
    print(f"max_abs_diff={np.nanmax(np.abs(difference)):.6g}")


if __name__ == "__main__":
    main()

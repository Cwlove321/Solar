#!/usr/bin/env python3
"""Generate visual-test figures for AIA preprocessing."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from sunpy.map import Map

from solarutils.obs.maps import aia_prep, aia_psf_path, get_submap
from solarutils.plotting import plot_aia_map


DEFAULT_INPUT = Path(
    "cases/20241008/data/2300_2400/sdo/aia/"
    "aia.lev1_euv_12s.2024-10-08T230010Z.171.image_lev1.fits"
)
DEFAULT_OUTPUT = Path("cases/20241008/products/figures/observations/2300_2400/sdo/aia/171/preprocess_register_smoke.png")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--with-psf", action="store_true", help="Run PSF deconvolution; slow for full-disk maps")
    parser.add_argument("--psf-path", type=Path, default=None, help="Cached PSF file")
    parser.add_argument("--iterations", type=int, default=5, help="PSF deconvolution iterations")
    parser.add_argument("--with-pointing", action="store_true", help="Apply pointing correction; requires a pointing table")
    parser.add_argument(
        "--extent",
        nargs=4,
        type=float,
        default=None,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        help="Optional submap extent in arcsec for visual comparison",
    )
    args = parser.parse_args()

    level1 = Map(args.input)
    psf = args.psf_path
    if psf is None and args.with_psf:
        psf = aia_psf_path(int(level1.wavelength.to_value("Angstrom")))
    prepared = aia_prep(
        level1,
        correct_pointing=args.with_pointing,
        deconvolve_psf=args.with_psf,
        psf=psf,
        deconvolve_iterations=args.iterations,
        register_map=True,
    )
    left_map = level1
    right_map = prepared
    if args.extent is not None:
        left_map = get_submap(left_map, tuple(args.extent))
        right_map = get_submap(right_map, tuple(args.extent))

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.1), dpi=180)
    major_tick = 100 if args.extent is not None else 500
    plot_aia_map(left_map, ax=axes[0], major_tick=major_tick)
    axes[0].set_title(f"Level 1\n{axes[0].get_title()}")
    plot_aia_map(right_map, ax=axes[1], major_tick=major_tick)
    right_label = "PSF + Registered" if args.with_psf else "Registered"
    axes[1].set_title(f"{right_label}\n{axes[1].get_title()}")

    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.15, top=0.84, wspace=0.28)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()

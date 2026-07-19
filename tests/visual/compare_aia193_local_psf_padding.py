#!/usr/bin/env python3
"""Compare padded local AIA PSF deconvolution against full-frame deconvolution."""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from aiapy.psf import deconvolve
from sunpy.map import Map

from solarutils.obs.maps import deconvolve_submap_with_padding, ensure_even_shape, get_submap, load_aia_psf
from solarutils.plotting import plot_aia_map


DEFAULT_INPUT = Path(
    "cases/20241008/data/2300_2400/sdo/aia/"
    "aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits"
)
DEFAULT_OUTPUT = Path("cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/local_psf_padding_compare.png")
DEFAULT_EXTENT = (100.0, 800.0, 150.0, 500.0)
DEFAULT_PADDINGS = (0.0, 100.0, 200.0, 400.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--extent", nargs=4, type=float, default=DEFAULT_EXTENT)
    parser.add_argument("--paddings", nargs="+", type=float, default=list(DEFAULT_PADDINGS))
    parser.add_argument("--iterations", type=int, default=5)
    args = parser.parse_args()

    smap = Map(args.input)
    full_psf = load_aia_psf(193)

    full_start = time.perf_counter()
    full_decon = deconvolve(smap, psf=full_psf, iterations=args.iterations)
    full_seconds = time.perf_counter() - full_start
    reference = get_submap(full_decon, tuple(args.extent))
    reference = ensure_even_shape(reference)

    results = []
    for padding in args.paddings:
        start = time.perf_counter()
        local_decon = deconvolve_submap_with_padding(
            smap,
            tuple(args.extent),
            psf=full_psf,
            padding=padding,
            iterations=args.iterations,
        )
        seconds = time.perf_counter() - start
        local_target = ensure_even_shape(local_decon)
        stats = _difference_stats(local_target, reference)
        results.append((padding, seconds, local_target, stats))

    for padding, seconds, _, stats in results:
        print(
            f"padding={padding:5.0f} arcsec  "
            f"local={seconds:6.3f}s  full={full_seconds:6.3f}s  "
            f"speedup={full_seconds / seconds:5.1f}x  "
            f"mean={stats['mean_abs']:7.3f}  rms={stats['rms']:7.3f}  max={stats['max_abs']:7.3f}"
        )

    _plot_results(results, reference, args.output)


def _difference_stats(local_map, reference_map):
    ny = min(local_map.data.shape[0], reference_map.data.shape[0])
    nx = min(local_map.data.shape[1], reference_map.data.shape[1])
    difference = np.nan_to_num(local_map.data[:ny, :nx] - reference_map.data[:ny, :nx])
    return {
        "mean_abs": float(np.mean(np.abs(difference))),
        "rms": float(np.sqrt(np.mean(difference**2))),
        "max_abs": float(np.max(np.abs(difference))),
    }


def _plot_results(results, reference, output):
    fig, axes = plt.subplots(1, len(results) + 1, figsize=(3.2 * (len(results) + 1), 3.4), dpi=160)
    plot_aia_map(reference, ax=axes[0], major_tick=100)
    axes[0].set_title("Full-frame\n" + axes[0].get_title())
    for ax, (padding, _, local_map, stats) in zip(axes[1:], results):
        plot_aia_map(local_map, ax=ax, major_tick=100)
        ax.set_title(f"Pad {padding:.0f} arcsec\nmean |diff|={stats['mean_abs']:.2f}")
    fig.subplots_adjust(left=0.045, right=0.99, bottom=0.18, top=0.78, wspace=0.28)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(output)


if __name__ == "__main__":
    main()

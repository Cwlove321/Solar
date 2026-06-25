#!/usr/bin/env python3
"""Compare the maintained AIA 193 preprocessing with the 20241008 CME notebook flow."""

from __future__ import annotations

import argparse
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from aiapy.calibrate import register, update_pointing
from aiapy.calibrate.utils import get_pointing_table
import aiapy.psf
from sunpy.map import Map

from solarutils.obs.maps import aia_prep, aia_psf_path, get_submap
from solarutils.plotting import plot_aia_map


DEFAULT_INPUT = Path(
    "/Volumes/Seagate/case/20241009/data/aia/193/5min/"
    "aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits"
)
FALLBACK_INPUT = Path(
    "data/test/observations/20241008_2300_2400/sdo/aia/"
    "aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits"
)
DEFAULT_OUTPUT = Path("data/test/figures/observations/sdo/aia/193/notebook_flow_compare.png")
DEFAULT_EXTENT = (100.0, 800.0, 150.0, 500.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--extent", nargs=4, type=float, default=DEFAULT_EXTENT)
    args = parser.parse_args()

    input_path = args.input if args.input.exists() else FALLBACK_INPUT
    smap = Map(input_path)
    pointing_table = get_pointing_table("JSOC", time_range=(smap.date - 12 * u.h, smap.date + 12 * u.h))
    psf = aia_psf_path(193)

    notebook_map = register(
        aiapy.psf.deconvolve(
            update_pointing(smap, pointing_table=pointing_table),
            psf=np.load(psf)["psf"],
            iterations=args.iterations,
        )
    )
    maintained_map = aia_prep(
        smap,
        pointing_table=pointing_table,
        psf=psf,
        deconvolve_iterations=args.iterations,
    )

    difference = np.nanmax(np.abs(notebook_map.data - maintained_map.data))
    print(f"max_abs_difference={difference:.6g}")

    notebook_sub = get_submap(notebook_map, tuple(args.extent))
    maintained_sub = get_submap(maintained_map, tuple(args.extent))

    fig, axes = plt.subplots(1, 2, figsize=(8.2, 4.0), dpi=180)
    plot_aia_map(notebook_sub, ax=axes[0], major_tick=100)
    axes[0].set_title(f"Notebook flow\n{axes[0].get_title()}")
    plot_aia_map(maintained_sub, ax=axes[1], major_tick=100)
    axes[1].set_title(f"solarutils aia_prep\n{axes[1].get_title()}")
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.16, top=0.82, wspace=0.28)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.output, bbox_inches="tight")
    plt.close(fig)
    print(args.output)


if __name__ == "__main__":
    main()

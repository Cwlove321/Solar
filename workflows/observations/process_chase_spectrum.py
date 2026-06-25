#!/usr/bin/env python3
"""Plot CHASE Doppler channels and calculate a centroid Doppler velocity map."""

from __future__ import annotations

import argparse
import os
from pathlib import Path

_cache_dir = Path.cwd() / ".mpl_config"
(_cache_dir / "cache").mkdir(parents=True, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", str(_cache_dir))
os.environ.setdefault("XDG_CACHE_HOME", str(_cache_dir / "cache"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from solarutils.obs import chase_doppler_velocity_map, chase_map_at_velocity
from solarutils.plotting import plot_chase_channel, plot_chase_doppler_velocity


DEFAULT_INPUT = Path("/Volumes/Seagate/case/20241009/data/chase/RSM20241008T232012_0003_HA.fits")
DEFAULT_OUTPUT_DIR = Path("data/test/figures/observations/20241008_2300_2400/chase")
DEFAULT_EXTENT = (100.0, 800.0, 150.0, 500.0)
DEFAULT_VELOCITIES = (-40.0, -20.0, 0.0, 20.0, 40.0)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--extent", nargs=4, type=float, default=DEFAULT_EXTENT)
    parser.add_argument("--velocity", nargs="+", type=float, default=list(DEFAULT_VELOCITIES))
    parser.add_argument("--center-extent", nargs=4, type=float, default=(-100.0, 100.0, -100.0, 100.0))
    parser.add_argument("--line-window", nargs=2, type=float, default=None)
    parser.add_argument(
        "--no-solar-rotation-correction",
        action="store_true",
        help="Do not subtract the differential-rotation LOS velocity.",
    )
    parser.add_argument("--doppler-vmax", type=float, default=15.0)
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    extent = tuple(args.extent)
    center_extent = tuple(args.center_extent)
    line_window = None if args.line_window is None else tuple(args.line_window)

    first_map = chase_map_at_velocity(
        args.input,
        args.velocity[0],
        extent=extent,
        center_extent=center_extent,
        line_window=line_window,
    )
    reference_line_center = float(first_map.meta["linecen"])
    print(f"reference_line_center={reference_line_center:.6f} Angstrom", flush=True)

    for velocity in args.velocity:
        output = args.output_dir / f"{args.input.stem}_chase_{velocity:g}kms.png"
        if output.exists() and not args.overwrite:
            continue
        smap = chase_map_at_velocity(
            args.input,
            velocity,
            extent=extent,
            reference_line_center=reference_line_center,
            center_extent=center_extent,
            line_window=line_window,
        )
        plot_chase_channel(smap, output=output, dpi=args.dpi)
        plt.close("all")
        print(output, flush=True)

    velocity_map, line_center_map = chase_doppler_velocity_map(
        args.input,
        extent=extent,
        reference_line_center=reference_line_center,
        center_extent=center_extent,
        line_window=line_window,
        method="centroid",
        remove_solar_rotation=not args.no_solar_rotation_correction,
    )
    correction_suffix = "raw" if args.no_solar_rotation_correction else "rotcorr"
    data_output = args.output_dir / f"{args.input.stem}_doppler_centroid_{correction_suffix}.npz"
    fig_output = args.output_dir / f"{args.input.stem}_doppler_centroid_{correction_suffix}.png"
    np.savez_compressed(
        data_output,
        velocity_kms=velocity_map.data,
        line_center_angstrom=line_center_map,
        reference_line_center_angstrom=reference_line_center,
        solar_rotation_corrected=not args.no_solar_rotation_correction,
        extent_arcsec=np.asarray(extent, dtype=float),
        input=str(args.input),
    )
    plot_chase_doppler_velocity(velocity_map, output=fig_output, vmax=args.doppler_vmax, dpi=args.dpi)
    plt.close("all")
    print(data_output, flush=True)
    print(fig_output, flush=True)


if __name__ == "__main__":
    main()

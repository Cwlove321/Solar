#!/usr/bin/env python3
"""Plot processed AIA derotated submap FITS files as PNG quicklooks."""

from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sunpy.map import Map

from solarutils.obs.maps import multiscale_sharpen_map, sharpen_map
from solarutils.plotting import plot_aia_map


DEFAULT_INPUT_DIR = Path("cases/20241008/products/processed/observations/2300_2400/sdo/aia")
DEFAULT_OUTPUT_DIR = Path("cases/20241008/products/figures/observations/2300_2400/sdo/aia/processed")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--wavelengths", nargs="*", type=int, default=None)
    parser.add_argument(
        "--derotate-method",
        choices=("center", "pixel"),
        default="pixel",
        help="Processed derotation method to plot",
    )
    parser.add_argument(
        "--norm-source",
        choices=("first", "all"),
        default="first",
        help="Use the first frame or all frames to calculate percentile vmin/vmax",
    )
    parser.add_argument(
        "--percentiles",
        nargs=2,
        type=float,
        default=(1.0, 99.5),
        metavar=("VMIN_PERCENTILE", "VMAX_PERCENTILE"),
        help="Percentiles used to set common vmin/vmax for each wavelength",
    )
    parser.add_argument(
        "--sharpen",
        choices=("none", "single", "multiscale"),
        default="none",
        help="Optional sharpening applied before calculating vmin/vmax and plotting",
    )
    parser.add_argument("--sigma", type=float, default=1.0, help="Gaussian scale for --sharpen single")
    parser.add_argument("--weight", type=float, default=1.5, help="Unsharp-mask weight for --sharpen single")
    parser.add_argument(
        "--sigmas",
        nargs="+",
        type=float,
        default=(1.0, 2.0, 4.0),
        help="Gaussian scales for --sharpen multiscale",
    )
    parser.add_argument(
        "--weights",
        nargs="+",
        type=float,
        default=(0.8, 0.35, 0.15),
        help="Weights for --sharpen multiscale",
    )
    parser.add_argument("--dpi", type=int, default=180)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    processed_dir = processed_dir_name(args.derotate_method)
    files = sorted(args.input_dir.glob(f"*/{processed_dir}/*.fits"))
    if not files:
        raise FileNotFoundError(f"No processed AIA FITS files found in {args.input_dir} using {processed_dir}")

    grouped = group_by_wavelength(files)
    if args.wavelengths is not None:
        grouped = {f"{wave:04d}": grouped[f"{wave:04d}"] for wave in args.wavelengths if f"{wave:04d}" in grouped}
    if not grouped:
        raise ValueError("No processed AIA files matched the requested wavelengths")

    output_root = args.output_dir if args.sharpen == "none" else args.output_dir / f"{args.sharpen}_sharpen"
    single_count = 0
    sequence_count = 0
    for wave, wave_files in sorted(grouped.items()):
        maps = [Map(path) for path in wave_files]
        plot_maps = apply_sharpening(
            maps,
            method=args.sharpen,
            sigma=args.sigma,
            weight=args.weight,
            sigmas=tuple(args.sigmas),
            weights=tuple(args.weights),
        )
        norm_min, norm_max = percentile_norm_limits(
            plot_maps,
            source=args.norm_source,
            percentiles=tuple(args.percentiles),
        )
        print(
            f"[{wave}] sharpen={args.sharpen} norm_source={args.norm_source} "
            f"percentiles=({args.percentiles[0]:g}, {args.percentiles[1]:g}) "
            f"vmin={norm_min:.6g} vmax={norm_max:.6g}",
            flush=True,
        )
        wave_output_dir = output_root / wave / processed_dir
        wave_output_dir.mkdir(parents=True, exist_ok=True)

        for path, smap in zip(wave_files, plot_maps, strict=True):
            output = wave_output_dir / f"{path.stem}.png"
            if output.exists() and not args.overwrite:
                continue
            plot_single(smap, output, norm_min=norm_min, norm_max=norm_max, dpi=args.dpi)
            single_count += 1
            print(output, flush=True)

        sequence_output = output_root / wave / f"{wave}_{processed_dir}_sequence.png"
        if args.overwrite or not sequence_output.exists():
            plot_sequence(
                plot_maps,
                sequence_output,
                norm_min=norm_min,
                norm_max=norm_max,
                dpi=args.dpi,
                sharpen=args.sharpen,
                derotate_method=args.derotate_method,
            )
            sequence_count += 1
            print(sequence_output, flush=True)

    print(f"single_png={single_count}", flush=True)
    print(f"sequence_png={sequence_count}", flush=True)


def group_by_wavelength(files: list[Path]) -> dict[str, list[Path]]:
    grouped = defaultdict(list)
    for path in files:
        grouped[path.parts[-3]].append(path)
    return dict(grouped)


def processed_dir_name(method: str) -> str:
    if method == "center":
        return "prep_psf_derotate_submap"
    return f"prep_psf_derotate_submap_{method}"


def apply_sharpening(
    maps: list,
    *,
    method: str,
    sigma: float = 1.0,
    weight: float = 1.5,
    sigmas: tuple[float, ...] = (1.0, 2.0, 4.0),
    weights: tuple[float, ...] = (0.8, 0.35, 0.15),
):
    if method == "none":
        return maps
    if method == "single":
        return [sharpen_map(smap, sigma=sigma, weight=weight) for smap in maps]
    if method == "multiscale":
        return [multiscale_sharpen_map(smap, sigmas=sigmas, weights=weights) for smap in maps]
    raise ValueError(f"Unknown sharpening method: {method}")


def percentile_norm_limits(
    maps: list,
    *,
    source: str = "first",
    percentiles: tuple[float, float] = (1.0, 99.5),
) -> tuple[float, float]:
    selected_maps = maps[:1] if source == "first" else maps
    samples = [smap.data[np.isfinite(smap.data)] for smap in selected_maps]
    if not samples:
        return 1.0, 1.0
    data = np.concatenate(samples)
    vmin, vmax = np.nanpercentile(data, percentiles)
    if vmin == vmax:
        vmax = vmin + 1.0
    return float(vmin), float(vmax)


def plot_single(smap, output: Path, *, norm_min: float, norm_max: float, dpi: int) -> None:
    fig, ax = plt.subplots(figsize=(5.6, 3.4), dpi=dpi)
    plot_aia_map(smap, ax=ax, norm_min=norm_min, norm_max=norm_max, major_tick=100)
    fig.subplots_adjust(left=0.14, right=0.98, bottom=0.18, top=0.86)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


def plot_sequence(
    maps: list,
    output: Path,
    *,
    norm_min: float,
    norm_max: float,
    dpi: int,
    sharpen: str = "none",
    derotate_method: str = "pixel",
) -> None:
    fig, axes = plt.subplots(2, 4, figsize=(12.2, 6.1), dpi=dpi)
    axes = axes.ravel()
    for ax, smap in zip(axes, maps, strict=False):
        plot_aia_map(smap, ax=ax, norm_min=norm_min, norm_max=norm_max, major_tick=200)
        ax.set_xlabel("")
        ax.set_ylabel("")
    for ax in axes[len(maps) :]:
        ax.axis("off")
    wave = int(maps[0].wavelength.to_value("Angstrom"))
    suffix = "" if sharpen == "none" else f" + {sharpen} sharpen"
    method_suffix = "" if derotate_method == "center" else f" ({derotate_method})"
    fig.suptitle(f"AIA {wave} Å prep + PSF + derotate submap{method_suffix}{suffix}", y=0.98)
    fig.subplots_adjust(left=0.055, right=0.99, bottom=0.07, top=0.91, wspace=0.18, hspace=0.28)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

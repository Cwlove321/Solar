#!/usr/bin/env python3
"""Create movies from processed AIA quicklook PNG sequences."""

from __future__ import annotations

import argparse
from pathlib import Path

from solarutils.plotting import MOVIE_QUALITY_PRESETS, make_movie_from_images


DEFAULT_INPUT_DIR = Path("cases/20241008/products/figures/observations/2300_2400/sdo/aia/processed")
DEFAULT_OUTPUT_DIR = Path("cases/20241008/products/movies/observations/2300_2400/sdo/aia")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--wavelengths", nargs="*", type=int, default=None)
    parser.add_argument("--variant", default="none", choices=("none", "single_sharpen", "multiscale_sharpen"))
    parser.add_argument("--fps", type=float, default=3.0)
    parser.add_argument("--backend", default="auto", choices=("auto", "ffmpeg", "imageio"))
    parser.add_argument(
        "--quality",
        default="simple",
        choices=tuple(MOVIE_QUALITY_PRESETS),
        help="Movie encoding preset: simple makes compact previews; high reduces compression.",
    )
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    input_root = args.input_dir if args.variant == "none" else args.input_dir / args.variant
    wave_dirs = sorted(path for path in input_root.iterdir() if path.is_dir() and path.name.isdigit())
    if args.wavelengths is not None:
        requested = {f"{wave:04d}" for wave in args.wavelengths}
        wave_dirs = [path for path in wave_dirs if path.name in requested]
    if not wave_dirs:
        raise FileNotFoundError(f"No wavelength quicklook directories found in {input_root}")

    count = 0
    for wave_dir in wave_dirs:
        frame_dir = wave_dir / "prep_psf_derotate_submap"
        frames = sorted(path for path in frame_dir.glob("*.png") if not path.name.startswith("._"))
        if not frames:
            print(f"skip {wave_dir.name}: no PNG frames", flush=True)
            continue
        variant_suffix = "" if args.variant == "none" else f"_{args.variant}"
        quality_suffix = "" if args.quality == "simple" else f"_{args.quality}"
        output = args.output_dir / wave_dir.name / f"{wave_dir.name}_prep_psf_derotate_submap{variant_suffix}{quality_suffix}.mp4"
        make_movie_from_images(
            frames,
            output,
            fps=args.fps,
            backend=args.backend,
            quality=args.quality,
            overwrite=args.overwrite,
        )
        count += 1
        print(output, flush=True)

    print(f"movies={count}", flush=True)


if __name__ == "__main__":
    main()

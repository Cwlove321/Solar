#!/usr/bin/env python3
"""Generate visual-test figures for AIA context/submap plotting."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from sunpy.map import Map

from solarutils.obs.maps import get_submap
from solarutils.plotting import plot_aia_context_submap


DATA_DIR = Path("data/test/observations/20241008_2300_2400/sdo/aia")
OUTPUT_DIR = Path("data/test/figures/observations/sdo/aia")

EXAMPLES = {
    "aia171_context_submap": {
        "input": DATA_DIR / "aia.lev1_euv_12s.2024-10-08T230010Z.171.image_lev1.fits",
        "extent": (-450.0, 250.0, -550.0, 150.0),
        "output": Path("171/context_submap.png"),
    },
    "aia193_context_rectangle": {
        "input": DATA_DIR / "aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits",
        "extent": (-650.0, 350.0, -450.0, 50.0),
        "output": Path("193/context_rectangle.png"),
    },
}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--case", choices=[*EXAMPLES.keys(), "all"], default="all")
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    args = parser.parse_args()

    names = EXAMPLES.keys() if args.case == "all" else [args.case]
    for name in names:
        example = EXAMPLES[name]
        output = args.output_dir / example["output"]
        smap = Map(example["input"])
        submap = get_submap(smap, example["extent"])
        fig, _ = plot_aia_context_submap(smap, submap, example["extent"], output=output)
        plt.close(fig)
        print(output)


if __name__ == "__main__":
    main()

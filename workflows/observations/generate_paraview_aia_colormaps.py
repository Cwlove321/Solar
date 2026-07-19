#!/usr/bin/env python3
"""Generate ParaView-importable SDO/AIA color preset JSON."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
import numpy as np
import sunpy.visualization.colormaps  # noqa: F401 - registers AIA colormaps


DEFAULT_WAVELENGTHS = (94, 131, 171, 193, 211, 304, 335)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--wavelengths", nargs="+", type=int, default=DEFAULT_WAVELENGTHS)
    parser.add_argument("--samples", type=int, default=256)
    return parser.parse_args()


def make_preset(wavelength: int, samples: int) -> dict[str, object]:
    if samples < 2:
        raise ValueError("samples must be at least 2")
    name = f"sdoaia{wavelength}"
    try:
        cmap = matplotlib.colormaps[name]
    except KeyError as exc:
        raise ValueError(f"SunPy does not provide colormap {name!r}") from exc

    rgb_points: list[float] = []
    for position in np.linspace(0.0, 1.0, samples):
        red, green, blue, _ = cmap(float(position))
        rgb_points.extend(
            [
                round(float(position), 8),
                round(float(red), 8),
                round(float(green), 8),
                round(float(blue), 8),
            ]
        )

    return {
        "ColorSpace": "RGB",
        "Name": f"SDO AIA {wavelength} Angstrom (SunPy)",
        "NanColor": [0.5, 0.5, 0.5],
        "RGBPoints": rgb_points,
    }


def main() -> None:
    args = parse_args()
    presets = [make_preset(wavelength, args.samples) for wavelength in args.wavelengths]
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(presets, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()

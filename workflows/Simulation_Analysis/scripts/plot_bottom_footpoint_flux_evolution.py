#!/usr/bin/env python3
"""Plot the bottom-boundary Bz flux of the traced y=0 footprint.

The positive and negative-polarity magnitudes are plotted separately.  Their
half unsigned sum is the one-copy paired flux appropriate for a closed flux
tube:

    Phi_pair = (Phi_positive + abs(Phi_negative)) / 2.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DEFAULT_DIR = Path(
    "/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/y0_qsl/enclosed_region/"
    "evolution_logq4/bottom_footpoint_flux"
)
DEFAULT_CSV = DEFAULT_DIR / "y0_mask_bottom_flux_t0010_t0095.csv"
DEFAULT_OUTPUT = DEFAULT_DIR / "bottom_boundary_bz_flux_evolution_t0010_t0095"

POSITIVE_COLOR = "#CC79A7"
NEGATIVE_COLOR = "#0072B2"
PAIR_COLOR = "#202020"
FLUX_SCALE_MX = 1.0e19


def read_flux_csv(path: Path) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    if not rows:
        raise ValueError(f"No rows in {path}")
    required = {
        "snapshot",
        "time_minutes",
        "phi_bz_positive_Mx",
        "phi_bz_negative_signed_Mx",
        "phi_bz_paired_Mx",
    }
    missing = required.difference(rows[0])
    if missing:
        raise KeyError(f"{path} is missing columns {sorted(missing)}")

    result = {
        "snapshot": np.asarray([int(row["snapshot"]) for row in rows], dtype=int),
        "time_minutes": np.asarray(
            [float(row["time_minutes"]) for row in rows], dtype=float
        ),
        "positive_Mx": np.asarray(
            [float(row["phi_bz_positive_Mx"]) for row in rows], dtype=float
        ),
        "negative_signed_Mx": np.asarray(
            [float(row["phi_bz_negative_signed_Mx"]) for row in rows], dtype=float
        ),
        "paired_Mx": np.asarray(
            [float(row["phi_bz_paired_Mx"]) for row in rows], dtype=float
        ),
    }
    size = len(result["snapshot"])
    if any(len(values) != size for values in result.values()):
        raise ValueError("Flux columns have inconsistent lengths")
    if any(not np.all(np.isfinite(values)) for values in result.values()):
        raise ValueError("Flux CSV contains non-finite values")
    if not np.all(np.diff(result["snapshot"]) == 1):
        raise ValueError("Snapshots are not continuous")
    if not np.all(np.diff(result["time_minutes"]) > 0.0):
        raise ValueError("Times are not strictly increasing")
    if np.any(result["positive_Mx"] <= 0.0):
        raise ValueError("Positive-polarity flux must be positive")
    if np.any(result["negative_signed_Mx"] >= 0.0):
        raise ValueError("Negative-polarity flux must be negative")

    recomputed_pair = 0.5 * (
        result["positive_Mx"] + np.abs(result["negative_signed_Mx"])
    )
    if not np.allclose(recomputed_pair, result["paired_Mx"], rtol=1.0e-12):
        raise ValueError("Stored paired flux does not equal the half unsigned sum")
    return result


def plot_flux(data: dict[str, np.ndarray], output_base: Path, dpi: int) -> None:
    time = data["time_minutes"]
    positive = data["positive_Mx"] / FLUX_SCALE_MX
    negative = np.abs(data["negative_signed_Mx"]) / FLUX_SCALE_MX
    paired = data["paired_Mx"] / FLUX_SCALE_MX

    figure, axis = plt.subplots(figsize=(10.2, 6.2), constrained_layout=True)
    figure.patch.set_facecolor("white")
    axis.set_facecolor("white")
    axis.fill_between(
        time,
        np.minimum(positive, negative),
        np.maximum(positive, negative),
        color="#777777",
        alpha=0.16,
        linewidth=0.0,
        label="Positive-negative difference",
        zorder=1,
    )
    axis.plot(
        time,
        positive,
        color=POSITIVE_COLOR,
        linewidth=2.4,
        label=r"$\Phi_{B_z,+}$",
        zorder=3,
    )
    axis.plot(
        time,
        negative,
        color=NEGATIVE_COLOR,
        linewidth=2.4,
        linestyle="--",
        label=r"$|\Phi_{B_z,-}|$",
        zorder=3,
    )
    axis.plot(
        time,
        paired,
        color=PAIR_COLOR,
        linewidth=2.8,
        label=r"$\Phi_{B_z,\mathrm{pair}}=(\Phi_+ + |\Phi_-|)/2$",
        zorder=4,
    )

    axis.set_xlabel("Time [min]", fontsize=16)
    axis.set_ylabel(r"Bottom-boundary magnetic flux [$10^{19}$ Mx]", fontsize=16)
    axis.set_xlim(float(time.min()), float(time.max()))
    axis.set_ylim(0.0, 1.07 * float(max(positive.max(), negative.max())))
    axis.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        labelsize=13,
        length=6,
        width=1.1,
    )
    axis.grid(color="#BBBBBB", alpha=0.28, linewidth=0.8)
    axis.legend(
        loc="upper right",
        frameon=False,
        fontsize=13,
        handlelength=3.0,
    )
    for spine in axis.spines.values():
        spine.set_linewidth(1.1)

    output_base.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(output_base.with_suffix(".png"), dpi=dpi, facecolor="white")
    figure.savefig(output_base.with_suffix(".pdf"), facecolor="white")
    plt.close(figure)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--output-base", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--dpi", type=int, default=240)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.dpi <= 0:
        raise SystemExit("--dpi must be positive")
    data = read_flux_csv(args.input_csv)
    plot_flux(data, args.output_base, args.dpi)
    print(args.output_base.with_suffix(".png"))
    print(args.output_base.with_suffix(".pdf"))


if __name__ == "__main__":
    main()

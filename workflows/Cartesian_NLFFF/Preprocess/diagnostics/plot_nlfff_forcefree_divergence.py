#!/usr/bin/env python3
"""Plot NLFFF force-free and divergence metrics from an AMRVAC mflog CSV."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    return parser.parse_args()


def read_mflog(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, skipinitialspace=True)
        rows = list(reader)
    iteration = np.asarray([float(row["itmf"]) for row in rows])
    divergence = np.asarray([float(row["<f_i>"]) for row in rows])
    force_free = np.asarray([float(row["<CW sin theta>"]) for row in rows])
    return iteration, force_free, divergence


def main() -> None:
    args = parse_args()
    iteration, force_free, divergence = read_mflog(args.csv)

    fig, axes = plt.subplots(2, 1, figsize=(9.2, 6.8), sharex=True)
    force_color = "#2878B5"
    divergence_color = "#D65F0E"

    axes[0].plot(iteration, force_free, color=force_color, linewidth=1.25)
    axes[0].scatter([iteration[0], iteration[-1]], [force_free[0], force_free[-1]], color=force_color, s=20, zorder=3)
    axes[0].set_ylabel(r"Force-free factor $\langle CW\sin\theta\rangle$")
    axes[0].grid(alpha=0.25, linewidth=0.6)
    axes[0].annotate(f"{force_free[0]:.4f}", (iteration[0], force_free[0]), xytext=(8, -14), textcoords="offset points")
    axes[0].annotate(f"{force_free[-1]:.4f}", (iteration[-1], force_free[-1]), xytext=(-44, 8), textcoords="offset points")

    axes[1].plot(iteration, divergence, color=divergence_color, linewidth=1.0)
    axes[1].scatter([iteration[0], iteration[-1]], [divergence[0], divergence[-1]], color=divergence_color, s=20, zorder=3)
    axes[1].set_ylabel(r"Divergence factor $\langle f_i\rangle$")
    axes[1].set_xlabel("Magnetic-friction iteration")
    axes[1].ticklabel_format(axis="y", style="sci", scilimits=(-3, 3), useMathText=True)
    axes[1].grid(alpha=0.25, linewidth=0.6)
    axes[1].annotate(f"{divergence[0]:.3e}", (iteration[0], divergence[0]), xytext=(8, 8), textcoords="offset points")
    axes[1].annotate(f"{divergence[-1]:.3e}", (iteration[-1], divergence[-1]), xytext=(-72, 8), textcoords="offset points")

    fig.suptitle("NLFFF relaxation diagnostics")
    fig.tight_layout()
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

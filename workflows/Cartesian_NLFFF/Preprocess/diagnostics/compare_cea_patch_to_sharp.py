#!/usr/bin/env python
"""Compare a generated CEA patch with official SHARP CEA components."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from sunpy.map import Map


COMPONENTS = ("Br", "Bt", "Bp")


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    patch_dir = Path(args.patch_dir)
    official_dir = Path(args.official_dir)
    output_dir = Path(args.output_dir) if args.output_dir else patch_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    official = {component: read_component(find_component_file(official_dir, component)) for component in COMPONENTS}
    generated = {component: read_component(patch_dir / f"{component}.fits") for component in COMPONENTS}

    assert_matching_shapes(official, generated)
    strong_mask = np.isfinite(official["Br"]) & (np.abs(official["Br"]) >= args.strong_threshold)

    metrics = {
        "official_dir": str(official_dir),
        "patch_dir": str(patch_dir),
        "shape": list(official["Br"].shape),
        "strong_threshold_g": args.strong_threshold,
        "components": {},
    }
    for component in COMPONENTS:
        metrics["components"][component] = {
            "all_finite": component_metrics(official[component], generated[component]),
            "strong_br_region": component_metrics(official[component], generated[component], extra_mask=strong_mask),
            "best_integer_shift": best_integer_shift(official[component], generated[component], max_shift=args.max_shift),
        }

    metrics_path = output_dir / "official_sharp_comparison_metrics.json"
    metrics_path.write_text(json.dumps(metrics, indent=2), encoding="utf-8")
    plot_component_comparison(output_dir / "official_vs_make_cea_patch_components.png", official, generated, vmax=args.vmax)
    plot_scatter(output_dir / "official_vs_make_cea_patch_scatter_strong.png", official, generated, strong_mask)

    print(json.dumps(metrics, indent=2))
    print(f"Saved comparison products in {output_dir}")
    return 0


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patch-dir", required=True, type=Path, help="Directory containing generated Br/Bt/Bp FITS.")
    parser.add_argument("--official-dir", required=True, type=Path, help="Directory containing official SHARP Br/Bt/Bp FITS.")
    parser.add_argument("--output-dir", type=Path, default=None, help="Directory for comparison outputs; defaults to patch-dir.")
    parser.add_argument("--strong-threshold", type=float, default=100.0, help="Use |official Br| above this value for strong-region metrics.")
    parser.add_argument("--max-shift", type=int, default=5, help="Maximum integer pixel shift searched in each direction.")
    parser.add_argument("--vmax", type=float, default=1500.0, help="Symmetric image color limit in G.")
    return parser.parse_args(argv)


def find_component_file(directory: Path, component: str) -> Path:
    direct = directory / f"{component}.fits"
    if direct.exists():
        return direct
    matches = sorted(directory.glob(f"*{component}.fits"))
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected exactly one {component}.fits in {directory}, found {len(matches)}")
    return matches[0]


def read_component(path: Path) -> np.ndarray:
    return np.asarray(Map(path).data, dtype=float)


def assert_matching_shapes(official: dict[str, np.ndarray], generated: dict[str, np.ndarray]) -> None:
    for component in COMPONENTS:
        if official[component].shape != generated[component].shape:
            raise ValueError(
                f"{component} shape mismatch: official {official[component].shape}, "
                f"generated {generated[component].shape}"
            )


def component_metrics(reference, candidate, *, extra_mask=None):
    reference = np.asarray(reference, dtype=float)
    candidate = np.asarray(candidate, dtype=float)
    mask = np.isfinite(reference) & np.isfinite(candidate)
    if extra_mask is not None:
        mask &= np.asarray(extra_mask, dtype=bool)
    if not np.any(mask):
        return {"count": 0}

    ref = reference[mask]
    cand = candidate[mask]
    diff = cand - ref
    corr = np.corrcoef(ref, cand)[0, 1] if ref.size > 1 else np.nan
    slope, intercept = np.polyfit(ref, cand, deg=1) if ref.size > 1 else (np.nan, np.nan)
    return {
        "count": int(ref.size),
        "corr": finite_float(corr),
        "rmse_g": finite_float(np.sqrt(np.mean(diff**2))),
        "mae_g": finite_float(np.mean(np.abs(diff))),
        "bias_g": finite_float(np.mean(diff)),
        "p95_abs_g": finite_float(np.percentile(np.abs(diff), 95.0)),
        "slope_generated_vs_official": finite_float(slope),
        "intercept_generated_vs_official_g": finite_float(intercept),
    }


def best_integer_shift(reference, candidate, *, max_shift: int):
    best = {"dy": 0, "dx": 0, "corr": -np.inf}
    for dy in range(-max_shift, max_shift + 1):
        for dx in range(-max_shift, max_shift + 1):
            ref_view, cand_view = overlapping_views(reference, candidate, dy, dx)
            mask = np.isfinite(ref_view) & np.isfinite(cand_view)
            if np.count_nonzero(mask) < 2:
                continue
            corr = np.corrcoef(ref_view[mask], cand_view[mask])[0, 1]
            if np.isfinite(corr) and corr > best["corr"]:
                best = {"dy": int(dy), "dx": int(dx), "corr": float(corr)}
    return best


def overlapping_views(reference, candidate, dy: int, dx: int):
    y_ref = slice(max(0, dy), reference.shape[0] + min(0, dy))
    x_ref = slice(max(0, dx), reference.shape[1] + min(0, dx))
    y_cand = slice(max(0, -dy), candidate.shape[0] + min(0, -dy))
    x_cand = slice(max(0, -dx), candidate.shape[1] + min(0, -dx))
    return reference[y_ref, x_ref], candidate[y_cand, x_cand]


def plot_component_comparison(path: Path, official: dict[str, np.ndarray], generated: dict[str, np.ndarray], *, vmax: float):
    fig, axes = plt.subplots(3, 3, figsize=(11, 9), dpi=160, constrained_layout=True)
    for row, component in enumerate(COMPONENTS):
        panels = [
            (official[component], f"Official {component}", vmax),
            (generated[component], f"make_cea_patch {component}", vmax),
            (generated[component] - official[component], f"Difference {component}", vmax / 5.0),
        ]
        for ax, (data, title, limit) in zip(axes[row], panels, strict=True):
            im = ax.imshow(data, origin="lower", cmap="RdBu_r", vmin=-limit, vmax=limit)
            ax.set_title(title)
            ax.set_axis_off()
            fig.colorbar(im, ax=ax, shrink=0.72)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def plot_scatter(path: Path, official: dict[str, np.ndarray], generated: dict[str, np.ndarray], strong_mask):
    fig, axes = plt.subplots(1, 3, figsize=(11, 3.6), dpi=160, constrained_layout=True)
    for ax, component in zip(axes, COMPONENTS, strict=True):
        mask = strong_mask & np.isfinite(official[component]) & np.isfinite(generated[component])
        ref = official[component][mask]
        cand = generated[component][mask]
        ax.hexbin(ref, cand, gridsize=80, bins="log", cmap="viridis", mincnt=1)
        lo = np.nanpercentile(np.concatenate([ref, cand]), 1.0)
        hi = np.nanpercentile(np.concatenate([ref, cand]), 99.0)
        ax.plot([lo, hi], [lo, hi], color="crimson", linewidth=1.0)
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_title(component)
        ax.set_xlabel("Official SHARP [G]")
        ax.set_ylabel("make_cea_patch [G]")
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def finite_float(value):
    return float(value) if np.isfinite(value) else None


if __name__ == "__main__":
    raise SystemExit(main())

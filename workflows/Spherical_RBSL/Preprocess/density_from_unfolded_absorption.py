#!/usr/bin/env python
"""Map a single-channel unfolded EUV absorption image to relative density."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import gaussian_filter


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True)
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--background-y-min-mm", type=float, default=15.0)
    parser.add_argument(
        "--background-intensity-i0",
        type=float,
        help="Fixed external-region mean intensity. When provided, I0 is not recomputed from the unfolded map.",
    )
    parser.add_argument("--background-description", default="")
    parser.add_argument(
        "--threshold-fraction",
        type=float,
        help="Set tau_threshold to this fraction of tau_max instead of estimating it from unfolded sidebands.",
    )
    parser.add_argument("--smooth-mm", type=float, default=1.0)
    parser.add_argument("--maximum-density-rho0", type=float, default=50.0)
    args = parser.parse_args(argv)

    source = np.load(args.input)
    intensity = np.asarray(source["intensity"], dtype=float)
    x_mm = np.asarray(source["axis_distance_mm"], dtype=float)
    y_mm = np.asarray(source["normal_distance_mm"], dtype=float)

    outside = np.broadcast_to(np.abs(y_mm[:, None]) >= args.background_y_min_mm, intensity.shape)
    outside_values = intensity[outside & np.isfinite(intensity)]
    if args.background_intensity_i0 is None:
        i0, i0_median, i0_std = sigma_clipped_stats(outside_values, sigma=3.0, maxiters=5)
        i0 = float(i0)
        background_mode = "sigma-clipped mean recomputed from unfolded sidebands"
    else:
        i0 = float(args.background_intensity_i0)
        i0_median = None
        i0_std = None
        background_mode = "fixed external-region value; not recomputed from unfolded data"

    ratio = np.clip(intensity / i0, 0.05, None)
    tau_raw = np.maximum(-np.log(ratio), 0.0)
    dx = float(np.median(np.diff(x_mm)))
    dy = float(np.median(np.diff(y_mm)))
    tau = gaussian_filter(
        tau_raw,
        sigma=(args.smooth_mm / dy, args.smooth_mm / dx),
        mode="nearest",
    )

    tau_max = float(np.nanmax(tau))
    outside_tau = tau[outside]
    tau_bg_median = float(np.median(outside_tau))
    tau_bg_sigma = float(1.4826 * np.median(np.abs(outside_tau - tau_bg_median)))
    if args.threshold_fraction is None:
        tau_threshold = tau_bg_median + 3.0 * tau_bg_sigma
        threshold_mode = "outside-region median plus three robust standard deviations"
    else:
        tau_threshold = float(args.threshold_fraction) * tau_max
        threshold_mode = f"fixed fraction of tau_max: {args.threshold_fraction:g}"
    if tau_max <= tau_threshold:
        raise ValueError("No significant absorption above the outside-region threshold")

    q = np.clip((tau - tau_threshold) / (tau_max - tau_threshold), 0.0, 1.0)
    rho_ratio = 1.0 + (args.maximum_density_rho0 - 1.0) * q
    filament_mask = q > 0.0

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    data_path = output_dir / "euvi195_density_distribution.npz"
    np.savez_compressed(
        data_path,
        rho_over_rho0=rho_ratio,
        density_enhancement_over_rho0=rho_ratio - 1.0,
        normalized_absorption=q,
        filament_mask=filament_mask,
        optical_depth=tau,
        optical_depth_raw=tau_raw,
        intensity=intensity,
        axis_distance_mm=x_mm,
        normal_distance_mm=y_mm,
        background_intensity_i0_dn=np.asarray(i0),
        optical_depth_threshold=np.asarray(tau_threshold),
    )

    xx, yy = np.meshgrid(x_mm, y_mm)
    table = np.column_stack(
        (
            xx.ravel(),
            yy.ravel(),
            rho_ratio.ravel(),
            (rho_ratio - 1.0).ravel(),
            tau.ravel(),
            q.ravel(),
            filament_mask.astype(int).ravel(),
        )
    )
    text_path = output_dir / "euvi195_density_distribution.txt"
    np.savetxt(
        text_path,
        table,
        fmt=["%.6f", "%.6f", "%.8f", "%.8f", "%.8e", "%.8f", "%d"],
        header=(
            "x_axis_mm  y_normal_mm  rho_over_rho0  density_enhancement_over_rho0  "
            "optical_depth  normalized_absorption  filament_mask"
        ),
    )

    fig, ax = plt.subplots(figsize=(11.0, 4.2), constrained_layout=True)
    image = ax.imshow(
        rho_ratio,
        origin="lower",
        aspect="auto",
        extent=[x_mm[0], x_mm[-1], y_mm[0], y_mm[-1]],
        cmap="magma",
        vmin=1.0,
        vmax=args.maximum_density_rho0,
        interpolation="nearest",
    )
    ax.contour(x_mm, y_mm, filament_mask.astype(float), levels=[0.5], colors="cyan", linewidths=0.8)
    ax.axhline(0.0, color="white", linestyle="--", linewidth=0.9, alpha=0.7)
    ax.set_xlabel("Distance along flux-rope axis [Mm]")
    ax.set_ylabel("Normal distance [Mm]")
    ax.set_title("EUVI 195 Å absorption-constrained density distribution")
    colorbar = fig.colorbar(image, ax=ax, pad=0.015)
    colorbar.set_label(r"Density $\rho/\rho_0$")
    figure_path = output_dir / "euvi195_density_distribution.png"
    fig.savefig(figure_path, dpi=200)
    plt.close(fig)

    metadata = {
        "source": str(Path(args.input)),
        "array_shape_y_x": list(rho_ratio.shape),
        "axis_distance_mm": [float(x_mm[0]), float(x_mm[-1])],
        "normal_distance_mm": [float(y_mm[0]), float(y_mm[-1])],
        "background_mode": background_mode,
        "background_description": args.background_description,
        "background_intensity_i0_dn": i0,
        "background_region_in_unfolded_map": None if args.background_intensity_i0 is not None else f"abs(y) >= {args.background_y_min_mm:g} Mm",
        "background_sample_count": None if args.background_intensity_i0 is not None else int(outside_values.size),
        "background_arithmetic_mean_dn": None if args.background_intensity_i0 is not None else float(np.mean(outside_values)),
        "background_median_dn": None if args.background_intensity_i0 is not None else float(np.median(outside_values)),
        "background_sigma_clipped_std_dn": None if i0_std is None else float(i0_std),
        "optical_depth_formula": "tau=max(-ln(I/I0), 0), then Gaussian-smoothed",
        "smoothing_mm": args.smooth_mm,
        "outside_tau_median": tau_bg_median,
        "outside_tau_robust_sigma": tau_bg_sigma,
        "threshold_mode": threshold_mode,
        "threshold_fraction": args.threshold_fraction,
        "optical_depth_threshold": tau_threshold,
        "maximum_optical_depth": tau_max,
        "density_formula": "rho/rho0 = 1 + (rho_max/rho0 - 1) * clip((tau-tau_threshold)/(tau_max-tau_threshold), 0, 1)",
        "maximum_density_rho0": args.maximum_density_rho0,
        "filament_pixel_count": int(filament_mask.sum()),
        "outputs": {
            "npz": str(data_path),
            "text_table": str(text_path),
            "figure": str(figure_path),
        },
    }
    metadata_path = output_dir / "euvi195_density_distribution.json"
    metadata_path.write_text(json.dumps(metadata, indent=2), encoding="utf-8")

    paragraph = (
        "The EUVI 195 Å density template was derived from the axis-unfolded intensity map using the fixed "
        f"external-region reference intensity I0={i0:.6f} DN; I0 was not recomputed from any part of the "
        f"unfolded filament region. The external background definition was: {args.background_description}. "
        "The optical depth was calculated as tau=max[-ln(I/I0),0] and smoothed over "
        f"{args.smooth_mm:g} Mm. The filament boundary used {threshold_mode} "
        f"(tau_th={tau_threshold:.5f}). Inside that boundary, "
        f"the density was mapped linearly from rho0 at tau_th to {args.maximum_density_rho0:g} rho0 at "
        f"the maximum absorption (tau_max={tau_max:.5f}); outside it the density is rho0. Arrays use shape "
        f"(ny,nx)={tuple(rho_ratio.shape)}, so rho_over_rho0[iy,ix] corresponds to "
        "normal_distance_mm[iy] and axis_distance_mm[ix]. This is a relative simulation density template "
        "under the assumption of a spatially constant effective line-of-sight depth, not an absolute "
        "single-view density inversion."
    )
    (output_dir / "README_density_distribution.txt").write_text(paragraph + "\n", encoding="utf-8")
    print(json.dumps(metadata, indent=2))
    print(paragraph)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

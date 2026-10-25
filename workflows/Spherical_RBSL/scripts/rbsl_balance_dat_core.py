#!/usr/bin/env python3
"""Estimate the thin-torus RBSL balance current from an axis and AMRVAC DAT."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np


def _add_bundled_simesh() -> None:
    repository = Path(__file__).resolve().parents[3]
    source = repository / "Packages" / "simesh" / "src"
    if source.is_dir():
        sys.path.insert(0, str(source))


def read_axis(path: Path) -> np.ndarray:
    """Read ASCII or raw float64 ``lon_deg, lat_deg, r/Rsun`` axis data."""
    try:
        values = np.loadtxt(path)
        if values.ndim == 2 and values.shape[1] == 3:
            return np.asarray(values, dtype=float)
    except (UnicodeDecodeError, ValueError):
        pass
    values = np.fromfile(path, dtype="<f8")
    if values.size == 0 or values.size % 3:
        raise ValueError(f"{path} is neither 3-column ASCII nor raw float64 triples")
    return values.reshape(-1, 3)


def spherical_basis(theta: np.ndarray, phi: np.ndarray) -> tuple[np.ndarray, ...]:
    sin_theta, cos_theta = np.sin(theta), np.cos(theta)
    sin_phi, cos_phi = np.sin(phi), np.cos(phi)
    e_r = np.column_stack((sin_theta * cos_phi, sin_theta * sin_phi, cos_theta))
    e_theta = np.column_stack((cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta))
    e_phi = np.column_stack((-sin_phi, cos_phi, np.zeros_like(phi)))
    return e_r, e_theta, e_phi


def fit_axis_circle(position: np.ndarray) -> tuple[np.ndarray, float, float]:
    """Fit the axis plane and its best-fit circle; return centre, radius, RMS."""
    if len(position) < 5:
        raise ValueError("At least five axis points are required to fit a circle")
    centroid = position.mean(axis=0)
    _, singular_values, vectors = np.linalg.svd(position - centroid, full_matrices=False)
    if singular_values[1] <= np.finfo(float).eps * singular_values[0]:
        raise ValueError("Axis points are nearly collinear; cannot fit a rope circle")
    plane_basis = vectors[:2]
    projected = (position - centroid) @ plane_basis.T
    design = np.column_stack((2.0 * projected[:, 0], 2.0 * projected[:, 1], np.ones(len(position))))
    rhs = np.sum(projected**2, axis=1)
    centre_2d, _, rank, _ = np.linalg.lstsq(design, rhs, rcond=None)
    if rank < 3:
        raise ValueError("Axis circle fit is rank deficient")
    centre = centroid + centre_2d[:2] @ plane_basis
    distances = np.linalg.norm(position - centre, axis=1)
    radius = float(np.mean(distances))
    return centre, radius, float(np.sqrt(np.mean((distances - radius) ** 2)))


def resolve_magnetic_indices(dataset) -> list[int]:
    names = {
        (name.decode() if isinstance(name, bytes) else str(name)).lower(): index
        for index, name in enumerate(dataset.wnames)
    }
    aliases = (("b1", "br"), ("b2", "btheta", "bth"), ("b3", "bphi", "bph"))
    indices = []
    for candidates in aliases:
        index = next((names[name] for name in candidates if name in names), None)
        if index is None:
            raise KeyError(f"DAT is missing magnetic component; tried {candidates}, has {list(names)}")
        indices.append(index)
    return indices


def current_from_flux(flux_mx: float, minor_radius_mm: float, unit_length_mm: float,
                      unit_b_gauss: float, helicity_sign: int) -> float:
    flux_unit_mx = unit_b_gauss * (unit_length_mm * 1.0e8) ** 2
    flux_code = flux_mx / flux_unit_mx
    radius_code = minor_radius_mm / unit_length_mm
    return helicity_sign * 5.0 * np.sqrt(2.0) * flux_code / (12.0 * np.pi * radius_code)


def flux_from_current(current: np.ndarray | float, minor_radius_mm: float,
                      unit_length_mm: float, unit_b_gauss: float) -> np.ndarray:
    radius_code = minor_radius_mm / unit_length_mm
    scale = unit_b_gauss * (unit_length_mm * 1.0e8) ** 2
    scale *= 12.0 * np.pi * radius_code / (5.0 * np.sqrt(2.0))
    return np.abs(np.asarray(current)) * scale


def build_parser(mode: str) -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Infer the thin-torus RBSL balance current from an axis and a background DAT."
    )
    parser.add_argument("--dat", type=Path, required=True, help="Background AMRVAC DAT before rope insertion.")
    parser.add_argument("--path-dat", type=Path, required=True, help="RBSL axis: ASCII or raw float64 lon/lat/r triples.")
    parser.add_argument("--minor-radius-mm", type=float, required=True, help="RBSL minor radius a [Mm].")
    parser.add_argument("--solar-radius-mm", type=float, default=695.5,
                        help="Code length unit / solar radius [Mm] (default: 695.5).")
    parser.add_argument("--major-radius-mm", type=float,
                        help="Use this major radius R [Mm] instead of fitting it from the axis.")
    parser.add_argument("--unit-b-gauss", type=float,
                        help="Magnetic-field code unit [G]; enables equivalent axial-flux output.")
    parser.add_argument("--current-flux-mx", type=float,
                        help="Optional trial axial flux to compare with the inferred current.")
    parser.add_argument("--helicity", choices=("negative", "positive"), default="negative")
    parser.add_argument("--central-fraction", nargs=2, type=float, default=(0.1, 0.9), metavar=("LOW", "HIGH"),
                        help="Arc-length interval used for the global mean (default: 0.1 0.9).")
    parser.add_argument("--output-dir", type=Path, help="Default: <DAT parent>/balance_current.")
    parser.add_argument("--prefix", default="rbsl_balance")
    if mode == "stretched":
        parser.add_argument("--grid-config", type=Path, required=True,
                            help="Original amrvac.par containing stretch_dim and qstretch_baselevel.")
    return parser


def main(mode: str) -> None:
    args = build_parser(mode).parse_args()
    lower_fraction, upper_fraction = args.central_fraction
    if not 0.0 <= lower_fraction < upper_fraction <= 1.0:
        raise ValueError("--central-fraction must satisfy 0 <= LOW < HIGH <= 1")
    if args.minor_radius_mm <= 0.0 or args.solar_radius_mm <= 0.0:
        raise ValueError("Rope radius and solar radius must be positive")
    if args.current_flux_mx is not None and args.unit_b_gauss is None:
        raise ValueError("--current-flux-mx requires --unit-b-gauss")

    _add_bundled_simesh()
    from simesh.amrvac import open_dataset

    dat_path = args.dat.expanduser().resolve()
    axis_path = args.path_dat.expanduser().resolve()
    grid_config = args.grid_config.expanduser().resolve() if mode == "stretched" else None
    dataset = open_dataset(str(dat_path), ghost_width=0,
                           grid_config=str(grid_config) if grid_config else None)
    if int(dataset.ndim) != 3 or "spherical" not in str(dataset.geometry).lower():
        raise ValueError("The DAT must be a three-dimensional spherical AMRVAC snapshot")
    if mode == "stretched" and not dataset.uses_stretched_coordinates:
        raise ValueError("The supplied grid configuration did not enable stretched coordinates")
    if mode == "uniform" and dataset.uses_stretched_coordinates:
        raise ValueError("Use the stretched-grid entry point for a stretched DAT")

    axis = read_axis(axis_path)
    phi = np.radians(axis[:, 0])
    theta = 0.5 * np.pi - np.radians(axis[:, 1])
    radius_code = axis[:, 2]
    samples = np.column_stack((radius_code, theta, np.mod(phi, 2.0 * np.pi)))
    e_r, e_theta, e_phi = spherical_basis(theta, phi)
    magnetic_indices = resolve_magnetic_indices(dataset)
    dataset.load_data(field_indices=magnetic_indices)
    b_spherical = dataset.sample_points(samples, field_indices=magnetic_indices,
                                        interpolation="linear")

    position = radius_code[:, None] * args.solar_radius_mm * e_r
    circle_center, fitted_major_radius, fit_rms = fit_axis_circle(position)
    major_radius = fitted_major_radius if args.major_radius_mm is None else args.major_radius_mm
    if major_radius <= args.minor_radius_mm / 8.0:
        raise ValueError("Major radius is too small for a thin-torus estimate")
    denominator = np.log(8.0 * major_radius / args.minor_radius_mm) - 1.0
    if denominator <= 0.0:
        raise ValueError("ln(8 R/a)-1 must be positive")

    segment = np.linalg.norm(np.diff(position, axis=0), axis=1)
    arc_length = np.r_[0.0, np.cumsum(segment)]
    tangent = np.gradient(position, arc_length, axis=0, edge_order=2)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    inward = circle_center - position
    inward -= np.sum(inward * tangent, axis=1)[:, None] * tangent
    inward /= np.linalg.norm(inward, axis=1)[:, None]
    binormal = np.cross(tangent, inward)
    binormal /= np.linalg.norm(binormal, axis=1)[:, None]

    b_cartesian = (b_spherical[:, 0, None] * e_r + b_spherical[:, 1, None] * e_theta
                   + b_spherical[:, 2, None] * e_phi)
    b_strap = -np.sum(b_cartesian * binormal, axis=1)
    balance_current = (major_radius / args.solar_radius_mm) * b_strap / denominator
    fraction = arc_length / arc_length[-1]
    valid = (np.all(np.isfinite(b_spherical), axis=1) & np.isfinite(balance_current)
             & (fraction >= lower_fraction) & (fraction <= upper_fraction))
    if not np.any(valid):
        raise RuntimeError("No finite samples in the requested central axis interval")
    weights = np.gradient(arc_length)[valid]
    global_current = float(np.average(balance_current[valid], weights=weights))
    helicity_sign = -1 if args.helicity == "negative" else 1

    output_dir = (args.output_dir or dat_path.parent / "balance_current").expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    profile = {
        "index": np.arange(len(axis)), "arc_fraction": fraction, "arc_length_mm": arc_length,
        "longitude_deg": axis[:, 0], "latitude_deg": axis[:, 1], "radius_rsun": radius_code,
        "br_code": b_spherical[:, 0], "btheta_code": b_spherical[:, 1], "bphi_code": b_spherical[:, 2],
        "b_strap_code": b_strap, "i_balance_code": balance_current,
    }
    summary = {
        "dat": str(dat_path), "axis": str(axis_path), "sampling": "simesh true-coordinate DAT linear interpolation",
        "grid_mode": mode, "grid_config": str(grid_config) if grid_config else None,
        "dat_time_code": float(dataset.metadata.get("time", 0.0)), "dat_iteration": int(dataset.metadata.get("it", 0)),
        "amr_leaf_blocks": int(dataset.nleafs), "amr_max_level": int(dataset.levmax),
        "axis_points": int(len(axis)), "central_fraction": [lower_fraction, upper_fraction],
        "valid_central_points": int(np.count_nonzero(valid)), "solar_radius_mm": args.solar_radius_mm,
        "minor_radius_mm": args.minor_radius_mm, "major_radius_mm": major_radius,
        "major_radius_fitted_mm": fitted_major_radius, "axis_circle_fit_rms_mm": fit_rms,
        "r_over_a": major_radius / args.minor_radius_mm, "log_denominator": denominator,
        "b_strap_code_mean": float(np.average(b_strap[valid], weights=weights)),
        "i_balance_code_mean": global_current, "i_balance_code_median": float(np.median(balance_current[valid]),
        ), "i_balance_code_p16": float(np.percentile(balance_current[valid], 16.0)),
        "i_balance_code_p84": float(np.percentile(balance_current[valid], 84.0)),
        "matching_helicity_fraction": float(np.mean(np.sign(balance_current[valid]) == helicity_sign)),
    }
    if args.unit_b_gauss is not None:
        profile["br_gauss"] = b_spherical[:, 0] * args.unit_b_gauss
        profile["btheta_gauss"] = b_spherical[:, 1] * args.unit_b_gauss
        profile["bphi_gauss"] = b_spherical[:, 2] * args.unit_b_gauss
        profile["b_strap_gauss"] = b_strap * args.unit_b_gauss
        profile["flux_balance_local_mx"] = flux_from_current(
            balance_current, args.minor_radius_mm, args.solar_radius_mm, args.unit_b_gauss
        )
        summary["unit_b_gauss"] = args.unit_b_gauss
        summary["flux_balance_global_mx"] = float(flux_from_current(
            global_current, args.minor_radius_mm, args.solar_radius_mm, args.unit_b_gauss
        ))
    if args.current_flux_mx is not None:
        trial_current = current_from_flux(args.current_flux_mx, args.minor_radius_mm,
                                          args.solar_radius_mm, args.unit_b_gauss, helicity_sign)
        summary.update({"current_flux_mx": args.current_flux_mx, "helicity_sign": helicity_sign,
                        "i_current_code": float(trial_current),
                        "i_current_over_balance_mean": float(trial_current / global_current),
                        "current_flux_over_balance_mean": float(
                            args.current_flux_mx / summary["flux_balance_global_mx"]),})

    columns = list(profile)
    np.savetxt(output_dir / f"{args.prefix}_profile.csv", np.column_stack([profile[key] for key in columns]),
               delimiter=",", header=",".join(columns), comments="")
    (output_dir / f"{args.prefix}_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

    import matplotlib.pyplot as plt
    panel_count = 3 if args.unit_b_gauss is not None else 2
    figure, axes = plt.subplots(panel_count, 1, figsize=(10, 3.2 * panel_count), sharex=True, constrained_layout=True)
    axes = np.atleast_1d(axes)
    axes[0].plot(fraction, profile.get("b_strap_gauss", b_strap), lw=1.2)
    axes[0].set_ylabel(r"$B_{strap}$ [G]" if args.unit_b_gauss is not None else r"$B_{strap}$ [code]")
    axes[1].plot(fraction, balance_current, lw=1.2, label=r"$I_{balance}(s)$")
    axes[1].axhline(global_current, color="black", ls="--", lw=1.0, label="central weighted mean")
    if args.current_flux_mx is not None:
        axes[1].axhline(summary["i_current_code"], color="tab:red", ls="--", lw=1.0, label="trial flux")
    axes[1].set_ylabel("current [code]")
    axes[1].legend()
    if args.unit_b_gauss is not None:
        axes[2].plot(fraction, profile["flux_balance_local_mx"] / 1.0e20, lw=1.2)
        axes[2].set_ylabel(r"balance flux [$10^{20}$ Mx]")
    axes[-1].set_xlabel("axis arc-length fraction")
    for axis_plot in axes:
        axis_plot.axvspan(lower_fraction, upper_fraction, color="0.9", zorder=-2)
        axis_plot.grid(alpha=0.25)
    figure.savefig(output_dir / f"{args.prefix}_profile.png", dpi=180)
    plt.close(figure)
    print(json.dumps(summary, indent=2))
    print(output_dir)

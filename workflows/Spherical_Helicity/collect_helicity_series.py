#!/usr/bin/env python3
"""Collect per-frame spherical-helicity JSON files into a time series."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def _reference(result: dict) -> tuple[str, dict, dict]:
    preferred = result["quality"].get("preferred_reference", "top")
    code = result["integrals_code_units"]["references"][preferred]
    physical = result["physical_units"]["references"][preferred]
    return preferred, code, physical


def _runtime_seconds(frame_dir: Path) -> float:
    path = frame_dir / "run.log"
    if not path.exists():
        return float("nan")
    match = re.search(r"^real_seconds\s+([0-9.]+)$", path.read_text(errors="ignore"), re.M)
    return float(match.group(1)) if match else float("nan")


def collect(root: Path, unit_time_s: float) -> list[dict]:
    rows: list[dict] = []
    for path in sorted(root.glob("frame_[0-9][0-9][0-9][0-9]/spherical_helicity.json")):
        result = json.loads(path.read_text())
        preferred, code, physical = _reference(result)
        quality = result["quality"]
        sol = result["solenoidal_projection"]
        rows.append(
            {
                "frame": int(path.parent.name.removeprefix("frame_")),
                "dat": result["input"]["dat"],
                "time_code": result["input"]["time_code"],
                "time_s": result["input"]["time_code"] * unit_time_s,
                "iteration": result["input"]["iteration"],
                "nleafs": result["input"]["nleafs"],
                "preferred_reference": preferred,
                "quality_status": quality["status"],
                "warnings": " | ".join(quality.get("warnings", [])),
                "H_V_code": code["H_V"],
                "H_J_code": code["H_J"],
                "H_PJ_code": code["H_PJ"],
                "H_V_mx2": physical["H_V_mx2"],
                "H_J_mx2": physical["H_J_mx2"],
                "H_PJ_mx2": physical["H_PJ_mx2"],
                "H_V_over_flux2": code["H_V_over_inner_unsigned_flux_squared"],
                "inner_unsigned_flux_mx": result["physical_units"]["inner_unsigned_flux_mx"],
                "boundary_flux_relative_imbalance": result["boundary_flux"]["relative_imbalance"],
                "solenoidal_correction_l2": sol["relative_field_correction_l2"],
                "cleaned_flux_residual_l2": sol["cleaned_flux_residual_l2"],
                "curl_A_relative_l2": code["curl_A_relative_l2_interior"],
                "curl_Ap_relative_l2": code["curl_Ap_relative_l2_interior"],
                "gauge_spread_H_V": quality["gauge_spread"]["H_V"],
                "decomposition_closure": code["decomposition_closure"],
                "runtime_s": _runtime_seconds(path.parent),
            }
        )
    return rows


def write_csv(rows: list[dict], path: Path) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_summary(rows: list[dict], path: Path, expected_frames: int | None) -> dict:
    frames = [row["frame"] for row in rows]
    statuses = {name: sum(row["quality_status"] == name for row in rows) for name in sorted({row["quality_status"] for row in rows})}
    summary = {
        "frame_count": len(rows),
        "expected_frames": expected_frames,
        "complete": expected_frames is None or frames == list(range(expected_frames)),
        "frames": frames,
        "time_s_bounds": [rows[0]["time_s"], rows[-1]["time_s"]],
        "quality_status_counts": statuses,
        "H_V_mx2_bounds": [min(row["H_V_mx2"] for row in rows), max(row["H_V_mx2"] for row in rows)],
        "H_J_mx2_bounds": [min(row["H_J_mx2"] for row in rows), max(row["H_J_mx2"] for row in rows)],
        "H_PJ_mx2_bounds": [min(row["H_PJ_mx2"] for row in rows), max(row["H_PJ_mx2"] for row in rows)],
        "max_solenoidal_correction_l2": max(row["solenoidal_correction_l2"] for row in rows),
        "max_curl_A_relative_l2": max(row["curl_A_relative_l2"] for row in rows),
        "max_gauge_spread_H_V": max(row["gauge_spread_H_V"] for row in rows),
        "total_runtime_s": float(np.nansum([row["runtime_s"] for row in rows])),
    }
    path.write_text(json.dumps(summary, indent=2) + "\n")
    return summary


def plot(rows: list[dict], path: Path) -> None:
    time_min = np.array([row["time_s"] for row in rows]) / 60.0
    scale = 1.0e42
    fig, (ax0, ax1) = plt.subplots(2, 1, figsize=(8.2, 7.2), sharex=True, constrained_layout=True)
    for key, label, color in (
        ("H_V_mx2", r"$H_V$", "black"),
        ("H_J_mx2", r"$H_J$", "tab:red"),
        ("H_PJ_mx2", r"$H_{PJ}$", "tab:blue"),
    ):
        ax0.plot(time_min, np.array([row[key] for row in rows]) / scale, marker="o", ms=3.5, lw=1.4, label=label, color=color)
    ax0.axhline(0.0, color="0.75", lw=0.8)
    ax0.set_ylabel(r"Helicity ($10^{42}$ Mx$^2$)")
    ax0.legend(frameon=False, ncol=3)
    ax0.grid(alpha=0.25)

    ax1.plot(time_min, 100.0 * np.array([row["solenoidal_correction_l2"] for row in rows]), marker="o", ms=3.5, label="Solenoidal correction")
    ax1.plot(time_min, 100.0 * np.array([row["curl_A_relative_l2"] for row in rows]), marker="s", ms=3.2, label=r"curl($A$) error")
    ax1.plot(time_min, 100.0 * np.array([row["curl_Ap_relative_l2"] for row in rows]), marker="^", ms=3.2, label=r"curl($A_p$) error")
    ax1.axhline(5.0, color="tab:red", ls="--", lw=1.0, label="5% threshold")
    ax1.set_xlabel("Simulation time (min)")
    ax1.set_ylabel("Quality metric (%)")
    ax1.grid(alpha=0.25)
    ax1.legend(frameon=False, ncol=2)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("--unit-time-s", type=float, required=True)
    parser.add_argument("--expected-frames", type=int)
    args = parser.parse_args()
    rows = collect(args.root, args.unit_time_s)
    if not rows:
        raise SystemExit(f"no frame results found under {args.root}")
    write_csv(rows, args.root / "helicity_evolution.csv")
    summary = write_summary(rows, args.root / "helicity_evolution_summary.json", args.expected_frames)
    plot(rows, args.root / "helicity_evolution.png")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

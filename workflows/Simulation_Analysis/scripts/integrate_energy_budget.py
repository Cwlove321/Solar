#!/usr/bin/env python3
"""Integrate energy and selected source terms over all spherical AMR leaf cells."""

from __future__ import annotations

import argparse
import csv
import gc
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from simulation_analysis.config import data_paths, load_config
from simulation_analysis.geometry import block_volume, padded_block_centers
from simulation_analysis.physics import derive_state, gamma_from_dataset, resolve_state_fields


GHOST_WIDTH = 2


def spherical_divergence(vector, r, theta, phi):
    vr, vtheta, vphi = vector
    rr = r[:, None, None]
    sintheta = np.sin(theta)[None, :, None]
    term_r = np.gradient(rr**2 * vr, r, axis=0, edge_order=2) / rr**2
    term_theta = np.gradient(sintheta * vtheta, theta, axis=1, edge_order=2) / (rr * sintheta)
    term_phi = np.gradient(vphi, phi, axis=2, edge_order=2) / (rr * sintheta)
    return term_r + term_theta + term_phi


def spherical_curl(vector, r, theta, phi):
    br, btheta, bphi = vector
    rr = r[:, None, None]
    sintheta = np.sin(theta)[None, :, None]
    curl_r = (
        np.gradient(sintheta * bphi, theta, axis=1, edge_order=2)
        - np.gradient(btheta, phi, axis=2, edge_order=2)
    ) / (rr * sintheta)
    curl_theta = (
        np.gradient(br, phi, axis=2, edge_order=2) / sintheta
        - np.gradient(rr * bphi, r, axis=0, edge_order=2)
    ) / rr
    curl_phi = (
        np.gradient(rr * btheta, r, axis=0, edge_order=2)
        - np.gradient(br, theta, axis=1, edge_order=2)
    ) / rr
    return curl_r, curl_theta, curl_phi


def integrate_file(path: Path, config: dict) -> dict:
    from simesh.amrvac import open_dataset

    dataset = open_dataset(str(path), ghost_width=GHOST_WIDTH)
    resolved = resolve_state_fields(dataset, config)
    dataset.load_data(field_indices=[resolved[key][1] for key in resolved])
    gamma = gamma_from_dataset(dataset, config)
    padded = dataset.blocks(include_ghosts=True)
    interior = tuple(slice(GHOST_WIDTH, -GHOST_WIDTH) for _ in range(3))
    levels = np.asarray(dataset.tree_info[0], dtype=int)
    eta = float(config["physics"]["eta_code"])
    totals = {
        "volume_code": 0.0,
        "internal_energy_code": 0.0,
        "magnetic_energy_code": 0.0,
        "kinetic_energy_code": 0.0,
        "joule_power_code": 0.0,
        "compression_power_code": 0.0,
        "compression_heating_positive_code": 0.0,
        "expansion_cooling_negative_code": 0.0,
    }

    for block in range(int(dataset.nleafs)):
        level = int(levels[block])
        dv = block_volume(dataset, block, level, config)
        r, theta, phi = padded_block_centers(dataset, block, level, GHOST_WIDTH, config)
        raw = {
            logical: padded[block, position]
            for position, logical in enumerate(resolved)
        }
        state = derive_state(raw, gamma, config)
        velocity = tuple(state[name] for name in ("v1_code", "v2_code", "v3_code"))
        magnetic = tuple(raw[name] for name in ("b1", "b2", "b3"))
        divv = spherical_divergence(velocity, r, theta, phi)[interior]
        current = spherical_curl(magnetic, r, theta, phi)
        j2 = sum(component[interior] ** 2 for component in current)
        rho = raw["rho"][interior]
        pressure = state["pressure_code"][interior]
        internal = state["internal_energy_code"][interior]
        b2 = sum(component[interior] ** 2 for component in magnetic)
        v2 = sum(component[interior] ** 2 for component in velocity)
        compression = -pressure * divv

        totals["volume_code"] += float(np.sum(dv))
        totals["internal_energy_code"] += float(np.sum(internal * dv))
        totals["magnetic_energy_code"] += float(np.sum(0.5 * b2 * dv))
        totals["kinetic_energy_code"] += float(np.sum(0.5 * rho * v2 * dv))
        totals["joule_power_code"] += float(np.sum(eta * j2 * dv))
        totals["compression_power_code"] += float(np.sum(compression * dv))
        totals["compression_heating_positive_code"] += float(np.sum(np.maximum(compression, 0.0) * dv))
        totals["expansion_cooling_negative_code"] += float(np.sum(np.minimum(compression, 0.0) * dv))

    time_code = float(dataset.metadata.get("time", 0.0))
    result = {
        "file": path.name,
        "time_code": time_code,
        "time_s": time_code * float(config["units"]["time_s"]),
        "nleafs": int(dataset.nleafs),
        **totals,
    }
    energy_unit = float(config["units"]["pressure_erg_cm3"]) * float(config["units"]["length_cm"]) ** 3
    power_unit = energy_unit / float(config["units"]["time_s"])
    for name in ("internal_energy", "magnetic_energy", "kinetic_energy"):
        result[f"{name}_erg"] = result[f"{name}_code"] * energy_unit
    for name in (
        "joule_power", "compression_power", "compression_heating_positive",
        "expansion_cooling_negative"
    ):
        result[f"{name}_erg_s"] = result[f"{name}_code"] * power_unit
    del dataset, padded
    gc.collect()
    return result


def add_cumulative(rows: list[dict], config: dict) -> None:
    initial = float(rows[0]["internal_energy_code"])
    cumulative_joule = 0.0
    cumulative_compression = 0.0
    for index, row in enumerate(rows):
        if index:
            previous = rows[index - 1]
            dt = float(row["time_code"]) - float(previous["time_code"])
            cumulative_joule += 0.5 * dt * (
                float(previous["joule_power_code"]) + float(row["joule_power_code"])
            )
            cumulative_compression += 0.5 * dt * (
                float(previous["compression_power_code"]) + float(row["compression_power_code"])
            )
        delta = float(row["internal_energy_code"]) - initial
        row["cumulative_joule_energy_code"] = cumulative_joule
        row["cumulative_compression_energy_code"] = cumulative_compression
        row["internal_energy_change_code"] = delta
        row["unresolved_energy_change_code"] = delta - cumulative_joule - cumulative_compression

    energy_unit = float(config["units"]["pressure_erg_cm3"]) * float(config["units"]["length_cm"]) ** 3
    for row in rows:
        for name in (
            "cumulative_joule_energy", "cumulative_compression_energy",
            "internal_energy_change", "unresolved_energy_change"
        ):
            row[f"{name}_erg"] = float(row[f"{name}_code"]) * energy_unit


def plot_budget(rows: list[dict], config: dict, output: Path) -> None:
    time = np.asarray([float(row["time_code"]) for row in rows])
    figure, axes = plt.subplots(3, 1, figsize=(11, 12), constrained_layout=True, sharex=True)
    for label, key in (
        ("Internal", "internal_energy_code"),
        ("Magnetic", "magnetic_energy_code"),
        ("Kinetic", "kinetic_energy_code"),
    ):
        axes[0].plot(time, [row[key] for row in rows], marker="o", label=label)
    axes[0].set_yscale("log")
    axes[0].set_ylabel("Integrated energy [code units]")
    axes[0].legend(ncol=3)

    for label, key in (
        (r"$\Delta E_{int}$", "internal_energy_change_code"),
        (r"$\int dt\int \eta J^2\,dV$", "cumulative_joule_energy_code"),
        (r"$\int dt\int -p\nabla\cdot v\,dV$", "cumulative_compression_energy_code"),
        ("Unresolved/boundary", "unresolved_energy_change_code"),
    ):
        axes[1].plot(time, [row[key] for row in rows], marker="o", label=label)
    axes[1].axhline(0.0, color="0.4", linewidth=0.8)
    axes[1].set_ylabel("Cumulative energy [code units]")
    axes[1].legend(ncol=2)

    for label, key in (
        (r"Joule $\int \eta J^2\,dV$", "joule_power_code"),
        (r"Compression $\int -p\nabla\cdot v\,dV$", "compression_power_code"),
        ("Positive compression", "compression_heating_positive_code"),
    ):
        axes[2].plot(time, [row[key] for row in rows], marker="o", label=label)
    axes[2].axhline(0.0, color="0.4", linewidth=0.8)
    axes[2].set_ylabel("Integrated power [code units]")
    axes[2].set_xlabel("Simulation time [code units]")
    axes[2].legend(ncol=2)
    for axis in axes:
        axis.grid(alpha=0.22)
    top = axes[0].secondary_xaxis(
        "top",
        functions=(
            lambda value: value * float(config["units"]["time_s"]) / 60.0,
            lambda value: value * 60.0 / float(config["units"]["time_s"]),
        ),
    )
    top.set_xlabel("Physical time [min]")
    figure.savefig(output, dpi=int(config["plot"]["dpi"]), bbox_inches="tight")
    plt.close(figure)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", required=True)
    parser.add_argument("--first", type=int, default=0)
    parser.add_argument("--count", type=int, default=0)
    parser.add_argument("--stride", type=int, default=1)
    args = parser.parse_args()
    config = load_config(args.config)
    paths = data_paths(config)[args.first :: args.stride]
    if args.count > 0:
        paths = paths[: args.count]
    if not paths:
        raise SystemExit("No DAT files selected")

    rows = []
    for position, path in enumerate(paths, start=1):
        print(f"[{position}/{len(paths)}] {path.name}", flush=True)
        rows.append(integrate_file(path, config))
    rows.sort(key=lambda row: float(row["time_code"]))
    add_cumulative(rows, config)

    output = Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    csv_path = output / "integrated_energy_budget.csv"
    with csv_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    json_path = output / "integrated_energy_budget.json"
    eta_code = float(config["physics"]["eta_code"])
    eta_km2_s = (
        eta_code
        * float(config["units"]["length_cm"])
        * float(config["units"]["velocity_km_s"])
        * 1.0e5
        / 1.0e10
    )
    payload = {
        "notes": {
            "joule": "Diagnostic integral eta*J^2; it is a heating source only if the model adds it to the energy equation.",
            "compression": "Signed integral -p*div(v).",
            "unresolved": "Contains conduction, radiation, background heating, energy advection, boundary fluxes, and numerical residuals.",
            "derivatives": "Second-order spherical derivatives on each leaf block using two exchanged ghost cells.",
        },
        "eta_code": eta_code,
        "eta_km2_s": eta_km2_s,
        "rows": rows,
    }
    json_path.write_text(json.dumps(payload, indent=2) + "\n")
    plot_path = output / "integrated_energy_budget.png"
    plot_budget(rows, config, plot_path)
    print(json.dumps(rows[-1], indent=2), flush=True)
    print(plot_path, flush=True)


if __name__ == "__main__":
    main()

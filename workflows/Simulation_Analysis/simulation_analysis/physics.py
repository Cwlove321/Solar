"""State-variable interpretation and common MHD diagnostics."""

from __future__ import annotations

from typing import Mapping

import numpy as np


STATE_KEYS = ("rho", "m1", "m2", "m3", "energy", "b1", "b2", "b3")


def gamma_from_dataset(dataset, config: dict) -> float:
    configured = config["physics"].get("gamma")
    if configured is not None:
        return float(configured)
    names = [str(value).lower() for value in dataset.metadata.get("param_names", [])]
    if "gamma" not in names:
        raise KeyError("gamma is absent from DAT metadata and physics.gamma is null")
    return float(dataset.metadata["params"][names.index("gamma")])


def resolve_state_fields(dataset, config: dict) -> dict[str, tuple[str, int]]:
    """Resolve logical state names to actual DAT names and original indices."""
    names = [str(value) for value in dataset.wnames]
    lower = {name.lower(): (name, index) for index, name in enumerate(names)}
    resolved = {}
    for logical in STATE_KEYS:
        candidates = config["fields"][logical]
        if isinstance(candidates, str):
            candidates = [candidates]
        for candidate in candidates:
            if candidate.lower() in lower:
                resolved[logical] = lower[candidate.lower()]
                break
        else:
            raise KeyError(f"Cannot resolve {logical}; tried {candidates}, DAT fields={names}")
    return resolved


def velocity_from_momentum(
    rho: np.ndarray,
    momentum: tuple[np.ndarray, np.ndarray, np.ndarray],
    magnetic: tuple[np.ndarray, np.ndarray, np.ndarray],
    config: dict,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    floor = 1.0e-30
    if config["physics"]["momentum_model"] == "classical":
        return tuple(component / np.maximum(rho, floor) for component in momentum)

    b2 = sum(component**2 for component in magnetic)
    mdotb = sum(m * b for m, b in zip(momentum, magnetic))
    coefficient = np.divide(mdotb, b2, out=np.zeros_like(mdotb), where=b2 > floor)
    parallel = tuple(coefficient * b for b in magnetic)
    velocity_unit = float(config["units"]["velocity_km_s"])
    reduced_c = float(config["physics"]["reduced_c_km_s"]) / velocity_unit
    perpendicular_inertia = rho + b2 / reduced_c**2
    return tuple(
        m_parallel / np.maximum(rho, floor)
        + (m - m_parallel) / np.maximum(perpendicular_inertia, floor)
        for m, m_parallel in zip(momentum, parallel)
    )


def derive_state(raw: Mapping[str, np.ndarray], gamma: float, config: dict) -> dict[str, np.ndarray]:
    """Derive physical fields from normalized AMRVAC state arrays."""
    rho = np.asarray(raw["rho"])
    momentum = tuple(np.asarray(raw[name]) for name in ("m1", "m2", "m3"))
    magnetic = tuple(np.asarray(raw[name]) for name in ("b1", "b2", "b3"))
    velocity = velocity_from_momentum(rho, momentum, magnetic, config)
    b2 = sum(component**2 for component in magnetic)
    v2 = sum(component**2 for component in velocity)
    energy = np.asarray(raw["energy"])
    if config["physics"]["energy_kind"] == "internal":
        internal = energy
    else:
        internal = energy - 0.5 * rho * v2 - 0.5 * b2
    pressure = (gamma - 1.0) * internal

    density_floor = np.maximum(rho, 1.0e-30)
    velocity_unit = float(config["units"]["velocity_km_s"])
    magnetic_unit = float(config["units"]["magnetic_field_G"])
    temperature_unit = float(config["units"]["temperature_K"])
    number_density_unit = float(config["units"]["number_density_cm3"])
    sound_code = np.sqrt(np.maximum(gamma * pressure / density_floor, 0.0))
    alfven_code = np.sqrt(np.maximum(b2 / density_floor, 0.0))
    reduced_c_code = float(config["physics"]["reduced_c_km_s"]) / velocity_unit
    alfven_reduced_code = alfven_code / np.sqrt(1.0 + (alfven_code / reduced_c_code) ** 2)

    result = dict(raw)
    result.update(
        {
            "rho_code": rho,
            "number_density_cm3": rho * number_density_unit,
            "internal_energy_code": internal,
            "pressure_code": pressure,
            "temperature_K": pressure / density_floor * temperature_unit,
            "entropy_proxy": pressure / density_floor**gamma,
            "bmag_code": np.sqrt(b2),
            "bmag_G": np.sqrt(b2) * magnetic_unit,
            "br_G": magnetic[0] * magnetic_unit,
            "btheta_G": magnetic[1] * magnetic_unit,
            "bphi_G": magnetic[2] * magnetic_unit,
            "v1_code": velocity[0],
            "v2_code": velocity[1],
            "v3_code": velocity[2],
            "vr_km_s": velocity[0] * velocity_unit,
            "vtheta_km_s": velocity[1] * velocity_unit,
            "vphi_km_s": velocity[2] * velocity_unit,
            "vmag_km_s": np.sqrt(v2) * velocity_unit,
            "sound_speed_km_s": sound_code * velocity_unit,
            "alfven_speed_km_s": alfven_code * velocity_unit,
            "alfven_reduced_km_s": alfven_reduced_code * velocity_unit,
            "fast_speed_upper_km_s": np.sqrt(sound_code**2 + alfven_code**2)
            * velocity_unit,
            "beta": np.divide(2.0 * pressure, b2, out=np.full_like(b2, np.nan), where=b2 > 0.0),
            "kinetic_energy_density_code": 0.5 * rho * v2,
            "magnetic_energy_density_code": 0.5 * b2,
        }
    )
    return result

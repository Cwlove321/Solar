"""Configuration loading and validation for the analysis workflow."""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any


DEFAULT_CONFIG: dict[str, Any] = {
    "run_dir": ".",
    "data_glob": "data/*_[0-9][0-9][0-9][0-9].dat",
    "output_dir": "python/analysis",
    "grid": {
        "radial_type": "geometric",
        "qstretch_base": 1.041541581503301,
        "surface_radius_code": 1.0,
        "radius_unit_Mm": 695.5,
    },
    "physics": {
        "energy_kind": "internal",
        "momentum_model": "semirelativistic",
        "reduced_c_km_s": 6000.0,
        "eta_code": 2.4693930682273452e-6,
        "gamma": None,
    },
    "units": {
        "length_cm": 6.955e10,
        "number_density_cm3": 1.0e9,
        "temperature_K": 1.0e6,
        "velocity_km_s": 116.448846777562,
        "magnetic_field_G": 1.9977023903979652,
        "time_s": 5972.4772,
        "pressure_erg_cm3": 0.31754927,
    },
    "fields": {
        "rho": ["rho"],
        "m1": ["m1", "mom1"],
        "m2": ["m2", "mom2"],
        "m3": ["m3", "mom3"],
        "energy": ["e", "energy", "eint"],
        "b1": ["b1", "br"],
        "b2": ["b2", "btheta", "bth"],
        "b3": ["b3", "bphi", "bph"],
    },
    "slice": {
        "axis": "phi",
        "coordinate": 176.5,
        "angular_unit": "degree",
        "target_level": None,
        "max_height_Mm": 60.0,
        "angular_bounds_degree": None,
        "reference_angle_degree": None,
        "equal_physical_scale": True,
    },
    "plot": {
        "temperature_K": [8000.0, 4.0e6],
        "number_density_cm3": [2.0e9, 6.0e13],
        "velocity_km_s": [-250.0, 250.0],
        "magnetic_field_G": [-100.0, 100.0],
        "dpi": 180,
    },
    "log": {
        "pattern": "data/*_.log",
        "columns": [
            "it", "time", "dt", "rho_max", "m1_extreme", "m2_extreme",
            "m3_extreme", "e_max", "b1_extreme", "b2_extreme", "b3_extreme",
            "c1", "c2", "c3", "c4", "n1", "n2", "n3", "n4"
        ],
    },
}


def _merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    result = copy.deepcopy(base)
    for key, value in override.items():
        if isinstance(value, dict) and isinstance(result.get(key), dict):
            result[key] = _merge(result[key], value)
        else:
            result[key] = value
    return result


def load_config(path: str | Path) -> dict[str, Any]:
    """Load a JSON config, merge defaults, and resolve run-relative paths."""
    source = Path(path).expanduser().resolve()
    user = json.loads(source.read_text(encoding="utf-8"))
    config = _merge(DEFAULT_CONFIG, user)

    run_dir = Path(config["run_dir"]).expanduser()
    if not run_dir.is_absolute():
        run_dir = (source.parent / run_dir).resolve()
    config["run_dir"] = str(run_dir)

    output_dir = Path(config["output_dir"]).expanduser()
    if not output_dir.is_absolute():
        output_dir = run_dir / output_dir
    config["output_dir"] = str(output_dir.resolve())
    config["config_path"] = str(source)

    radial_type = config["grid"]["radial_type"]
    if radial_type not in {"uniform", "geometric"}:
        raise ValueError("grid.radial_type must be 'uniform' or 'geometric'")
    if config["physics"]["energy_kind"] not in {"internal", "total"}:
        raise ValueError("physics.energy_kind must be 'internal' or 'total'")
    if config["physics"]["momentum_model"] not in {"classical", "semirelativistic"}:
        raise ValueError(
            "physics.momentum_model must be 'classical' or 'semirelativistic'"
        )
    return config


def data_paths(config: dict[str, Any]) -> list[Path]:
    """Return sorted DAT paths selected by the run-relative glob."""
    return sorted(Path(config["run_dir"]).glob(config["data_glob"]))

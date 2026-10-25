"""Shared atmosphere schema for case generation and static validation.

JSON keys use physical units; FORTRAN_KEYS is the explicit mapping to the
``&atmosphere_list`` namelist. Defaults reproduce the corrected RBSL15
analytic atmosphere normalized by n_H at 10 Mm. A pressure anchor at h=0
requires both normalization='pressure' and reference_height_mm=0 explicitly.
"""

from __future__ import annotations

import math
from typing import Any, Mapping


DEFAULT_ATMOSPHERE: dict[str, Any] = {
    "model": "analytic_tr",
    "normalization": "number_density",
    "reference_height_mm": 10.0,
    "reference_nh_cm3": 6.2e8,
    "reference_pressure_cgs": 138.365,
    "chromosphere_k": 8000.0,
    "tanh_top_k": 1.5e6,
    "transition_height_mm": 2.0,
    "transition_width_mm": 0.25,
    "transition_k": 1.6e5,
    "conductive_flux_cgs": 2.0e5,
    "floor_height_mm": -0.2,
    "isothermal_k": 1.0e6,
}

FORTRAN_KEYS = {
    "model": "atmosphere_model",
    "normalization": "atmosphere_normalization",
    "reference_height_mm": "atmosphere_reference_height_Mm",
    "reference_nh_cm3": "atmosphere_reference_nH_cgs",
    "reference_pressure_cgs": "atmosphere_reference_pressure_cgs",
    "chromosphere_k": "atmosphere_chromosphere_K",
    "tanh_top_k": "atmosphere_tanh_top_K",
    "transition_height_mm": "atmosphere_transition_height_Mm",
    "transition_width_mm": "atmosphere_transition_width_Mm",
    "transition_k": "atmosphere_transition_K",
    "conductive_flux_cgs": "atmosphere_conductive_flux_cgs",
    "floor_height_mm": "atmosphere_floor_height_Mm",
    "isothermal_k": "atmosphere_isothermal_K",
}
_MISSING = object()


def normalize_atmosphere(overrides: object = _MISSING) -> dict[str, Any]:
    """Reject invalid values before writing files; fill omitted fields safely."""
    if overrides is _MISSING:
        overrides = {}
    if not isinstance(overrides, dict):
        raise ValueError("atmosphere must be a JSON object")
    unknown = set(overrides) - set(DEFAULT_ATMOSPHERE)
    if unknown:
        raise ValueError(f"unknown atmosphere keys: {sorted(unknown)}")
    result = {**DEFAULT_ATMOSPHERE, **overrides}
    for key, allowed in (
        ("model", {"analytic_tr", "isothermal"}),
        ("normalization", {"number_density", "pressure"}),
    ):
        value = result[key]
        if not isinstance(value, str) or value.strip().lower() not in allowed:
            raise ValueError(f"atmosphere.{key} must be one of {sorted(allowed)}")
        result[key] = value.strip().lower()
    if result["normalization"] == "pressure":
        required = {"reference_height_mm", "reference_pressure_cgs"}
        if not required.issubset(overrides):
            raise ValueError("pressure normalization requires explicit reference_height_mm and reference_pressure_cgs")
    for key in DEFAULT_ATMOSPHERE:
        if key in {"model", "normalization"}:
            continue
        value = result[key]
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"atmosphere.{key} must be a finite number")
        try:
            value = float(value)
        except (OverflowError, ValueError) as exc:
            raise ValueError(f"atmosphere.{key} must be a finite number") from exc
        if not math.isfinite(value):
            raise ValueError(f"atmosphere.{key} must be a finite number")
        result[key] = value
        if key not in {"reference_height_mm", "floor_height_mm", "transition_height_mm", "conductive_flux_cgs"} and value <= 0.0:
            raise ValueError(f"atmosphere.{key} must be positive")
    if result["conductive_flux_cgs"] < 0.0:
        raise ValueError("atmosphere.conductive_flux_cgs must be non-negative")
    if result["model"] == "analytic_tr":
        if not result["chromosphere_k"] < result["transition_k"] < result["tanh_top_k"]:
            raise ValueError("analytic_tr requires chromosphere_k < transition_k < tanh_top_k")
        if result["floor_height_mm"] >= result["transition_height_mm"]:
            raise ValueError("analytic_tr requires floor_height_mm < transition_height_mm")
    return result


def atmosphere_from_namelist(values: Mapping[str, str]) -> dict[str, Any]:
    """Read complete parsed assignments; compare numeric values, not formatting."""
    prefix = "atmosphere_list."
    reverse = {fortran.lower(): key for key, fortran in FORTRAN_KEYS.items()}
    actual = {key[len(prefix):]: value for key, value in values.items() if key.startswith(prefix)}
    missing, unknown = set(reverse) - set(actual), set(actual) - set(reverse)
    if missing or unknown:
        raise ValueError(f"atmosphere_list missing={sorted(missing)}, unknown={sorted(unknown)}")
    parsed: dict[str, Any] = {}
    for name, value in actual.items():
        key = reverse[name]
        if key in {"model", "normalization"}:
            if len(value) < 2 or value[0] not in "\"'" or value[-1] != value[0]:
                raise ValueError(f"{name} must be a quoted Fortran string")
            parsed[key] = value[1:-1]
        else:
            try:
                parsed[key] = float(value.lower().replace("d", "e"))
            except ValueError as exc:
                raise ValueError(f"{name} must be a finite scalar number") from exc
    return normalize_atmosphere(parsed)

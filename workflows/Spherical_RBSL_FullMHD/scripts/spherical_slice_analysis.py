#!/usr/bin/env python3
"""Unified spherical-AMR slice extraction, diagnostics, and plotting.

The primary slice modes are fixed ``phi``, fixed ``theta``, and fixed ``r``.
An arbitrary Cartesian plane is also supported.  All modes sample the original
MPI-AMRVAC DAT through simesh and write one reusable, uniform NPZ cache.

Derivative quantities (curl B, div B, div v, eta J^2, and J/B) can be
computed either by AMR-aware multi-scale local fits around a single central
surface or by the retained uniform three-plane slab backend.  Neither method
treats a single two-dimensional image as a complete 3-D field.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import math
import re
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable

import numpy as np


VERSION = "1.2.0"

AMR_LEVEL_COLORS = {
    1: "#355070",
    2: "#6B7D3E",
    3: "#E28B36",
    4: "#B84A62",
}

DEFAULT_UNITS = {
    "unit_length_cm": 6.955e10,
    "unit_numberdensity_cm3": 1.0e9,
    "unit_temperature_K": 1.0e6,
    "unit_density_g_cm3": 2.341668e-15,
    "unit_pressure_cgs": 0.317538,
    "unit_magneticfield_G": 1.99757357615242,
    "unit_velocity_cm_s": 1.16448846777562e7,
}
DEFAULT_UNITS["unit_time_s"] = (
    DEFAULT_UNITS["unit_length_cm"] / DEFAULT_UNITS["unit_velocity_cm_s"]
)

HE_ABUNDANCE = 0.1
LOG10_T_JCCORONA = np.array(
    [
        4.00000, 4.14230, 4.21995, 4.29761, 4.37528, 4.45294,
        4.53061, 4.60827, 4.68593, 4.76359, 4.79705, 4.83049,
        4.86394, 4.89739, 4.93084, 4.96428, 4.99773, 5.03117,
        5.06461, 5.17574, 5.28684, 5.39796, 5.50907, 5.62018,
        5.73129, 5.84240, 5.95351, 6.06461, 6.17574, 6.28684,
        6.39796, 6.50907, 6.62018, 6.73129, 6.84240, 6.95351,
        7.06461, 7.17574, 7.28684, 7.39796, 7.50907, 7.62018,
        7.73129, 7.84240, 7.95351,
    ],
    dtype=float,
)
LOG10_LAMBDA_JCCORONA = np.array(
    [
        -200.18883, -100.78630, -30.60384, -22.68481, -21.76445,
        -21.67936, -21.54218, -21.37958, -21.25172, -21.17584,
        -21.15783, -21.14491, -21.13527, -21.12837, -21.12485,
        -21.12439, -21.12642, -21.12802, -21.12548, -21.08965,
        -21.08812, -21.19542, -21.34582, -21.34839, -21.31701,
        -21.29072, -21.28900, -21.34104, -21.43122, -21.62448,
        -21.86694, -22.02897, -22.08051, -22.06057, -22.01973,
        -22.00000, -22.05161, -22.22175, -22.41452, -22.52581,
        -22.56914, -22.57486, -22.56151, -22.53969, -22.51490,
    ],
    dtype=float,
)


@dataclass
class SampleGrid:
    mode: str
    points_spherical: np.ndarray
    u: np.ndarray
    v: np.ndarray
    w: np.ndarray
    r: np.ndarray
    theta: np.ndarray
    phi: np.ndarray
    axis_map: dict[str, int] | None
    basis_u: np.ndarray | None
    basis_v: np.ndarray | None
    basis_w: np.ndarray | None
    metadata: dict[str, Any]


def _strip_fortran_comment(line: str) -> str:
    return line.split("!", 1)[0].strip()


def _parse_fortran_scalar(text: str) -> Any:
    value = text.strip().rstrip(",")
    lower = value.lower()
    if lower in {".true.", "true"}:
        return True
    if lower in {".false.", "false"}:
        return False
    if (value.startswith("'") and value.endswith("'")) or (
        value.startswith('"') and value.endswith('"')
    ):
        return value[1:-1]
    number = re.sub(r"[dD]", "e", value)
    try:
        if re.fullmatch(r"[+-]?\d+", number):
            return int(number)
        return float(number)
    except ValueError:
        return value


def parse_amrvac_par(path: Path) -> dict[str, Any]:
    """Parse scalar assignments used by the slice physics.

    This is intentionally not a complete Fortran namelist parser.  It records
    scalar values by their lower-case key and keeps namelist-qualified copies
    to make duplicated keys auditable.
    """

    result: dict[str, Any] = {}
    section = ""
    assignment = re.compile(r"^\s*([A-Za-z_]\w*(?:\(\d+\))?)\s*=\s*(.*?)\s*$")
    for raw_line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = _strip_fortran_comment(raw_line)
        if not line:
            continue
        if line.startswith("&"):
            section = line[1:].strip().lower()
            continue
        if line == "/":
            section = ""
            continue
        match = assignment.match(line)
        if not match:
            continue
        key = match.group(1).lower()
        value = _parse_fortran_scalar(match.group(2))
        result[key] = value
        if section:
            result[f"{section}.{key}"] = value
    return result


def _as_float(value: Any, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return float(default)


def _as_bool(value: Any, default: bool) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        if value.lower() in {".true.", "true"}:
            return True
        if value.lower() in {".false.", "false"}:
            return False
    return bool(default)


def _parse_vector(text: str, name: str) -> np.ndarray:
    try:
        values = np.array([float(item) for item in text.split(",")], dtype=float)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"{name} must be x,y,z") from exc
    if values.shape != (3,):
        raise argparse.ArgumentTypeError(f"{name} must contain exactly three values")
    return values


def _parse_float_list(text: str, name: str) -> np.ndarray:
    try:
        values = np.array([float(item) for item in text.split(",")], dtype=float)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"{name} must be a comma-separated list of numbers"
        ) from exc
    if values.size < 2 or np.any(~np.isfinite(values)) or np.any(values <= 0.0):
        raise argparse.ArgumentTypeError(
            f"{name} requires at least two finite positive values"
        )
    if np.any(np.diff(values) <= 0.0):
        raise argparse.ArgumentTypeError(f"{name} values must increase")
    return values


def _normalize(vector: np.ndarray, name: str) -> np.ndarray:
    norm = float(np.linalg.norm(vector))
    if not np.isfinite(norm) or norm <= 0.0:
        raise ValueError(f"{name} must be a finite non-zero vector")
    return vector / norm


def spherical_to_cartesian_points(
    radius: np.ndarray, theta: np.ndarray, phi: np.ndarray
) -> np.ndarray:
    sin_theta = np.sin(theta)
    return np.stack(
        (
            radius * sin_theta * np.cos(phi),
            radius * sin_theta * np.sin(phi),
            radius * np.cos(theta),
        ),
        axis=-1,
    )


def cartesian_to_spherical_points(points: np.ndarray) -> tuple[np.ndarray, ...]:
    x, y, z = np.moveaxis(points, -1, 0)
    radius = np.sqrt(x * x + y * y + z * z)
    with np.errstate(invalid="ignore", divide="ignore"):
        theta = np.arccos(np.clip(z / radius, -1.0, 1.0))
    phi = np.mod(np.arctan2(y, x), 2.0 * np.pi)
    return radius, theta, phi


def spherical_vectors_to_cartesian(
    radial: np.ndarray,
    polar: np.ndarray,
    azimuthal: np.ndarray,
    theta: np.ndarray,
    phi: np.ndarray,
) -> np.ndarray:
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    sin_phi = np.sin(phi)
    cos_phi = np.cos(phi)
    x = (
        radial * sin_theta * cos_phi
        + polar * cos_theta * cos_phi
        - azimuthal * sin_phi
    )
    y = (
        radial * sin_theta * sin_phi
        + polar * cos_theta * sin_phi
        + azimuthal * cos_phi
    )
    z = radial * cos_theta - polar * sin_theta
    return np.stack((x, y, z), axis=-1)


def _domain_bounds(dataset: Any) -> tuple[np.ndarray, np.ndarray]:
    bounds = np.asarray(dataset.physical_domain, dtype=float)
    return bounds[0], bounds[1]


def _interior_limits(low: float, high: float) -> tuple[float, float]:
    epsilon = max(abs(high - low) * 1.0e-8, 1.0e-10)
    return low + epsilon, high - epsilon


def _odd_slab_points(value: int) -> int:
    if value < 1 or value % 2 != 1:
        raise ValueError("--slab-points must be a positive odd integer")
    return value


def build_sample_grid(args: argparse.Namespace, dataset: Any) -> SampleGrid:
    """Construct a tensor-product sampling surface or thin slab."""

    mode = args.mode
    slab_points = _odd_slab_points(args.slab_points)
    xmin, xmax = _domain_bounds(dataset)
    rsun_mm = args.unit_length_cm / 1.0e8

    if mode in {"phi", "theta", "r"}:
        if args.value is None:
            raise ValueError(f"--value is required for mode={mode}")

        r_low, r_high = _interior_limits(float(xmin[0]), float(xmax[0]))
        theta_low, theta_high = float(xmin[1]), float(xmax[1])
        phi_low, phi_high = float(xmin[2]), float(xmax[2])
        # Match the established workflow27 convention exactly.  The 1e-5 deg
        # inset avoids asking a cell-centred interpolator to extrapolate at a
        # physical angular boundary.
        lat_low = 90.0 - math.degrees(theta_high) + 1.0e-5
        lat_high = 90.0 - math.degrees(theta_low) - 1.0e-5
        lon_low = math.degrees(phi_low) + 1.0e-5
        lon_high = math.degrees(phi_high) - 1.0e-5
        height_low = max(0.05, (r_low - 1.0) * rsun_mm)
        height_high = (r_high - 1.0) * rsun_mm
        if args.max_height_mm is not None:
            height_high = min(height_high, args.max_height_mm)

        if mode == "phi":
            u = np.linspace(
                lat_low if args.u_min is None else args.u_min,
                lat_high if args.u_max is None else args.u_max,
                args.nu,
            )
            v = np.linspace(
                height_low if args.v_min is None else args.v_min,
                height_high if args.v_max is None else args.v_max,
                args.nv,
            )
            center = math.radians(args.value)
            half = math.radians(
                0.1 if args.slab_half_width is None else args.slab_half_width
            )
            w = (
                np.array([center])
                if slab_points == 1
                else np.linspace(center - half, center + half, slab_points)
            )
            ww, vv, uu = np.meshgrid(w, v, u, indexing="ij")
            radius = 1.0 + vv / rsun_mm
            theta = 0.5 * np.pi - np.radians(uu)
            phi = ww
            axis_map = {"phi": 0, "r": 1, "theta": 2}
            xlabel, ylabel = "Latitude [deg]", "Height [Mm]"
            value_description = f"phi={args.value:g} deg"
            slab_unit = "deg"
        elif mode == "theta":
            u = np.linspace(
                lon_low if args.u_min is None else args.u_min,
                lon_high if args.u_max is None else args.u_max,
                args.nu,
            )
            v = np.linspace(
                height_low if args.v_min is None else args.v_min,
                height_high if args.v_max is None else args.v_max,
                args.nv,
            )
            center = math.radians(args.value)
            half = math.radians(
                0.1 if args.slab_half_width is None else args.slab_half_width
            )
            w = (
                np.array([center])
                if slab_points == 1
                else np.linspace(center - half, center + half, slab_points)
            )
            ww, vv, uu = np.meshgrid(w, v, u, indexing="ij")
            radius = 1.0 + vv / rsun_mm
            theta = ww
            phi = np.radians(uu)
            axis_map = {"theta": 0, "r": 1, "phi": 2}
            xlabel, ylabel = "Longitude [deg]", "Height [Mm]"
            value_description = f"theta={args.value:g} deg"
            slab_unit = "deg"
        else:
            u = np.linspace(
                lon_low if args.u_min is None else args.u_min,
                lon_high if args.u_max is None else args.u_max,
                args.nu,
            )
            v = np.linspace(
                lat_low if args.v_min is None else args.v_min,
                lat_high if args.v_max is None else args.v_max,
                args.nv,
            )
            center = float(args.value)
            half = (
                0.1 if args.slab_half_width is None else args.slab_half_width
            ) / rsun_mm
            w = (
                np.array([center])
                if slab_points == 1
                else np.linspace(center - half, center + half, slab_points)
            )
            ww, vv, uu = np.meshgrid(w, v, u, indexing="ij")
            radius = ww
            theta = 0.5 * np.pi - np.radians(vv)
            phi = np.radians(uu)
            axis_map = {"r": 0, "theta": 1, "phi": 2}
            xlabel, ylabel = "Longitude [deg]", "Latitude [deg]"
            value_description = f"r={args.value:g} R_sun"
            slab_unit = "Mm"

        points = np.stack((radius, theta, phi), axis=-1)
        if (
            np.nanmin(radius) < xmin[0]
            or np.nanmax(radius) > xmax[0]
            or np.nanmin(theta) < xmin[1]
            or np.nanmax(theta) > xmax[1]
            or np.nanmin(phi) < xmin[2]
            or np.nanmax(phi) > xmax[2]
        ):
            raise ValueError(
                "Requested fixed-coordinate slab extends outside the DAT domain. "
                "Reduce --slab-half-width or adjust the slice ranges."
            )
        return SampleGrid(
            mode=mode,
            points_spherical=points,
            u=u,
            v=v,
            w=w,
            r=radius,
            theta=theta,
            phi=phi,
            axis_map=axis_map,
            basis_u=None,
            basis_v=None,
            basis_w=None,
            metadata={
                "xlabel": xlabel,
                "ylabel": ylabel,
                "value_description": value_description,
                "slab_half_width": 0.0
                if slab_points == 1
                else float(args.slab_half_width or 0.1),
                "slab_half_width_unit": slab_unit,
            },
        )

    if mode != "plane":
        raise ValueError(f"Unknown mode: {mode}")
    if args.plane_origin_spherical is not None:
        origin_spherical = _parse_vector(
            args.plane_origin_spherical, "--plane-origin-spherical"
        )
        radius0, latitude0, longitude0 = origin_spherical
        theta0 = math.radians(90.0 - latitude0)
        phi0 = math.radians(longitude0)
        origin = spherical_to_cartesian_points(
            np.asarray(radius0), np.asarray(theta0), np.asarray(phi0)
        )
    elif args.plane_origin_cartesian is not None:
        origin = _parse_vector(args.plane_origin_cartesian, "--plane-origin-cartesian")
    else:
        raise ValueError(
            "mode=plane requires --plane-origin-spherical=r,lat_deg,lon_deg "
            "or --plane-origin-cartesian=x,y,z in R_sun"
        )
    if args.plane_normal is None:
        raise ValueError("mode=plane requires --plane-normal=x,y,z")
    normal = _normalize(_parse_vector(args.plane_normal, "--plane-normal"), "normal")
    up_guess = (
        np.array([0.0, 0.0, 1.0])
        if args.plane_up is None
        else _parse_vector(args.plane_up, "--plane-up")
    )
    basis_v = up_guess - np.dot(up_guess, normal) * normal
    if np.linalg.norm(basis_v) < 1.0e-10:
        fallback = np.array([0.0, 1.0, 0.0])
        basis_v = fallback - np.dot(fallback, normal) * normal
    basis_v = _normalize(basis_v, "in-plane up")
    basis_u = _normalize(np.cross(basis_v, normal), "in-plane horizontal")
    basis_w = normal
    u = np.linspace(-args.plane_width_mm / 2.0, args.plane_width_mm / 2.0, args.nu)
    v = np.linspace(-args.plane_height_mm / 2.0, args.plane_height_mm / 2.0, args.nv)
    half = 0.1 if args.slab_half_width is None else args.slab_half_width
    w = (
        np.array([0.0])
        if slab_points == 1
        else np.linspace(-half, half, slab_points)
    )
    ww, vv, uu = np.meshgrid(w, v, u, indexing="ij")
    cartesian = (
        origin[None, None, None, :]
        + (uu / rsun_mm)[..., None] * basis_u
        + (vv / rsun_mm)[..., None] * basis_v
        + (ww / rsun_mm)[..., None] * basis_w
    )
    radius, theta, phi = cartesian_to_spherical_points(cartesian)
    points = np.stack((radius, theta, phi), axis=-1)
    return SampleGrid(
        mode=mode,
        points_spherical=points,
        u=u,
        v=v,
        w=w,
        r=radius,
        theta=theta,
        phi=phi,
        axis_map=None,
        basis_u=basis_u,
        basis_v=basis_v,
        basis_w=basis_w,
        metadata={
            "xlabel": "Plane u [Mm]",
            "ylabel": "Plane v [Mm]",
            "value_description": "arbitrary Cartesian plane",
            "origin_cartesian_Rsun": origin.tolist(),
            "basis_u": basis_u.tolist(),
            "basis_v": basis_v.tolist(),
            "basis_normal": basis_w.tolist(),
            "slab_half_width": 0.0 if slab_points == 1 else float(half),
            "slab_half_width_unit": "Mm",
        },
    )


def _dataset_parameter(dataset: Any, name: str, default: float) -> float:
    names = list(dataset.metadata.get("param_names", []))
    params = np.asarray(dataset.metadata.get("params", []))
    if name in names:
        return float(params[names.index(name)])
    return float(default)


def sample_raw_fields(
    dataset: Any,
    grid: SampleGrid,
    interpolation: str,
    chunk_size: int,
) -> tuple[dict[str, np.ndarray], dict[str, int]]:
    names = {str(name): index for index, name in enumerate(dataset.wnames)}
    required = ["rho", "e", "b1", "b2", "b3"]
    missing_required = [name for name in required if name not in names]
    if missing_required:
        raise KeyError(f"DAT is missing required fields: {missing_required}")
    momentum = [name for name in ("m1", "m2", "m3") if name in names]
    primitive_velocity = [name for name in ("v1", "v2", "v3") if name in names]
    requested_names = required + momentum + primitive_velocity
    requested_indices = [names[name] for name in requested_names]
    dataset.load_data(field_indices=requested_indices)

    points = np.ascontiguousarray(grid.points_spherical.reshape(-1, 3))
    values = np.empty((points.shape[0], len(requested_names)), dtype=float)
    linear_count = 0
    nearest_count = 0
    invalid_count = 0
    for start in range(0, points.shape[0], chunk_size):
        stop = min(start + chunk_size, points.shape[0])
        chunk = dataset.sample_points(
            points[start:stop],
            field_indices=requested_indices,
            interpolation=interpolation,
        )
        bad = ~np.all(np.isfinite(chunk), axis=1)
        linear_count += int(np.count_nonzero(~bad))
        if np.any(bad) and interpolation != "nearest":
            nearest = dataset.sample_points(
                points[start:stop][bad],
                field_indices=requested_indices,
                interpolation="nearest",
            )
            chunk[bad] = nearest
            nearest_count += int(np.count_nonzero(np.all(np.isfinite(nearest), axis=1)))
        invalid_count += int(np.count_nonzero(~np.all(np.isfinite(chunk), axis=1)))
        values[start:stop] = chunk

    shape = grid.points_spherical.shape[:-1]
    raw = {
        name: values[:, column].reshape(shape)
        for column, name in enumerate(requested_names)
    }
    return raw, {
        "linear_or_requested_samples": linear_count,
        "nearest_fallback_samples": nearest_count,
        "invalid_samples": invalid_count,
        "total_samples": int(points.shape[0]),
    }


def reconstruct_velocity(
    raw: dict[str, np.ndarray],
    par: dict[str, Any],
    units: dict[str, float],
) -> tuple[np.ndarray, str]:
    if all(name in raw for name in ("v1", "v2", "v3")):
        return (
            np.stack((raw["v1"], raw["v2"], raw["v3"]), axis=-1),
            "primitive velocity fields v1/v2/v3",
        )
    if not all(name in raw for name in ("m1", "m2", "m3")):
        shape = raw["rho"].shape + (3,)
        return np.full(shape, np.nan), "velocity unavailable"
    momentum = np.stack((raw["m1"], raw["m2"], raw["m3"]), axis=-1)
    rho = raw["rho"]
    if not _as_bool(par.get("mhd_semirelativistic"), False):
        with np.errstate(divide="ignore", invalid="ignore"):
            return momentum / rho[..., None], "classical conserved momentum/rho"
    magnetic = np.stack((raw["b1"], raw["b2"], raw["b3"]), axis=-1)
    reduced_c = _as_float(par.get("mhd_reduced_c"), 6.0e8)
    inv_c2_code = (units["unit_velocity_cm_s"] / reduced_c) ** 2
    b2 = np.sum(magnetic * magnetic, axis=-1)
    mdotb = np.sum(momentum * magnetic, axis=-1)
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        factor = 1.0 / (rho * (rho + b2 * inv_c2_code))
        velocity = factor[..., None] * (
            momentum * rho[..., None]
            + magnetic * (mdotb * inv_c2_code)[..., None]
        )
    return velocity, "AMRVAC semi-relativistic conserved momentum inversion"


def locate_owner_leaves(dataset: Any, points: np.ndarray) -> np.ndarray:
    """Locate active AMR leaves without scanning every leaf for every query.

    The DAT tree stores each active leaf's AMR level and integer block id.
    Reconstructing the global nodes at each level lets us map a physical point
    to the corresponding block id, then resolve that id through sorted integer
    keys.  Finest active leaves are tested first.
    """

    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 3:
        raise ValueError("points must have shape (n, 3)")
    cache = getattr(dataset, "_spherical_slice_owner_lookup", None)
    if cache is None:
        coordinates = dataset.coordinates
        levels = np.asarray(coordinates.block_levels, dtype=int)
        block_ids = np.asarray(coordinates.block_ids, dtype=int)
        levmax = int(np.max(levels))
        base_blocks = np.asarray(dataset.domain_nx, dtype=int) // np.asarray(
            dataset.block_nx, dtype=int
        )
        key_bases = base_blocks * 2 ** (levmax - 1) + 1

        def encode(ids: np.ndarray) -> np.ndarray:
            return (
                (ids[:, 0] * key_bases[1] + ids[:, 1]) * key_bases[2]
                + ids[:, 2]
            )

        by_level: dict[int, dict[str, Any]] = {}
        for level in sorted(set(levels), reverse=True):
            leaves = np.flatnonzero(levels == level)
            keys = encode(block_ids[leaves])
            order = np.argsort(keys)
            by_level[int(level)] = {
                "keys": keys[order],
                "leaves": leaves[order],
                "nodes": tuple(
                    coordinates._level_nodes(int(level), dim) for dim in range(3)
                ),
            }
        cache = {
            "levels": levels,
            "key_bases": key_bases,
            "by_level": by_level,
        }
        dataset._spherical_slice_owner_lookup = cache

    low, high = _domain_bounds(dataset)
    inside = np.all(np.isfinite(points), axis=1)
    inside &= np.all(points >= low[None, :], axis=1)
    inside &= np.all(points <= high[None, :], axis=1)
    owner = np.full(len(points), -1, dtype=np.int32)
    key_bases = cache["key_bases"]
    for level, item in cache["by_level"].items():
        target = np.flatnonzero(inside & (owner < 0))
        if target.size == 0:
            break
        block_id = np.empty((target.size, 3), dtype=np.int64)
        for dim, nodes in enumerate(item["nodes"]):
            cell = np.searchsorted(nodes, points[target, dim], side="right") - 1
            cell = np.clip(cell, 0, len(nodes) - 2)
            block_id[:, dim] = cell // int(dataset.block_nx[dim])
        keys = (
            (block_id[:, 0] * key_bases[1] + block_id[:, 1]) * key_bases[2]
            + block_id[:, 2]
        )
        position = np.searchsorted(item["keys"], keys)
        valid = position < len(item["keys"])
        matched = np.zeros_like(valid)
        matched[valid] = item["keys"][position[valid]] == keys[valid]
        owner[target[matched]] = item["leaves"][position[matched]]
    return owner


def owner_cell_widths(
    dataset: Any, points: np.ndarray, owner: np.ndarray
) -> np.ndarray:
    """Return native (dr, dtheta, dphi) widths of each point's owner cell."""

    points = np.asarray(points, dtype=float)
    owner = np.asarray(owner, dtype=int)
    widths = np.full((len(points), 3), np.nan, dtype=float)
    nodes = dataset.block_node_coordinates
    for leaf in np.unique(owner[owner >= 0]):
        target = np.flatnonzero(owner == leaf)
        for dim in range(3):
            local_nodes = nodes[dim][leaf]
            cell = np.searchsorted(
                local_nodes, points[target, dim], side="right"
            ) - 1
            cell = np.clip(cell, 0, len(local_nodes) - 2)
            widths[target, dim] = np.diff(local_nodes)[cell]
    return widths


def sample_points_with_known_owner(
    dataset: Any,
    points: np.ndarray,
    owner: np.ndarray,
    field_indices: list[int],
    interpolation: str = "linear",
) -> np.ndarray:
    """Use simesh interpolation after resolving owner leaves once."""

    points = np.asarray(points, dtype=float)
    owner = np.asarray(owner, dtype=int)
    loaded_map = dataset._loaded_field_map()
    missing = [int(index) for index in field_indices if int(index) not in loaded_map]
    if missing:
        raise ValueError(f"Requested field indices are not loaded: {missing}")
    columns = [loaded_map[int(index)] for index in field_indices]
    result = np.full((len(points), len(columns)), np.nan, dtype=float)
    centers = dataset.block_center_coordinates
    for leaf in np.unique(owner[owner >= 0]):
        target = np.flatnonzero(owner == leaf)
        if interpolation == "nearest":
            cell_index = [
                np.abs(
                    centers[dim][leaf, None, :] - points[target, dim, None]
                ).argmin(axis=1)
                for dim in range(3)
            ]
            values = dataset.data[leaf, columns]
            result[target] = values[
                :, cell_index[0], cell_index[1], cell_index[2]
            ].T
            continue
        centre_inside = np.ones(len(target), dtype=bool)
        for dim in range(3):
            centre_inside &= (
                points[target, dim] >= centers[dim][leaf, 0]
            ) & (points[target, dim] <= centers[dim][leaf, -1])
        interior = target[centre_inside]
        boundary = target[~centre_inside]
        if interior.size:
            dataset._sample_leaf_centres_linear(
                int(leaf), points, interior, columns, result
            )
        if boundary.size:
            dataset._sample_leaf_boundary_linear(
                int(leaf), points, boundary, columns, result
            )
    return result


def spherical_curl_div_from_partials(
    vector: np.ndarray,
    partials: np.ndarray,
    radius: np.ndarray,
    theta: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Curl and divergence from partial derivatives of spherical components.

    ``partials[..., coordinate, component]`` uses coordinate order
    ``(r, theta, phi)`` and component order ``(r, theta, phi)``.
    """

    radial, polar, azimuthal = np.moveaxis(vector, -1, 0)
    sin_theta = np.sin(theta)
    cos_theta = np.cos(theta)
    d_dr = partials[..., 0, :]
    d_dtheta = partials[..., 1, :]
    d_dphi = partials[..., 2, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        divergence = (
            d_dr[..., 0]
            + 2.0 * radial / radius
            + d_dtheta[..., 1] / radius
            + cos_theta * polar / (radius * sin_theta)
            + d_dphi[..., 2] / (radius * sin_theta)
        )
        curl_r = (
            sin_theta * d_dtheta[..., 2]
            + cos_theta * azimuthal
            - d_dphi[..., 1]
        ) / (radius * sin_theta)
        curl_theta = (
            d_dphi[..., 0] / (radius * sin_theta)
            - d_dr[..., 2]
            - azimuthal / radius
        )
        curl_phi = (
            d_dr[..., 1]
            + polar / radius
            - d_dtheta[..., 0] / radius
        )
    return np.stack((curl_r, curl_theta, curl_phi), axis=-1), divergence


def _offset_state(
    dataset: Any,
    points: np.ndarray,
    owner: np.ndarray,
    requested_names: list[str],
    requested_indices: list[int],
    par: dict[str, Any],
    units: dict[str, float],
) -> tuple[np.ndarray, np.ndarray]:
    values = sample_points_with_known_owner(
        dataset, points, owner, requested_indices, interpolation="linear"
    )
    fields = {
        name: values[:, column]
        for column, name in enumerate(requested_names)
    }
    magnetic = np.stack(
        (fields["b1"], fields["b2"], fields["b3"]), axis=-1
    )
    velocity, _ = reconstruct_velocity(fields, par, units)
    return magnetic, velocity


def local_fit_derivatives(
    dataset: Any,
    grid: SampleGrid,
    raw: dict[str, np.ndarray],
    par: dict[str, Any],
    units: dict[str, float],
    fit_scales: np.ndarray,
    residual_max: float,
    sensitivity_max: float,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Fit 3-D derivatives around central slice points on native cell scales."""

    center_index = grid.points_spherical.shape[0] // 2
    points = np.ascontiguousarray(
        grid.points_spherical[center_index].reshape(-1, 3)
    )
    shape = grid.points_spherical.shape[1:3]
    owner = locate_owner_leaves(dataset, points)
    widths = owner_cell_widths(dataset, points, owner)
    levels_by_leaf = np.asarray(dataset.coordinates.block_levels, dtype=int)
    owner_level = np.full(len(points), -1, dtype=int)
    valid_owner = owner >= 0
    owner_level[valid_owner] = levels_by_leaf[owner[valid_owner]]

    names = {str(name): index for index, name in enumerate(dataset.wnames)}
    requested_names = [
        name
        for name in ("rho", "e", "m1", "m2", "m3", "b1", "b2", "b3", "v1", "v2", "v3")
        if name in raw
    ]
    requested_indices = [names[name] for name in requested_names]
    nscale = len(fit_scales)
    npoint = len(points)
    pair_b = np.full((nscale, npoint, 3, 3), np.nan, dtype=float)
    pair_v = np.full_like(pair_b, np.nan)
    deltas = np.full((nscale, npoint, 3), np.nan, dtype=float)
    owner_at_largest: list[np.ndarray] = []

    for coordinate in range(3):
        for scale_index, scale in enumerate(fit_scales):
            delta = widths[:, coordinate] * float(scale)
            deltas[scale_index, :, coordinate] = delta
            plus = points.copy()
            minus = points.copy()
            plus[:, coordinate] += delta
            minus[:, coordinate] -= delta
            owner_plus = locate_owner_leaves(dataset, plus)
            owner_minus = locate_owner_leaves(dataset, minus)
            if scale_index == nscale - 1:
                owner_at_largest.extend((owner_plus, owner_minus))
            magnetic_plus, velocity_plus = _offset_state(
                dataset,
                plus,
                owner_plus,
                requested_names,
                requested_indices,
                par,
                units,
            )
            magnetic_minus, velocity_minus = _offset_state(
                dataset,
                minus,
                owner_minus,
                requested_names,
                requested_indices,
                par,
                units,
            )
            with np.errstate(divide="ignore", invalid="ignore"):
                pair_b[scale_index, :, coordinate, :] = (
                    magnetic_plus - magnetic_minus
                ) / (2.0 * delta[:, None])
                pair_v[scale_index, :, coordinate, :] = (
                    velocity_plus - velocity_minus
                ) / (2.0 * delta[:, None])

    def weighted_fit(pair: np.ndarray) -> np.ndarray:
        weight = deltas[..., None] ** 2
        finite = np.isfinite(pair)
        numerator = np.nansum(np.where(finite, weight * pair, 0.0), axis=0)
        denominator = np.sum(np.where(finite, weight, 0.0), axis=0)
        with np.errstate(divide="ignore", invalid="ignore"):
            return numerator / denominator

    gradient_b = weighted_fit(pair_b)
    gradient_v = weighted_fit(pair_v)
    center_b = np.stack(
        (
            _central(raw["b1"]),
            _central(raw["b2"]),
            _central(raw["b3"]),
        ),
        axis=-1,
    ).reshape(npoint, 3)
    velocity_full, _ = reconstruct_velocity(raw, par, units)
    center_v = _central(velocity_full).reshape(npoint, 3)
    radius = points[:, 0]
    theta = points[:, 1]
    current, divb = spherical_curl_div_from_partials(
        center_b, gradient_b, radius, theta
    )
    _, divv = spherical_curl_div_from_partials(
        center_v, gradient_v, radius, theta
    )

    current_scales = np.full((nscale, npoint, 3), np.nan)
    divb_scales = np.full((nscale, npoint), np.nan)
    divv_scales = np.full((nscale, npoint), np.nan)
    for scale_index in range(nscale):
        current_scales[scale_index], divb_scales[scale_index] = (
            spherical_curl_div_from_partials(
                center_b, pair_b[scale_index], radius, theta
            )
        )
        _, divv_scales[scale_index] = spherical_curl_div_from_partials(
            center_v, pair_v[scale_index], radius, theta
        )

    def fit_residual(pair: np.ndarray, fitted: np.ndarray, amplitude: np.ndarray) -> np.ndarray:
        actual = 2.0 * deltas[..., None] * pair
        predicted = 2.0 * deltas[..., None] * fitted[None, ...]
        numerator = np.nansum((actual - predicted) ** 2, axis=(0, 2, 3))
        denominator = np.nansum(actual**2, axis=(0, 2, 3))
        denominator += (
            nscale * 3 * 3 * (1.0e-6 * np.maximum(amplitude, 1.0e-30)) ** 2
        )
        with np.errstate(divide="ignore", invalid="ignore"):
            return np.sqrt(numerator / denominator)

    bmag = np.linalg.norm(center_b, axis=-1)
    vmag = np.linalg.norm(center_v, axis=-1)
    b_residual = fit_residual(pair_b, gradient_b, bmag)
    v_residual = fit_residual(pair_v, gradient_v, vmag)
    physical_widths = np.column_stack(
        (
            widths[:, 0],
            radius * widths[:, 1],
            radius * np.sin(theta) * widths[:, 2],
        )
    )
    characteristic_width = np.cbrt(np.prod(physical_widths, axis=1))
    current_magnitude = np.linalg.norm(current, axis=-1)
    current_difference = np.linalg.norm(
        current_scales - current[None, :, :], axis=-1
    )
    current_floor = 0.05 * bmag / np.maximum(
        characteristic_width, np.finfo(float).tiny
    )
    def nanmax_without_warning(values: np.ndarray) -> np.ndarray:
        finite = np.any(np.isfinite(values), axis=0)
        maximum = np.max(np.where(np.isfinite(values), values, -np.inf), axis=0)
        return np.where(finite, maximum, np.nan)

    with np.errstate(divide="ignore", invalid="ignore"):
        b_sensitivity = nanmax_without_warning(
            current_difference
            / (current_magnitude[None, :] + current_floor[None, :])
        )
        v_scale = np.maximum(
            np.abs(divv),
            0.05 * vmag / np.maximum(characteristic_width, np.finfo(float).tiny),
        )
        v_scale = np.maximum(v_scale, 1.0e-12)
        v_sensitivity = nanmax_without_warning(
            np.abs(divv_scales - divv[None, :]) / v_scale[None, :]
        )

    cross_leaf = np.zeros(npoint, dtype=bool)
    cross_level = np.zeros(npoint, dtype=bool)
    outside_stencil = np.zeros(npoint, dtype=bool)
    for offset_owner in owner_at_largest:
        valid = offset_owner >= 0
        outside_stencil |= ~valid
        cross_leaf |= valid & (offset_owner != owner)
        offset_level = np.full(npoint, -1, dtype=int)
        offset_level[valid] = levels_by_leaf[offset_owner[valid]]
        cross_level |= valid & (offset_level != owner_level)

    magnetic_reliable = (
        valid_owner
        & np.isfinite(current_magnitude)
        & np.isfinite(divb)
        & np.isfinite(b_residual)
        & np.isfinite(b_sensitivity)
        & ~outside_stencil
        & (b_residual <= residual_max)
        & (b_sensitivity <= sensitivity_max)
    )
    magnetic_strict_reliable = magnetic_reliable & ~cross_level
    velocity_reliable = (
        valid_owner
        & np.isfinite(divv)
        & np.isfinite(v_residual)
        & np.isfinite(v_sensitivity)
        & ~outside_stencil
        & (v_residual <= residual_max)
        & (v_sensitivity <= sensitivity_max)
    )
    velocity_strict_reliable = velocity_reliable & ~cross_level
    rsun_mm = units["unit_length_cm"] / 1.0e8

    def plane(values: np.ndarray) -> np.ndarray:
        return values.reshape(shape + values.shape[1:])

    arrays = {
        "owner_leaf_id": plane(owner),
        "owner_amr_level": plane(owner_level),
        "cross_leaf_stencil": plane(cross_leaf),
        "cross_level_stencil": plane(cross_level),
        "stencil_outside_domain": plane(outside_stencil),
        "local_cell_dr_Mm": plane(widths[:, 0] * rsun_mm),
        "local_cell_r_dtheta_Mm": plane(physical_widths[:, 1] * rsun_mm),
        "local_cell_r_sintheta_dphi_Mm": plane(physical_widths[:, 2] * rsun_mm),
        "derivative_fit_residual": plane(b_residual),
        "derivative_scale_sensitivity": plane(b_sensitivity),
        "velocity_derivative_fit_residual": plane(v_residual),
        "velocity_derivative_scale_sensitivity": plane(v_sensitivity),
        "derivative_reliable_mask": plane(magnetic_reliable),
        "derivative_strict_reliable_mask": plane(magnetic_strict_reliable),
        "velocity_derivative_reliable_mask": plane(velocity_reliable),
        "velocity_derivative_strict_reliable_mask": plane(
            velocity_strict_reliable
        ),
        "j_r_code": plane(current[:, 0]),
        "j_theta_code": plane(current[:, 1]),
        "j_phi_code": plane(current[:, 2]),
        "current_magnitude_code": plane(current_magnitude),
        "divb_code": plane(divb),
        "divv_code": plane(divv),
        "divv_per_s": plane(divv / units["unit_time_s"]),
    }
    current_cartesian = spherical_vectors_to_cartesian(
        current[:, 0], current[:, 1], current[:, 2], theta, points[:, 2]
    )
    arrays["j_x_code"] = plane(current_cartesian[:, 0])
    arrays["j_y_code"] = plane(current_cartesian[:, 1])
    arrays["j_z_code"] = plane(current_cartesian[:, 2])
    metadata = {
        "backend": "local-fit",
        "fit_scales_in_owner_cells": [float(value) for value in fit_scales],
        "offset_samples": int(npoint * 6 * nscale),
        "residual_max": float(residual_max),
        "sensitivity_max": float(sensitivity_max),
        "owner_level_counts": {
            str(level): int(np.count_nonzero(owner_level == level))
            for level in np.unique(owner_level)
        },
        "cross_leaf_fraction": float(np.mean(cross_leaf)),
        "cross_level_fraction": float(np.mean(cross_level)),
        "outside_domain_stencil_fraction": float(np.mean(outside_stencil)),
        "magnetic_reliable_fraction": float(np.mean(magnetic_reliable)),
        "magnetic_strict_reliable_fraction": float(
            np.mean(magnetic_strict_reliable)
        ),
        "velocity_reliable_fraction": float(np.mean(velocity_reliable)),
        "velocity_strict_reliable_fraction": float(
            np.mean(velocity_strict_reliable)
        ),
    }
    return arrays, metadata


def _gradient(values: np.ndarray, coordinates: np.ndarray, axis: int) -> np.ndarray:
    if len(coordinates) < 3:
        raise ValueError("At least three coordinates are required for a derivative")
    return np.gradient(values, coordinates, axis=axis, edge_order=2)


def spherical_vector_derivatives(
    vector: np.ndarray,
    grid: SampleGrid,
) -> tuple[np.ndarray, np.ndarray]:
    """Return spherical curl and divergence on a tensor-product spherical slab."""

    if grid.axis_map is None:
        raise ValueError("Spherical derivatives require a fixed-coordinate grid")
    axes = grid.axis_map
    def coordinate_line(values: np.ndarray, axis: int) -> np.ndarray:
        selector = [0, 0, 0]
        selector[axis] = slice(None)
        return np.asarray(values[tuple(selector)], dtype=float)

    r_vector = coordinate_line(grid.r, axes["r"])
    theta_vector = coordinate_line(grid.theta, axes["theta"])
    phi_vector = coordinate_line(grid.phi, axes["phi"])
    radial, polar, azimuthal = np.moveaxis(vector, -1, 0)
    radius = grid.r
    sin_theta = np.sin(grid.theta)

    d_r2_br_dr = _gradient(radius * radius * radial, r_vector, axes["r"])
    d_sin_btheta_dtheta = _gradient(
        sin_theta * polar, theta_vector, axes["theta"]
    )
    d_bphi_dphi = _gradient(azimuthal, phi_vector, axes["phi"])
    with np.errstate(divide="ignore", invalid="ignore"):
        divergence = (
            d_r2_br_dr / (radius * radius)
            + d_sin_btheta_dtheta / (radius * sin_theta)
            + d_bphi_dphi / (radius * sin_theta)
        )

    d_sin_bphi_dtheta = _gradient(
        sin_theta * azimuthal, theta_vector, axes["theta"]
    )
    d_btheta_dphi = _gradient(polar, phi_vector, axes["phi"])
    d_br_dphi = _gradient(radial, phi_vector, axes["phi"])
    d_r_bphi_dr = _gradient(radius * azimuthal, r_vector, axes["r"])
    d_r_btheta_dr = _gradient(radius * polar, r_vector, axes["r"])
    d_br_dtheta = _gradient(radial, theta_vector, axes["theta"])
    with np.errstate(divide="ignore", invalid="ignore"):
        curl_r = (d_sin_bphi_dtheta - d_btheta_dphi) / (radius * sin_theta)
        curl_theta = d_br_dphi / (radius * sin_theta) - d_r_bphi_dr / radius
        curl_phi = d_r_btheta_dr / radius - d_br_dtheta / radius
    return np.stack((curl_r, curl_theta, curl_phi), axis=-1), divergence


def cartesian_plane_vector_derivatives(
    vector_cartesian: np.ndarray,
    grid: SampleGrid,
    rsun_mm: float,
) -> tuple[np.ndarray, np.ndarray]:
    if grid.basis_u is None or grid.basis_v is None or grid.basis_w is None:
        raise ValueError("Cartesian plane derivatives require plane basis vectors")
    u = grid.u / rsun_mm
    v = grid.v / rsun_mm
    w = grid.w / rsun_mm
    component_u = np.einsum("...i,i->...", vector_cartesian, grid.basis_u)
    component_v = np.einsum("...i,i->...", vector_cartesian, grid.basis_v)
    component_w = np.einsum("...i,i->...", vector_cartesian, grid.basis_w)
    d_bu_du = _gradient(component_u, u, 2)
    d_bu_dv = _gradient(component_u, v, 1)
    d_bu_dw = _gradient(component_u, w, 0)
    d_bv_du = _gradient(component_v, u, 2)
    d_bv_dv = _gradient(component_v, v, 1)
    d_bv_dw = _gradient(component_v, w, 0)
    d_bw_du = _gradient(component_w, u, 2)
    d_bw_dv = _gradient(component_w, v, 1)
    d_bw_dw = _gradient(component_w, w, 0)
    curl_basis = np.stack(
        (d_bw_dv - d_bv_dw, d_bu_dw - d_bw_du, d_bv_du - d_bu_dv),
        axis=-1,
    )
    curl_cartesian = (
        curl_basis[..., 0, None] * grid.basis_u
        + curl_basis[..., 1, None] * grid.basis_v
        + curl_basis[..., 2, None] * grid.basis_w
    )
    divergence = d_bu_du + d_bv_dv + d_bw_dw
    return curl_cartesian, divergence


def _central(values: np.ndarray) -> np.ndarray:
    return values[values.shape[0] // 2]


def derive_fields(
    raw: dict[str, np.ndarray],
    grid: SampleGrid,
    dataset: Any,
    par: dict[str, Any],
    units: dict[str, float],
    derivative_backend: str,
    fit_scales: np.ndarray,
    fit_residual_max: float,
    fit_sensitivity_max: float,
) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    gamma = _dataset_parameter(
        dataset, "gamma", _as_float(par.get("mhd_gamma"), 5.0 / 3.0)
    )
    rho = raw["rho"]
    internal_energy = raw["e"]
    magnetic_spherical = np.stack(
        (raw["b1"], raw["b2"], raw["b3"]), axis=-1
    )
    velocity_spherical, velocity_method = reconstruct_velocity(raw, par, units)
    magnetic_cartesian = spherical_vectors_to_cartesian(
        raw["b1"], raw["b2"], raw["b3"], grid.theta, grid.phi
    )
    velocity_cartesian = spherical_vectors_to_cartesian(
        velocity_spherical[..., 0],
        velocity_spherical[..., 1],
        velocity_spherical[..., 2],
        grid.theta,
        grid.phi,
    )
    central_index = rho.shape[0] // 2
    radius = grid.r[central_index]
    theta = grid.theta[central_index]
    phi = grid.phi[central_index]
    rho2 = _central(rho)
    e2 = _central(internal_energy)
    b_sph2 = _central(magnetic_spherical)
    b_cart2 = _central(magnetic_cartesian)
    v_sph2 = _central(velocity_spherical)
    v_cart2 = _central(velocity_cartesian)
    bmag = np.linalg.norm(b_sph2, axis=-1)
    speed = np.linalg.norm(v_sph2, axis=-1)
    pressure_code = (gamma - 1.0) * e2
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        temperature = pressure_code / rho2 * units["unit_temperature_K"]
        number_density = rho2 * units["unit_numberdensity_cm3"]
        mass_density = rho2 * units["unit_density_g_cm3"]
        pressure_cgs = pressure_code * units["unit_pressure_cgs"]
        entropy_proxy = pressure_code / rho2**gamma
        speed_km_s = speed * units["unit_velocity_cm_s"] / 1.0e5
        alfven_classical = (
            bmag / np.sqrt(rho2) * units["unit_velocity_cm_s"] / 1.0e5
        )
        plasma_beta = 2.0 * pressure_code / (bmag * bmag)

    reduced_c_cm_s = _as_float(par.get("mhd_reduced_c"), 6.0e8)
    with np.errstate(divide="ignore", invalid="ignore"):
        alfven_limited = 1.0 / np.sqrt(
            1.0 / np.maximum(alfven_classical, 1.0e-300) ** 2
            + 1.0 / (reduced_c_cm_s / 1.0e5) ** 2
        )

    arrays: dict[str, np.ndarray] = {
        "u": grid.u,
        "v": grid.v,
        "radius_Rsun": radius,
        "theta_rad": theta,
        "phi_rad": phi,
        "valid_mask": np.isfinite(rho2),
        "rho_code": rho2,
        "internal_energy_code": e2,
        "temperature_K": temperature,
        "number_density_cm3": number_density,
        "mass_density_g_cm3": mass_density,
        "pressure_code": pressure_code,
        "pressure_cgs": pressure_cgs,
        "entropy_proxy": entropy_proxy,
        "b_r_code": b_sph2[..., 0],
        "b_theta_code": b_sph2[..., 1],
        "b_phi_code": b_sph2[..., 2],
        "b_x_code": b_cart2[..., 0],
        "b_y_code": b_cart2[..., 1],
        "b_z_code": b_cart2[..., 2],
        "b_magnitude_code": bmag,
        "b_magnitude_G": bmag * units["unit_magneticfield_G"],
        "v_r_code": v_sph2[..., 0],
        "v_theta_code": v_sph2[..., 1],
        "v_phi_code": v_sph2[..., 2],
        "v_x_code": v_cart2[..., 0],
        "v_y_code": v_cart2[..., 1],
        "v_z_code": v_cart2[..., 2],
        "speed_km_s": speed_km_s,
        "alfven_speed_classical_km_s": alfven_classical,
        "alfven_speed_boris_limited_km_s": alfven_limited,
        "plasma_beta": plasma_beta,
    }

    rsun_mm = units["unit_length_cm"] / 1.0e8
    if grid.mode == "phi":
        arrays["magnetic_u_plot"] = -b_sph2[..., 1] / radius * 180.0 / np.pi
        arrays["magnetic_v_plot"] = b_sph2[..., 0] * rsun_mm
        arrays["velocity_u_plot"] = -v_sph2[..., 1] / radius * 180.0 / np.pi
        arrays["velocity_v_plot"] = v_sph2[..., 0] * rsun_mm
    elif grid.mode == "theta":
        metric = radius * np.sin(theta)
        arrays["magnetic_u_plot"] = b_sph2[..., 2] / metric * 180.0 / np.pi
        arrays["magnetic_v_plot"] = b_sph2[..., 0] * rsun_mm
        arrays["velocity_u_plot"] = v_sph2[..., 2] / metric * 180.0 / np.pi
        arrays["velocity_v_plot"] = v_sph2[..., 0] * rsun_mm
    elif grid.mode == "r":
        metric = radius * np.sin(theta)
        arrays["magnetic_u_plot"] = b_sph2[..., 2] / metric * 180.0 / np.pi
        arrays["magnetic_v_plot"] = -b_sph2[..., 1] / radius * 180.0 / np.pi
        arrays["velocity_u_plot"] = v_sph2[..., 2] / metric * 180.0 / np.pi
        arrays["velocity_v_plot"] = -v_sph2[..., 1] / radius * 180.0 / np.pi
    else:
        arrays["magnetic_u_plot"] = np.einsum(
            "...i,i->...", b_cart2, grid.basis_u
        )
        arrays["magnetic_v_plot"] = np.einsum(
            "...i,i->...", b_cart2, grid.basis_v
        )
        arrays["velocity_u_plot"] = np.einsum(
            "...i,i->...", v_cart2, grid.basis_u
        )
        arrays["velocity_v_plot"] = np.einsum(
            "...i,i->...", v_cart2, grid.basis_v
        )

    height_mm = (radius - _as_float(par.get("xprobmin1"), 1.0)) * rsun_mm
    log_temperature = np.log10(np.maximum(temperature, 1.0))
    cooling_curve = str(par.get("coolcurve", "JCcorona"))
    if cooling_curve.lower() == "jccorona":
        log_lambda = np.interp(
            log_temperature, LOG10_T_JCCORONA, LOG10_LAMBDA_JCCORONA
        )
        lambda_cgs = 10.0**log_lambda
    else:
        lambda_cgs = np.full_like(temperature, np.nan)
    qrad_curve = (
        lambda_cgs * (1.0 + 2.0 * HE_ABUNDANCE) * number_density**2
    )
    cutoff_K = _as_float(par.get("tcool_cutoff"), 0.0) * units["unit_temperature_K"]
    qrad_curve = np.where(temperature > cutoff_K, qrad_curve, 0.0)
    if _as_bool(par.get("rad_damp"), False):
        damp_height_mm = _as_float(par.get("rad_damp_height"), 0.0) * rsun_mm
        damp_scale_mm = _as_float(par.get("rad_damp_scale"), 1.0) * rsun_mm
        damping = np.ones_like(height_mm)
        low = height_mm <= damp_height_mm
        damping[low] = np.exp(
            -((height_mm[low] - damp_height_mm) / damp_scale_mm) ** 2
        )
    else:
        damping = np.ones_like(height_mm)
    qrad_model = qrad_curve * damping
    radiation_active = _as_bool(par.get("mhd_radiative_cooling"), False)
    arrays["radiative_loss_curve_cgs"] = qrad_curve
    arrays["radiative_damping_factor"] = damping
    arrays["radiative_loss_model_cgs"] = qrad_model
    arrays["radiative_loss_active_estimate_cgs"] = (
        qrad_model if radiation_active else np.zeros_like(qrad_model)
    )
    heating0 = _as_float(par.get("background_heating_cgs"), 0.0)
    heating_scale_mm = 60.0
    arrays["background_heating_cgs"] = heating0 * np.exp(
        -height_mm / heating_scale_mm
    )
    with np.errstate(divide="ignore", invalid="ignore"):
        arrays["radiative_to_background_heating"] = (
            arrays["radiative_loss_active_estimate_cgs"]
            / arrays["background_heating_cgs"]
        )
        arrays["radiative_cooling_time_s"] = (
            e2
            * units["unit_pressure_cgs"]
            / arrays["radiative_loss_active_estimate_cgs"]
        )

    derivative_note = "not computed"
    derivative_details: dict[str, Any] = {"backend": derivative_backend}
    if derivative_backend == "uniform-slab":
        if len(grid.w) < 3:
            raise ValueError(
                "uniform-slab derivatives require --slab-points >= 3"
            )
        if grid.mode == "plane":
            current_full, divb_full = cartesian_plane_vector_derivatives(
                magnetic_cartesian, grid, rsun_mm
            )
            _, divv_full = cartesian_plane_vector_derivatives(
                velocity_cartesian, grid, rsun_mm
            )
            current2 = _central(current_full)
            current_cartesian2 = current2
        else:
            current_full, divb_full = spherical_vector_derivatives(
                magnetic_spherical, grid
            )
            _, divv_full = spherical_vector_derivatives(
                velocity_spherical, grid
            )
            current2 = _central(current_full)
            current_cartesian2 = spherical_vectors_to_cartesian(
                current2[..., 0],
                current2[..., 1],
                current2[..., 2],
                theta,
                phi,
            )
            arrays["j_r_code"] = current2[..., 0]
            arrays["j_theta_code"] = current2[..., 1]
            arrays["j_phi_code"] = current2[..., 2]
        divb = _central(divb_full)
        divv = _central(divv_full)
        current_magnitude = np.linalg.norm(current2, axis=-1)
        eta_code = _as_float(par.get("mhd_eta"), 0.0)
        arrays["j_x_code"] = current_cartesian2[..., 0]
        arrays["j_y_code"] = current_cartesian2[..., 1]
        arrays["j_z_code"] = current_cartesian2[..., 2]
        arrays["current_magnitude_code"] = current_magnitude
        arrays["divb_code"] = divb
        arrays["divv_code"] = divv
        arrays["divv_per_s"] = divv / units["unit_time_s"]
        derivative_note = (
            "computed from the uniform 3-D slab; divB is cell-centred "
            "post-processing and is not the exact CT face-centred get_divb diagnostic"
        )
    elif derivative_backend == "local-fit":
        local_arrays, derivative_details = local_fit_derivatives(
            dataset,
            grid,
            raw,
            par,
            units,
            fit_scales,
            fit_residual_max,
            fit_sensitivity_max,
        )
        arrays.update(local_arrays)
        derivative_note = (
            "computed by multi-scale symmetric local least-squares fits using "
            "each query point's true owner-cell dr/dtheta/dphi; quality and "
            "AMR-crossing masks are stored with the cache"
        )
    elif derivative_backend != "none":
        raise ValueError(f"Unknown derivative backend: {derivative_backend}")

    if "current_magnitude_code" in arrays:
        eta_code = _as_float(par.get("mhd_eta"), 0.0)
        arrays["j_over_b_per_Mm"] = (
            arrays["current_magnitude_code"] / bmag / rsun_mm
        )
        arrays["abs_divb_over_b_per_Mm"] = (
            np.abs(arrays["divb_code"]) / bmag / rsun_mm
        )
        arrays["ohmic_heating_code"] = (
            eta_code * arrays["current_magnitude_code"] ** 2
        )
        arrays["ohmic_heating_cgs"] = (
            arrays["ohmic_heating_code"]
            * units["unit_pressure_cgs"]
            / units["unit_time_s"]
        )
        if "derivative_reliable_mask" in arrays:
            reliable = arrays["derivative_reliable_mask"].astype(bool)
            strict_reliable = arrays.get(
                "derivative_strict_reliable_mask",
                arrays["derivative_reliable_mask"],
            ).astype(bool)
            velocity_reliable = arrays[
                "velocity_derivative_reliable_mask"
            ].astype(bool)
            velocity_strict_reliable = arrays.get(
                "velocity_derivative_strict_reliable_mask",
                arrays["velocity_derivative_reliable_mask"],
            ).astype(bool)
            for source, target in (
                ("current_magnitude_code", "current_magnitude_reliable_code"),
                ("j_over_b_per_Mm", "j_over_b_reliable_per_Mm"),
                ("divb_code", "divb_reliable_code"),
                ("abs_divb_over_b_per_Mm", "abs_divb_over_b_reliable_per_Mm"),
                ("ohmic_heating_cgs", "ohmic_heating_reliable_cgs"),
            ):
                arrays[target] = np.where(reliable, arrays[source], np.nan)
                arrays[target.replace("_reliable_", "_strict_reliable_")] = (
                    np.where(strict_reliable, arrays[source], np.nan)
                )
            arrays["divv_reliable_code"] = np.where(
                velocity_reliable, arrays["divv_code"], np.nan
            )
            arrays["divv_strict_reliable_code"] = np.where(
                velocity_strict_reliable, arrays["divv_code"], np.nan
            )

    physics = {
        "gamma": gamma,
        "velocity_reconstruction": velocity_method,
        "radiation_active": radiation_active,
        "coolcurve": cooling_curve,
        "tcool_cutoff_K": cutoff_K,
        "rad_damp": _as_bool(par.get("rad_damp"), False),
        "background_heating_cgs_at_bottom": heating0,
        "background_heating_scale_Mm": heating_scale_mm,
        "mhd_eta_code": _as_float(par.get("mhd_eta"), 0.0),
        "mhd_trac": _as_bool(par.get("mhd_trac"), False),
        "mhd_trac_type": par.get("mhd_trac_type"),
        "radiation_postprocess_note": (
            "JCcorona table, cutoff, and configured radial damping are reproduced. "
            "Dynamic TRAC cutoff/source bookkeeping is not stored in DAT, so this "
            "is a model estimate rather than the exact time-integrated solver source."
        ),
        "derivative_note": derivative_note,
        "derivative_details": derivative_details,
    }
    return arrays, physics


def _file_identity(path: Path, with_hash: bool) -> dict[str, Any]:
    stat = path.stat()
    identity: dict[str, Any] = {
        "path": str(path.resolve()),
        "size_bytes": int(stat.st_size),
        "mtime_ns": int(stat.st_mtime_ns),
    }
    if with_hash:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(16 * 1024 * 1024), b""):
                digest.update(block)
        identity["sha256"] = digest.hexdigest()
    return identity


def write_cache(
    output: Path, arrays: dict[str, np.ndarray], metadata: dict[str, Any]
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    payload = {name: np.asarray(value) for name, value in arrays.items()}
    payload["__metadata_json__"] = np.asarray(
        json.dumps(metadata, ensure_ascii=False, sort_keys=True)
    )
    np.savez_compressed(output, **payload)


def load_cache(path: Path) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    with np.load(path, allow_pickle=False) as archive:
        metadata = json.loads(str(archive["__metadata_json__"]))
        arrays = {
            key: np.asarray(archive[key])
            for key in archive.files
            if key != "__metadata_json__"
        }
    return arrays, metadata


def finite_summary(values: np.ndarray) -> dict[str, float | int]:
    finite = np.asarray(values)[np.isfinite(values)]
    if finite.size == 0:
        return {"count": 0, "min": math.nan, "median": math.nan, "p99": math.nan, "max": math.nan}
    return {
        "count": int(finite.size),
        "min": float(np.min(finite)),
        "median": float(np.median(finite)),
        "p99": float(np.percentile(finite, 99.0)),
        "max": float(np.max(finite)),
    }


FIELD_STYLES: dict[str, tuple[str, str, str]] = {
    "temperature_K": ("Temperature [K]", "inferno", "log"),
    "number_density_cm3": (r"Number density [cm$^{-3}$]", "viridis", "log"),
    "mass_density_g_cm3": (r"Mass density [g cm$^{-3}$]", "viridis", "log"),
    "pressure_cgs": (r"Pressure [dyn cm$^{-2}$]", "magma", "log"),
    "speed_km_s": ("Speed [km/s]", "turbo", "log"),
    "alfven_speed_classical_km_s": ("Classical Alfven speed [km/s]", "plasma", "log"),
    "alfven_speed_boris_limited_km_s": ("Boris-limited Alfven diagnostic [km/s]", "plasma", "log"),
    "plasma_beta": ("Plasma beta", "cividis", "log"),
    "b_magnitude_G": ("|B| [G]", "cividis", "log"),
    "radiative_loss_model_cgs": (r"Model radiative loss [erg cm$^{-3}$ s$^{-1}$]", "magma", "log"),
    "radiative_loss_active_estimate_cgs": (r"Active radiative-loss estimate [erg cm$^{-3}$ s$^{-1}$]", "magma", "log"),
    "background_heating_cgs": (r"Background heating [erg cm$^{-3}$ s$^{-1}$]", "YlGnBu_r", "log"),
    "radiative_to_background_heating": (r"$Q_{\rm rad}/Q_{\rm bg}$", "coolwarm", "log"),
    "radiative_cooling_time_s": ("Radiative cooling time [s]", "viridis", "log"),
    "current_magnitude_code": ("|J| [code]", "magma", "log"),
    "current_magnitude_reliable_code": ("Reliable |J| [code]", "magma", "log"),
    "current_magnitude_strict_reliable_code": ("Strict reliable |J| [code]", "magma", "log"),
    "j_over_b_per_Mm": (r"|J|/|B| [Mm$^{-1}$]", "magma", "log"),
    "j_over_b_reliable_per_Mm": (r"Reliable |J|/|B| [Mm$^{-1}$]", "magma", "log"),
    "j_over_b_strict_reliable_per_Mm": (r"Strict reliable |J|/|B| [Mm$^{-1}$]", "magma", "log"),
    "ohmic_heating_cgs": (r"Ohmic heating [erg cm$^{-3}$ s$^{-1}$]", "inferno", "log"),
    "ohmic_heating_reliable_cgs": (r"Reliable Ohmic heating [erg cm$^{-3}$ s$^{-1}$]", "inferno", "log"),
    "ohmic_heating_strict_reliable_cgs": (r"Strict reliable Ohmic heating [erg cm$^{-3}$ s$^{-1}$]", "inferno", "log"),
    "divb_code": (r"Cell-centred $\nabla\cdot B$ [code]", "RdBu_r", "symlog"),
    "divb_reliable_code": (r"Reliable cell-centred $\nabla\cdot B$ [code]", "RdBu_r", "symlog"),
    "divb_strict_reliable_code": (r"Strict reliable cell-centred $\nabla\cdot B$ [code]", "RdBu_r", "symlog"),
    "abs_divb_over_b_per_Mm": (r"$|\nabla\cdot B|/|B|$ [Mm$^{-1}$]", "magma", "log"),
    "abs_divb_over_b_reliable_per_Mm": (r"Reliable $|\nabla\cdot B|/|B|$ [Mm$^{-1}$]", "magma", "log"),
    "abs_divb_over_b_strict_reliable_per_Mm": (r"Strict reliable $|\nabla\cdot B|/|B|$ [Mm$^{-1}$]", "magma", "log"),
    "divv_code": (r"$\nabla\cdot v$ [code]", "RdBu_r", "symlog"),
    "divv_reliable_code": (r"Reliable $\nabla\cdot v$ [code]", "RdBu_r", "symlog"),
    "divv_strict_reliable_code": (r"Strict reliable $\nabla\cdot v$ [code]", "RdBu_r", "symlog"),
    "derivative_fit_residual": ("Magnetic derivative fit residual", "viridis", "linear"),
    "derivative_scale_sensitivity": ("Magnetic derivative scale sensitivity", "viridis", "linear"),
    "owner_amr_level": ("Owner AMR level", "viridis", "linear"),
    "v_r_code": (r"$v_r$ [code]", "RdBu_r", "symlog"),
    "v_theta_code": (r"$v_\theta$ [code]", "RdBu_r", "symlog"),
    "v_phi_code": (r"$v_\phi$ [code]", "RdBu_r", "symlog"),
}


def plot_native_amr_grid(
    dataset: Any,
    dat_path: Path,
    par_path: Path,
    output: Path,
    mode: str,
    value_deg: float,
    min_height_mm: float,
    max_height_mm: float,
    u_min: float | None,
    u_max: float | None,
    dpi: int,
    unit_length_cm: float,
    unit_time_s: float,
    marker_u: float | None = None,
    marker_v: float | None = None,
    marker_label: str = "reference point",
    hide_cell_lines: bool = False,
    hash_input: bool = False,
) -> dict[str, Any]:
    """Plot native AMR leaves and cells intersecting a fixed angular plane.

    Unlike the uniform slice renderer, this function does not interpolate the
    mesh.  It draws the true leaf-block and cell nodes reconstructed by simesh
    from the DAT tree together with ``amrvac.par``.
    """

    import matplotlib.pyplot as plt
    from matplotlib.collections import LineCollection, PatchCollection
    from matplotlib.lines import Line2D
    from matplotlib.patches import Rectangle

    if mode not in {"phi", "theta"}:
        raise ValueError("native AMR grid plots support mode=phi or mode=theta")
    if int(dataset.ndim) != 3:
        raise ValueError("a fixed-coordinate AMR grid plot requires a 3-D DAT")
    if not str(dataset.metadata.get("geometry", "")).lower().startswith("spherical"):
        raise ValueError("native fixed-phi/theta AMR grid plots require spherical geometry")
    if not np.isfinite(value_deg):
        raise ValueError("--value must be finite")
    if not 0.0 <= min_height_mm < max_height_mm:
        raise ValueError("height limits must satisfy 0 <= MIN < MAX")
    if (marker_u is None) != (marker_v is None):
        raise ValueError("--marker-u and --marker-v must be supplied together")

    levels = np.asarray(dataset.tree_info[0], dtype=int)
    if levels.ndim != 1:
        raise ValueError("unexpected AMR level array shape")
    r_nodes, theta_nodes, phi_nodes = (
        np.asarray(nodes, dtype=float)
        for nodes in dataset.block_node_coordinates
    )
    if not (len(levels) == len(r_nodes) == len(theta_nodes) == len(phi_nodes)):
        raise ValueError("AMR tree and block-node arrays have inconsistent lengths")

    rsun_mm = float(unit_length_cm) / 1.0e8
    height_nodes = (r_nodes - 1.0) * rsun_mm
    height_lo = height_nodes[:, 0]
    height_hi = height_nodes[:, -1]
    fixed_index = 2 if mode == "phi" else 1
    fixed_nodes = phi_nodes if mode == "phi" else theta_nodes
    fixed_value = math.radians(value_deg)
    fixed_domain = (
        float(dataset.metadata["xmin"][fixed_index]),
        float(dataset.metadata["xmax"][fixed_index]),
    )
    tolerance = max(abs(fixed_domain[1] - fixed_domain[0]) * 1.0e-12, 1.0e-14)
    if not fixed_domain[0] - tolerance <= fixed_value <= fixed_domain[1] + tolerance:
        raise ValueError(
            f"requested {mode}={value_deg:g} deg lies outside the DAT domain"
        )
    intersects_fixed = (fixed_value >= fixed_nodes[:, 0] - tolerance) & (
        fixed_value < fixed_nodes[:, -1] - tolerance
    )
    if math.isclose(fixed_value, fixed_domain[1], abs_tol=tolerance):
        intersects_fixed |= np.isclose(
            fixed_nodes[:, -1], fixed_domain[1], atol=tolerance, rtol=0.0
        )

    if mode == "phi":
        horizontal_nodes = 90.0 - np.degrees(theta_nodes)
        default_u_limits = (
            90.0 - math.degrees(float(dataset.metadata["xmax"][1])),
            90.0 - math.degrees(float(dataset.metadata["xmin"][1])),
        )
        xlabel = "Latitude [deg]"
        fixed_description = rf"fixed $\phi={value_deg:.4f}^\circ$"
    else:
        horizontal_nodes = np.degrees(phi_nodes)
        default_u_limits = (
            math.degrees(float(dataset.metadata["xmin"][2])),
            math.degrees(float(dataset.metadata["xmax"][2])),
        )
        xlabel = r"Longitude $\phi$ [deg]"
        latitude_deg = 90.0 - value_deg
        fixed_description = (
            rf"fixed $\theta={value_deg:.4f}^\circ$ "
            rf"(latitude ${latitude_deg:.4f}^\circ$)"
        )
    horizontal_nodes = np.sort(horizontal_nodes, axis=1)
    horizontal_lo = horizontal_nodes[:, 0]
    horizontal_hi = horizontal_nodes[:, -1]
    view_u_min = default_u_limits[0] if u_min is None else float(u_min)
    view_u_max = default_u_limits[1] if u_max is None else float(u_max)
    if not np.isfinite(view_u_min) or not np.isfinite(view_u_max):
        raise ValueError("horizontal limits must be finite")
    if view_u_min >= view_u_max:
        raise ValueError("horizontal limits must satisfy MIN < MAX")
    view_u_min = max(view_u_min, default_u_limits[0])
    view_u_max = min(view_u_max, default_u_limits[1])
    if view_u_min >= view_u_max:
        raise ValueError("horizontal limits do not overlap the DAT domain")

    selected = np.flatnonzero(
        intersects_fixed
        & (height_hi >= min_height_mm)
        & (height_lo <= max_height_mm)
        & (horizontal_hi >= view_u_min)
        & (horizontal_lo <= view_u_max)
    )
    if selected.size == 0:
        raise RuntimeError("no AMR leaf blocks intersect the requested plane and view")

    fig, axis = plt.subplots(figsize=(12.0, 7.5), constrained_layout=True)
    slice_counts: dict[str, int] = {}
    block_width_stats: dict[str, Any] = {}
    handles: list[Any] = []
    for level in sorted(np.unique(levels)):
        level_selected = selected[levels[selected] == level]
        level_key = str(int(level))
        slice_counts[level_key] = int(level_selected.size)
        if level_selected.size == 0:
            continue
        patches = []
        cell_segments = []
        for leaf in level_selected:
            u0 = max(view_u_min, float(horizontal_lo[leaf]))
            u1 = min(view_u_max, float(horizontal_hi[leaf]))
            h0 = max(min_height_mm, float(height_lo[leaf]))
            h1 = min(max_height_mm, float(height_hi[leaf]))
            if u1 <= u0 or h1 <= h0:
                continue
            patches.append(Rectangle((u0, h0), u1 - u0, h1 - h0))
            if not hide_cell_lines:
                for u_value in horizontal_nodes[leaf, 1:-1]:
                    if u0 < u_value < u1:
                        cell_segments.append(((float(u_value), h0), (float(u_value), h1)))
                for h_value in height_nodes[leaf, 1:-1]:
                    if h0 < h_value < h1:
                        cell_segments.append(((u0, float(h_value)), (u1, float(h_value))))

        color = AMR_LEVEL_COLORS.get(int(level), "#555555")
        axis.add_collection(
            PatchCollection(
                patches,
                facecolor=color,
                edgecolor=color,
                linewidth=0.48,
                alpha=0.23,
                zorder=2 + int(level),
            )
        )
        if cell_segments:
            axis.add_collection(
                LineCollection(
                    cell_segments,
                    colors="#30343B",
                    linewidths=0.16,
                    alpha=0.20,
                    zorder=7,
                )
            )
        handles.append(
            Line2D(
                [0],
                [0],
                color=color,
                linewidth=5,
                label=f"L{int(level)}: {len(level_selected)} leaf blocks",
            )
        )
        horizontal_width = horizontal_hi[level_selected] - horizontal_lo[level_selected]
        height_width = height_hi[level_selected] - height_lo[level_selected]
        block_width_stats[level_key] = {
            "horizontal_block_width_deg": [
                float(np.min(horizontal_width)),
                float(np.max(horizontal_width)),
            ],
            "height_block_width_Mm": [
                float(np.min(height_width)),
                float(np.max(height_width)),
            ],
        }

    marker = None
    if marker_u is not None and marker_v is not None:
        marker = {
            "u": float(marker_u),
            "v_Mm": float(marker_v),
            "label": marker_label,
        }
        axis.plot(
            marker_u,
            marker_v,
            marker="x",
            color="#00B7C7",
            markersize=9,
            markeredgewidth=2.0,
            linestyle="none",
            zorder=10,
        )
        handles.append(
            Line2D(
                [0],
                [0],
                marker="x",
                color="#00B7C7",
                markersize=8,
                markeredgewidth=1.8,
                linestyle="none",
                label=marker_label,
            )
        )

    axis.set_xlim(view_u_min, view_u_max)
    axis.set_ylim(min_height_mm, max_height_mm)
    axis.set_xlabel(xlabel)
    axis.set_ylabel("Height [Mm]")
    axis.set_facecolor("#FAFAF7")
    axis.grid(color="#60646C", alpha=0.15, linewidth=0.45, zorder=1)
    axis.legend(handles=handles, loc="upper right", framealpha=0.94, fontsize=9)

    time_code = float(dataset.metadata["time"])
    time_seconds = time_code * float(unit_time_s)
    cell_note = "leaf blocks and cells" if not hide_cell_lines else "leaf blocks"
    axis.set_title(
        f"Native AMR grid, {dat_path.stem[-4:]}, {fixed_description}\n"
        f"t={time_code:.6f} ({time_seconds:.2f} s); {cell_note}"
    )

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    total_counts = {
        str(int(level)): int(np.count_nonzero(levels == level))
        for level in sorted(np.unique(levels))
    }
    report = {
        "tool": {
            "name": "spherical_slice_analysis.py",
            "version": VERSION,
            "command": "plot-amr-grid",
            "created_utc": datetime.now(timezone.utc).isoformat(),
        },
        "dataset": {
            **_file_identity(dat_path, hash_input),
            "time_code": time_code,
            "time_seconds": time_seconds,
            "iteration": int(dataset.metadata["it"]),
            "geometry": str(dataset.metadata["geometry"]),
            "coordinate_source": dataset.metadata.get("coordinate_source"),
        },
        "grid_config": _file_identity(par_path, hash_input),
        "slice": {
            "mode": mode,
            "value_deg": float(value_deg),
            "horizontal_coordinate": "latitude_deg" if mode == "phi" else "longitude_deg",
            "horizontal_limits_deg": [view_u_min, view_u_max],
            "height_limits_Mm": [min_height_mm, max_height_mm],
            "marker": marker,
            "cell_lines_visible": not hide_cell_lines,
        },
        "block_nx": [int(value) for value in dataset.block_nx],
        "total_leaf_blocks": int(len(levels)),
        "total_leaf_counts_by_level": total_counts,
        "slice_leaf_blocks": int(selected.size),
        "slice_leaf_counts_by_level": slice_counts,
        "block_width_stats": block_width_stats,
        "output": _file_identity(output, True),
    }
    output.with_suffix(".json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


def command_plot_amr_grid(args: argparse.Namespace) -> int:
    from simesh.amrvac import open_dataset

    if not args.dat.is_file():
        raise FileNotFoundError(args.dat)
    if not args.amrvac_par.is_file():
        raise FileNotFoundError(args.amrvac_par)
    dataset = open_dataset(
        str(args.dat), ghost_width=0, grid_config=str(args.amrvac_par)
    )
    report = plot_native_amr_grid(
        dataset=dataset,
        dat_path=args.dat,
        par_path=args.amrvac_par,
        output=args.output,
        mode=args.mode,
        value_deg=args.value,
        min_height_mm=args.min_height_mm,
        max_height_mm=args.max_height_mm,
        u_min=args.u_min,
        u_max=args.u_max,
        dpi=args.dpi,
        unit_length_cm=args.unit_length_cm,
        unit_time_s=args.unit_time_s,
        marker_u=args.marker_u,
        marker_v=args.marker_v,
        marker_label=args.marker_label,
        hide_cell_lines=args.hide_cell_lines,
        hash_input=args.hash_input,
    )
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0


def _positive_limits(values: np.ndarray, requested: tuple[float, float] | None) -> tuple[float, float]:
    if requested is not None:
        return requested
    finite = values[np.isfinite(values) & (values > 0.0)]
    if finite.size == 0:
        return 1.0, 10.0
    low, high = np.percentile(finite, [1.0, 99.5])
    low = max(float(low), np.finfo(float).tiny)
    high = max(float(high), low * 1.001)
    return low, high


def plot_cache(
    cache_path: Path,
    output: Path,
    fields: list[str],
    overlay: str,
    limits: dict[str, tuple[float, float]],
    dpi: int,
) -> dict[str, Any]:
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm, SymLogNorm

    arrays, metadata = load_cache(cache_path)
    missing = [field for field in fields if field not in arrays]
    if missing:
        raise KeyError(
            f"Cache does not contain {missing}. Extract the cache with "
            "--derivative-backend local-fit, or use uniform-slab together "
            "with --slab-points >= 3."
        )
    columns = min(3, len(fields))
    rows = math.ceil(len(fields) / columns)
    fig, axes = plt.subplots(
        rows,
        columns,
        figsize=(6.4 * columns, 5.2 * rows),
        squeeze=False,
        constrained_layout=True,
    )
    u = arrays["u"]
    v = arrays["v"]
    metrics: dict[str, Any] = {}
    for axis, field in zip(axes.flat, fields):
        values = np.asarray(arrays[field])
        if values.dtype == np.bool_:
            values = values.astype(float)
        label, cmap, scale = FIELD_STYLES.get(field, (field, "viridis", "linear"))
        field_limits = limits.get(field)
        if scale == "log":
            low, high = _positive_limits(values, field_limits)
            norm = LogNorm(low, high)
            plotted = np.ma.masked_where(~np.isfinite(values) | (values <= 0.0), values)
        elif scale == "symlog":
            finite = np.abs(values[np.isfinite(values)])
            maximum = (
                float(np.percentile(finite, 99.5)) if finite.size else 1.0
            )
            if field_limits is not None:
                maximum = max(abs(field_limits[0]), abs(field_limits[1]))
            maximum = max(maximum, np.finfo(float).tiny)
            norm = SymLogNorm(linthresh=maximum * 1.0e-3, vmin=-maximum, vmax=maximum)
            plotted = np.ma.masked_invalid(values)
        else:
            norm = None
            plotted = np.ma.masked_invalid(values)
        image = axis.pcolormesh(
            u, v, plotted, shading="auto", cmap=cmap, norm=norm, rasterized=True
        )
        fig.colorbar(image, ax=axis, pad=0.02, label=label)
        if overlay != "none":
            prefix = "magnetic" if overlay == "magnetic" else "velocity"
            flow_u = np.asarray(arrays[f"{prefix}_u_plot"])
            flow_v = np.asarray(arrays[f"{prefix}_v_plot"])
            magnitude = np.hypot(flow_u, flow_v)
            valid = np.isfinite(magnitude) & (magnitude > 0.0)
            axis.streamplot(
                u,
                v,
                np.ma.masked_where(~valid, flow_u),
                np.ma.masked_where(~valid, flow_v),
                color="white",
                density=(1.2, 1.25),
                linewidth=0.65,
                arrowsize=0.7,
                minlength=0.1,
                maxlength=4.0,
                integration_direction="both",
                broken_streamlines=True,
            )
        axis.set_xlabel(metadata["surface"]["xlabel"])
        axis.set_ylabel(metadata["surface"]["ylabel"])
        axis.set_title(field)
        axis.set_xlim(float(u[0]), float(u[-1]))
        axis.set_ylim(float(v[0]), float(v[-1]))
        axis.grid(color="white", alpha=0.10, linewidth=0.35)
        metrics[field] = finite_summary(values)
    for axis in axes.flat[len(fields):]:
        axis.set_visible(False)
    time_code = metadata["dataset"]["time_code"]
    time_seconds = metadata["dataset"]["time_seconds"]
    fig.suptitle(
        f"{metadata['surface']['value_description']}; "
        f"t={time_code:.6g} ({time_seconds:.2f} s); overlay={overlay}",
        fontsize=13,
    )
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=dpi, bbox_inches="tight")
    plt.close(fig)
    report = {
        "cache": str(cache_path.resolve()),
        "output": str(output.resolve()),
        "fields": fields,
        "overlay": overlay,
        "metrics": metrics,
    }
    output.with_suffix(".json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


def _parse_limits(items: Iterable[str]) -> dict[str, tuple[float, float]]:
    result: dict[str, tuple[float, float]] = {}
    for item in items:
        parts = item.split(":")
        if len(parts) != 3:
            raise ValueError("--limit must be FIELD:MIN:MAX")
        result[parts[0]] = (float(parts[1]), float(parts[2]))
    return result


def _import_legacy_script(path: Path) -> Any:
    spec = importlib.util.spec_from_file_location("legacy_slice_script", path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load legacy script: {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def validate_fixed_phi_against_legacy(
    cache_path: Path,
    legacy_script: Path,
    dat_path: Path,
    par_path: Path,
    output: Path,
) -> dict[str, Any]:
    from simesh.amrvac import open_dataset

    arrays, metadata = load_cache(cache_path)
    if metadata["surface"]["mode"] != "phi":
        raise ValueError("Legacy validation currently applies only to mode=phi")
    legacy = _import_legacy_script(legacy_script)
    dataset = open_dataset(str(dat_path), ghost_width=0, grid_config=str(par_path))
    phi_deg = float(metadata["surface"]["value"])
    height = arrays["v"]
    latitude = arrays["u"]
    result = legacy.sample_slice(
        dataset,
        phi_deg,
        float(height[-1]),
        int(height.size),
        int(latitude.size),
    )
    legacy_latitude, legacy_height, _, legacy_fields, legacy_fallbacks = result
    gamma = _dataset_parameter(dataset, "gamma", 5.0 / 3.0)
    legacy_temperature = (
        (gamma - 1.0)
        * legacy_fields["e"]
        / legacy_fields["rho"]
        * DEFAULT_UNITS["unit_temperature_K"]
    )
    legacy_density = (
        legacy_fields["rho"] * DEFAULT_UNITS["unit_numberdensity_cm3"]
    )

    def comparison(reference: np.ndarray, candidate: np.ndarray) -> dict[str, float]:
        difference = candidate - reference
        scale = np.maximum(np.abs(reference), np.finfo(float).tiny)
        return {
            "max_abs": float(np.nanmax(np.abs(difference))),
            "rms_abs": float(np.sqrt(np.nanmean(difference**2))),
            "max_relative": float(np.nanmax(np.abs(difference) / scale)),
            "rms_relative": float(np.sqrt(np.nanmean((difference / scale) ** 2))),
        }

    report = {
        "cache": str(cache_path.resolve()),
        "legacy_script": str(legacy_script.resolve()),
        "coordinates": {
            "latitude_max_abs_deg": float(np.max(np.abs(arrays["u"] - legacy_latitude))),
            "height_max_abs_Mm": float(np.max(np.abs(arrays["v"] - legacy_height))),
        },
        "rho_code": comparison(legacy_fields["rho"], arrays["rho_code"]),
        "internal_energy_code": comparison(
            legacy_fields["e"], arrays["internal_energy_code"]
        ),
        "temperature_K": comparison(legacy_temperature, arrays["temperature_K"]),
        "number_density_cm3": comparison(
            legacy_density, arrays["number_density_cm3"]
        ),
        "legacy_nearest_fallback_samples": int(legacy_fallbacks),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return report


def make_image_comparison(old: Path, new: Path, output: Path) -> None:
    import matplotlib.image as mpimg
    import matplotlib.pyplot as plt

    images = [mpimg.imread(old), mpimg.imread(new)]
    fig, axes = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True)
    for axis, image, title in zip(axes, images, ("Legacy verified script", "Unified tool")):
        axis.imshow(image)
        axis.set_title(title, fontsize=14)
        axis.axis("off")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=160, bbox_inches="tight")
    plt.close(fig)


def add_extract_arguments(parser: argparse.ArgumentParser) -> None:
    parser.add_argument("dat", type=Path)
    parser.add_argument("amrvac_par", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--mode", choices=("phi", "theta", "r", "plane"), default="phi")
    parser.add_argument(
        "--value",
        type=float,
        help="phi/theta in degrees, or r in R_sun; not used for mode=plane",
    )
    parser.add_argument("--nu", type=int, default=420)
    parser.add_argument("--nv", type=int, default=360)
    parser.add_argument("--u-min", type=float)
    parser.add_argument("--u-max", type=float)
    parser.add_argument("--v-min", type=float)
    parser.add_argument("--v-max", type=float)
    parser.add_argument("--max-height-mm", type=float, default=120.0)
    parser.add_argument(
        "--slab-points",
        type=int,
        default=1,
        help="odd number; >=3 is required only by the uniform-slab backend",
    )
    parser.add_argument(
        "--slab-half-width",
        type=float,
        help="degrees for phi/theta; Mm for r/plane (default 0.1)",
    )
    parser.add_argument(
        "--derivative-backend",
        choices=("auto", "none", "uniform-slab", "local-fit"),
        default="auto",
        help=(
            "auto uses uniform-slab when slab-points>=3, otherwise none; "
            "local-fit is recommended for quantitative existing-DAT diagnostics"
        ),
    )
    parser.add_argument(
        "--fit-scales",
        default="1.5,2.0,2.5",
        help="local-fit symmetric offsets in owner-cell widths",
    )
    parser.add_argument("--fit-residual-max", type=float, default=0.35)
    parser.add_argument("--fit-sensitivity-max", type=float, default=0.35)
    parser.add_argument("--plane-origin-spherical", help="r_Rsun,latitude_deg,longitude_deg")
    parser.add_argument("--plane-origin-cartesian", help="x,y,z in R_sun")
    parser.add_argument("--plane-normal", help="Cartesian nx,ny,nz")
    parser.add_argument("--plane-up", help="preferred in-plane Cartesian up vector")
    parser.add_argument("--plane-width-mm", type=float, default=200.0)
    parser.add_argument("--plane-height-mm", type=float, default=160.0)
    parser.add_argument("--interpolation", choices=("linear", "nearest"), default="linear")
    parser.add_argument("--chunk-size", type=int, default=200000)
    parser.add_argument("--hash-input", action="store_true")
    parser.add_argument("--unit-length-cm", type=float, default=DEFAULT_UNITS["unit_length_cm"])


def command_extract(args: argparse.Namespace) -> int:
    from simesh.amrvac import open_dataset

    if not args.dat.is_file():
        raise FileNotFoundError(args.dat)
    if not args.amrvac_par.is_file():
        raise FileNotFoundError(args.amrvac_par)
    par = parse_amrvac_par(args.amrvac_par)
    units = dict(DEFAULT_UNITS)
    units["unit_length_cm"] = float(args.unit_length_cm)
    units["unit_time_s"] = (
        units["unit_length_cm"] / units["unit_velocity_cm_s"]
    )
    dataset = open_dataset(
        str(args.dat), ghost_width=0, grid_config=str(args.amrvac_par)
    )
    if not getattr(dataset, "uses_stretched_coordinates", False):
        coordinate_note = (
            "simesh reports a uniform mesh. Verify this is intentional; the DAT "
            "header itself does not store stretching parameters."
        )
    else:
        coordinate_note = "true stretched coordinates reconstructed from amrvac.par"
    grid = build_sample_grid(args, dataset)
    raw, sample_stats = sample_raw_fields(
        dataset, grid, args.interpolation, args.chunk_size
    )
    derivative_backend = args.derivative_backend
    if derivative_backend == "auto":
        derivative_backend = (
            "uniform-slab" if int(args.slab_points) >= 3 else "none"
        )
    fit_scales = _parse_float_list(args.fit_scales, "--fit-scales")
    arrays, physics = derive_fields(
        raw,
        grid,
        dataset,
        par,
        units,
        derivative_backend,
        fit_scales,
        args.fit_residual_max,
        args.fit_sensitivity_max,
    )
    value = args.value
    if args.mode == "plane":
        value = None
    metadata = {
        "tool": {
            "name": "spherical_slice_analysis.py",
            "version": VERSION,
            "created_utc": datetime.now(timezone.utc).isoformat(),
        },
        "dataset": {
            **_file_identity(args.dat, args.hash_input),
            "time_code": float(dataset.metadata["time"]),
            "time_seconds": float(dataset.metadata["time"]) * units["unit_time_s"],
            "iteration": int(dataset.metadata["it"]),
            "geometry": str(dataset.metadata["geometry"]),
            "coordinate_source": dataset.metadata.get("coordinate_source"),
            "coordinate_note": coordinate_note,
            "wnames": [str(name) for name in dataset.wnames],
        },
        "grid_config": _file_identity(args.amrvac_par, args.hash_input),
        "surface": {
            "mode": args.mode,
            "value": value,
            "nu": int(args.nu),
            "nv": int(args.nv),
            "slab_points": int(args.slab_points),
            **grid.metadata,
        },
        "sampling": {
            "interpolation": args.interpolation,
            **sample_stats,
        },
        "units": units,
        "physics": physics,
        "available_fields": sorted(arrays),
    }
    write_cache(args.output, arrays, metadata)
    summary = {
        "output": str(args.output.resolve()),
        "surface": metadata["surface"],
        "sampling": metadata["sampling"],
        "physics": physics,
        "available_fields": sorted(arrays),
        "temperature_K": finite_summary(arrays["temperature_K"]),
        "number_density_cm3": finite_summary(arrays["number_density_cm3"]),
    }
    args.output.with_suffix(".summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Uniform fixed-coordinate and arbitrary-plane diagnostics for spherical AMRVAC DAT"
    )
    parser.add_argument("--version", action="version", version=VERSION)
    subparsers = parser.add_subparsers(dest="command", required=True)

    extract = subparsers.add_parser("extract", help="sample DAT and write a reusable NPZ cache")
    add_extract_arguments(extract)
    extract.set_defaults(func=command_extract)

    plot = subparsers.add_parser("plot", help="plot one or more fields from a cache")
    plot.add_argument("cache", type=Path)
    plot.add_argument("output", type=Path)
    plot.add_argument(
        "--fields",
        nargs="+",
        default=["temperature_K", "number_density_cm3"],
    )
    plot.add_argument("--overlay", choices=("none", "magnetic", "velocity"), default="magnetic")
    plot.add_argument("--limit", action="append", default=[], help="FIELD:MIN:MAX")
    plot.add_argument("--dpi", type=int, default=220)
    plot.set_defaults(
        func=lambda args: (
            print(
                json.dumps(
                    plot_cache(
                        args.cache,
                        args.output,
                        args.fields,
                        args.overlay,
                        _parse_limits(args.limit),
                        args.dpi,
                    ),
                    ensure_ascii=False,
                    indent=2,
                )
            )
            or 0
        )
    )

    amr_grid = subparsers.add_parser(
        "plot-amr-grid",
        help="plot native AMR leaf blocks and cells on a fixed phi/theta plane",
    )
    amr_grid.add_argument("dat", type=Path)
    amr_grid.add_argument("amrvac_par", type=Path)
    amr_grid.add_argument("output", type=Path)
    amr_grid.add_argument("--mode", choices=("phi", "theta"), required=True)
    amr_grid.add_argument(
        "--value",
        type=float,
        required=True,
        help="fixed phi or theta in degrees",
    )
    amr_grid.add_argument("--min-height-mm", type=float, default=0.0)
    amr_grid.add_argument("--max-height-mm", type=float, default=200.0)
    amr_grid.add_argument(
        "--u-min",
        type=float,
        help="minimum latitude for mode=phi or longitude for mode=theta",
    )
    amr_grid.add_argument(
        "--u-max",
        type=float,
        help="maximum latitude for mode=phi or longitude for mode=theta",
    )
    amr_grid.add_argument("--marker-u", type=float)
    amr_grid.add_argument("--marker-v", type=float)
    amr_grid.add_argument("--marker-label", default="reference point")
    amr_grid.add_argument("--hide-cell-lines", action="store_true")
    amr_grid.add_argument("--dpi", type=int, default=200)
    amr_grid.add_argument("--hash-input", action="store_true")
    amr_grid.add_argument(
        "--unit-length-cm",
        type=float,
        default=DEFAULT_UNITS["unit_length_cm"],
    )
    amr_grid.add_argument(
        "--unit-time-s",
        type=float,
        default=DEFAULT_UNITS["unit_time_s"],
    )
    amr_grid.set_defaults(func=command_plot_amr_grid)

    validate = subparsers.add_parser(
        "validate-legacy",
        help="numerically compare a fixed-phi cache with the verified legacy sampler",
    )
    validate.add_argument("cache", type=Path)
    validate.add_argument("legacy_script", type=Path)
    validate.add_argument("dat", type=Path)
    validate.add_argument("amrvac_par", type=Path)
    validate.add_argument("output", type=Path)
    validate.set_defaults(
        func=lambda args: (
            print(
                json.dumps(
                    validate_fixed_phi_against_legacy(
                        args.cache,
                        args.legacy_script,
                        args.dat,
                        args.amrvac_par,
                        args.output,
                    ),
                    ensure_ascii=False,
                    indent=2,
                )
            )
            or 0
        )
    )

    images = subparsers.add_parser(
        "compare-images", help="make a side-by-side legacy/new image sheet"
    )
    images.add_argument("legacy_image", type=Path)
    images.add_argument("new_image", type=Path)
    images.add_argument("output", type=Path)
    images.set_defaults(
        func=lambda args: (
            make_image_comparison(args.legacy_image, args.new_image, args.output)
            or 0
        )
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return int(args.func(args))
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        if "--debug" in (argv or sys.argv[1:]):
            raise
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Trace the tracked y=0 mask to bottom footpoints and measure its flux.

For every snapshot, the solid ``log10(Q)=4`` y=0 mask is used in two
independent measurements:

1. FastQSL2 interpolates ``B_y`` at every selected y=0 sample and the flux is
   integrated directly through that cross-section.
2. The same samples are traced to obtain a safe bottom-boundary search window.
   A fixed, fine bottom grid is then traced back through the volume.  A bottom
   sample is selected when its field line crosses the y=0 mask at least once.
   Each bottom sample is therefore counted once even when a line crosses the
   plane repeatedly or is reached by many y=0 seeds.

Only ``rboundary == 11`` lines enter the paired positive/negative bottom
footprint.  Geometrically filled mask pixels whose lines reach another
boundary are retained in the full y=0 flux and reported separately.

The default bottom spacing is one quarter of a native magnetic cell, matching
the saved QSL-mask spacing.  Magnetic flux is written in maxwells using
``B0=2 G`` and ``L0=10 Mm``.

Example
-------
/opt/anaconda3/bin/python \
  /Users/xiaoxingyun/Desktop/solar_code/workflows/Simulation_Analysis/scripts/\
trace_y0_mask_footpoints_flux_evolution.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import tempfile
from collections import Counter
from dataclasses import dataclass
from io import BytesIO
from pathlib import Path
from typing import Any, Iterable

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import ListedColormap
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
from PIL import Image
from scipy.ndimage import binary_fill_holes, distance_transform_edt, label

from calculate_fixed_y_qsl_slices import (
    DOMAIN,
    LENGTH_UNIT_MM,
    TIME_MINUTES_PER_SNAPSHOT,
    cell_centres,
    fastqsl,
    read_binary_data,
)


MAGNETIC_DIR = Path("/Volumes/Seagate/w2_jet_mhd/out_path1/nlfff")
SLICE_DIR = Path("/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/y0_qsl")
MASK_DIR = SLICE_DIR / "enclosed_region" / "evolution_logq4"
OUTPUT_DIR = MASK_DIR / "bottom_footpoint_flux"

START_SNAPSHOT = 10
STOP_SNAPSHOT = 95
TIME_MINUTES = TIME_MINUTES_PER_SNAPSHOT
FIXED_Y_CODE = 0.0
BOTTOM_DATUM_Z_CODE = 0.1
B_UNIT_G = 2.0
THREADS = 8
BOTTOM_SUBSAMPLE = 4

# Trace the complete cell-centre bottom plane.  This avoids losing remote
# high-Q branches during the t=21--25 connectivity transition and avoids
# defining the quantitative footprint by a convex hull or display operation.
BOTTOM_TRACE_LIMITS_MM = (-40.0, 80.0, -60.0, 60.0)
Y0_VIEW_LIMITS_MM = (-15.0, 30.0, 0.0, 30.0)

MASK_COLOR = "#43C7CA"
MASK_EDGE_COLOR = "#006B70"
CLOSED_COLOR = "#6A3D9A"
Q4_COLOR = "#E69F00"
POSITIVE_COLOR = "#CC79A7"
NEGATIVE_COLOR = "#0072B2"
BZ0_COLOR = "#D62728"

SCHEMA_VERSION = 1


@dataclass(frozen=True)
class SourceMask:
    snapshot: int
    mask_path: Path
    slice_path: Path
    mask: np.ndarray
    seed_indices_iz_ix: np.ndarray
    seed_xyz_code: np.ndarray
    cached_rboundary: np.ndarray
    x_code: np.ndarray
    z_code: np.ndarray
    x_mm: np.ndarray
    height_mm: np.ndarray
    logq: np.ndarray
    slice_rboundary: np.ndarray
    q_delta_code: float


def source_paths(mask_dir: Path, slice_dir: Path, snapshot: int) -> tuple[Path, Path]:
    mask_path = mask_dir / f"y0_hybrid_mask_logq4_t{snapshot:04d}.npz"
    slice_path = slice_dir / (
        f"fixed_y_qsl_t{snapshot:04d}_yp0p000000_delta0.25.npz"
    )
    return mask_path, slice_path


def load_source_mask(mask_dir: Path, slice_dir: Path, snapshot: int) -> SourceMask:
    mask_path, slice_path = source_paths(mask_dir, slice_dir, snapshot)
    if not mask_path.is_file():
        raise FileNotFoundError(mask_path)
    if not slice_path.is_file():
        raise FileNotFoundError(slice_path)

    with np.load(mask_path, allow_pickle=False) as archive:
        required = {
            "mask",
            "selected_seed_xyz_amrvac",
            "selected_seed_grid_indices_iz_ix",
            "selected_seed_rboundary",
        }
        missing = required.difference(archive.files)
        if missing:
            raise KeyError(f"{mask_path} is missing {sorted(missing)}")
        mask = np.asarray(archive["mask"], dtype=bool)
        seed_xyz = np.asarray(archive["selected_seed_xyz_amrvac"], dtype=np.float32)
        seed_indices = np.asarray(
            archive["selected_seed_grid_indices_iz_ix"], dtype=np.int32
        )
        cached_rb = np.asarray(archive["selected_seed_rboundary"], dtype=np.int16)

    with np.load(slice_path, allow_pickle=False) as archive:
        required = {"seed", "logq", "rboundary", "x_mm", "height_mm", "delta"}
        missing = required.difference(archive.files)
        if missing:
            raise KeyError(f"{slice_path} is missing {sorted(missing)}")
        plane_seed = np.asarray(archive["seed"], dtype=np.float32)
        logq = np.asarray(archive["logq"], dtype=float)
        rb = np.asarray(archive["rboundary"], dtype=np.int16)
        x_mm_2d = np.asarray(archive["x_mm"], dtype=float)
        h_mm_2d = np.asarray(archive["height_mm"], dtype=float)
        q_delta = float(archive["delta"])

    if mask.shape != logq.shape or rb.shape != mask.shape:
        raise ValueError(f"Inconsistent source shapes at t={snapshot}")
    if plane_seed.shape != mask.shape + (3,):
        raise ValueError(f"Unexpected QSL seed shape at t={snapshot}: {plane_seed.shape}")
    if seed_indices.shape != (int(mask.sum()), 2):
        raise ValueError(f"Unexpected selected-index shape at t={snapshot}")
    if len(seed_xyz) != len(seed_indices) or len(cached_rb) != len(seed_indices):
        raise ValueError(f"Selected arrays have different lengths at t={snapshot}")
    if len(np.unique(seed_indices, axis=0)) != len(seed_indices):
        raise ValueError(f"Duplicate selected mask indices at t={snapshot}")
    if not np.all(mask[seed_indices[:, 0], seed_indices[:, 1]]):
        raise ValueError(f"Selected indices do not lie in the source mask at t={snapshot}")
    if not np.array_equal(
        np.argwhere(mask).astype(np.int32), seed_indices.astype(np.int32)
    ):
        raise ValueError(f"Selected-index ordering differs from np.argwhere at t={snapshot}")
    expected_seed = plane_seed[seed_indices[:, 0], seed_indices[:, 1]]
    if not np.allclose(expected_seed, seed_xyz, atol=1.0e-7, rtol=0.0):
        raise ValueError(f"Saved selected seeds disagree with the QSL plane at t={snapshot}")
    if not np.allclose(seed_xyz[:, 1], FIXED_Y_CODE, atol=1.0e-7):
        raise ValueError(f"Selected seeds are not on exact y=0 at t={snapshot}")

    x_code = plane_seed[0, :, 0].astype(float)
    z_code = plane_seed[:, 0, 2].astype(float)
    x_mm = x_mm_2d[0].astype(float)
    height_mm = h_mm_2d[:, 0].astype(float)
    if not np.allclose(x_mm, x_code * LENGTH_UNIT_MM, atol=1.0e-5):
        raise ValueError(f"x coordinate conversion mismatch at t={snapshot}")
    if not np.allclose(
        height_mm,
        (z_code - BOTTOM_DATUM_Z_CODE) * LENGTH_UNIT_MM,
        atol=2.0e-5,
    ):
        raise ValueError(f"height coordinate conversion mismatch at t={snapshot}")
    return SourceMask(
        snapshot=snapshot,
        mask_path=mask_path,
        slice_path=slice_path,
        mask=mask,
        seed_indices_iz_ix=seed_indices,
        seed_xyz_code=seed_xyz,
        cached_rboundary=cached_rb,
        x_code=x_code,
        z_code=z_code,
        x_mm=x_mm,
        height_mm=height_mm,
        logq=logq,
        slice_rboundary=rb,
        q_delta_code=q_delta,
    )


def coordinate_pixel_edges(values: np.ndarray, *, clip_to_centres: bool) -> np.ndarray:
    """Return pixel edges for a strictly increasing coordinate vector."""

    values = np.asarray(values, dtype=float)
    if values.ndim != 1 or len(values) < 2 or not np.all(np.diff(values) > 0.0):
        raise ValueError("Coordinate vector must be 1-D and strictly increasing")
    middle = 0.5 * (values[:-1] + values[1:])
    if clip_to_centres:
        first, last = values[0], values[-1]
    else:
        first = values[0] - 0.5 * (values[1] - values[0])
        last = values[-1] + 0.5 * (values[-1] - values[-2])
    return np.concatenate(([first], middle, [last]))


def source_pixel_areas_cm2(source: SourceMask) -> np.ndarray:
    """Pixel areas, clipped at the FastQSL cell-centre domain boundary."""

    x_edges = coordinate_pixel_edges(source.x_code, clip_to_centres=True)
    z_edges = coordinate_pixel_edges(source.z_code, clip_to_centres=True)
    widths = np.diff(x_edges) * LENGTH_UNIT_MM * 1.0e8
    heights = np.diff(z_edges) * LENGTH_UNIT_MM * 1.0e8
    area = heights[:, None] * widths[None, :]
    return area[source.mask]


def trace_fastqsl(
    bx: np.ndarray,
    by: np.ndarray,
    bz: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    seeds: np.ndarray,
    *,
    threads: int,
    endpoints: bool = False,
    paths: bool = False,
) -> dict[str, Any]:
    if seeds.ndim != 2 or seeds.shape[1] != 3 or not len(seeds):
        raise ValueError(f"Invalid FastQSL seed shape: {seeds.shape}")
    with tempfile.TemporaryDirectory(prefix="fastqsl-y0-footpoints-") as tmp_dir:
        return fastqsl(
            bx,
            by,
            bz,
            xa=x,
            ya=y,
            za=z,
            seed=np.asarray(seeds, dtype=np.float32),
            B_out=True,
            rF_out=endpoints,
            targetB_out=endpoints,
            path_out=paths,
            nthreads=threads,
            silent=True,
            tmp_dir=tmp_dir,
        )


def endpoint_polarities(result: dict[str, Any], bottom_z: float) -> dict[str, np.ndarray]:
    rb = np.asarray(result["rboundary"], dtype=np.int16)
    rfs = np.asarray(result["rFs"], dtype=float)
    rfe = np.asarray(result["rFe"], dtype=float)
    bs = np.asarray(result["Bs"], dtype=float)
    be = np.asarray(result["Be"], dtype=float)
    if any(array.shape != rfs.shape for array in (rfe, bs, be)):
        raise ValueError("Endpoint arrays have inconsistent shapes")
    if rfs.shape != (len(rb), 3):
        raise ValueError(f"Unexpected endpoint shape {rfs.shape}")

    both_bottom = (
        (rb == 11)
        & np.isclose(rfs[:, 2], bottom_z, atol=2.0e-5)
        & np.isclose(rfe[:, 2], bottom_z, atol=2.0e-5)
    )
    opposite = bs[:, 2] * be[:, 2] < 0.0
    valid = both_bottom & opposite
    positive = np.full_like(rfs, np.nan)
    negative = np.full_like(rfs, np.nan)
    positive_b = np.full_like(bs, np.nan)
    negative_b = np.full_like(bs, np.nan)
    choose_s_positive = valid & (bs[:, 2] > 0.0)
    choose_e_positive = valid & (be[:, 2] > 0.0)
    positive[choose_s_positive] = rfs[choose_s_positive]
    positive[choose_e_positive] = rfe[choose_e_positive]
    positive_b[choose_s_positive] = bs[choose_s_positive]
    positive_b[choose_e_positive] = be[choose_e_positive]
    negative[choose_s_positive] = rfe[choose_s_positive]
    negative[choose_e_positive] = rfs[choose_e_positive]
    negative_b[choose_s_positive] = be[choose_s_positive]
    negative_b[choose_e_positive] = bs[choose_e_positive]
    return {
        "rboundary": rb,
        "rFs": rfs,
        "rFe": rfe,
        "Bs": bs,
        "Be": be,
        "valid": valid,
        "both_bottom": both_bottom,
        "opposite": opposite,
        "positive_endpoint": positive,
        "negative_endpoint": negative,
        "positive_endpoint_B": positive_b,
        "negative_endpoint_B": negative_b,
    }


def flux_components(field_g: np.ndarray, area_cm2: np.ndarray) -> dict[str, float]:
    field_g = np.asarray(field_g, dtype=float)
    area_cm2 = np.asarray(area_cm2, dtype=float)
    if field_g.shape != area_cm2.shape:
        raise ValueError("Field and area arrays must have identical shapes")
    if not np.all(np.isfinite(field_g)) or not np.all(np.isfinite(area_cm2)):
        raise ValueError("Flux inputs contain non-finite values")
    weighted = field_g * area_cm2
    positive = float(weighted[field_g > 0.0].sum())
    negative = float(weighted[field_g < 0.0].sum())
    return {
        "signed_Mx": float(weighted.sum()),
        "positive_Mx": positive,
        "negative_signed_Mx": negative,
        "unsigned_Mx": positive + abs(negative),
    }


def points_in_source_mask(
    x_cross: np.ndarray,
    z_cross: np.ndarray,
    source: SourceMask,
) -> np.ndarray:
    x_edges = coordinate_pixel_edges(source.x_code, clip_to_centres=True)
    z_edges = coordinate_pixel_edges(source.z_code, clip_to_centres=True)
    ix = np.searchsorted(x_edges, x_cross, side="right") - 1
    iz = np.searchsorted(z_edges, z_cross, side="right") - 1
    # Include an intersection exactly on the upper/right outer edge.
    ix[np.isclose(x_cross, x_edges[-1], atol=1.0e-8)] = len(source.x_code) - 1
    iz[np.isclose(z_cross, z_edges[-1], atol=1.0e-8)] = len(source.z_code) - 1
    inside_domain = (
        (iz >= 0)
        & (iz < source.mask.shape[0])
        & (ix >= 0)
        & (ix < source.mask.shape[1])
    )
    selected = np.zeros(len(x_cross), dtype=bool)
    selected[inside_domain] = source.mask[iz[inside_domain], ix[inside_domain]]
    return selected


def path_crossing_count(path: np.ndarray, source: SourceMask, y_tolerance: float) -> int:
    """Count distinct y=0 crossing events that land inside the source mask."""

    path = np.asarray(path, dtype=float)
    if path.ndim != 2 or path.shape[1] != 3 or not len(path):
        return 0
    y = path[:, 1]
    candidates: list[np.ndarray] = []

    zero = np.abs(y) <= y_tolerance
    if np.any(zero):
        indices = np.flatnonzero(zero)
        breaks = np.flatnonzero(np.diff(indices) > 1) + 1
        for group in np.split(indices, breaks):
            candidates.append(path[int(group[len(group) // 2])])

    strict_cross = ((y[:-1] < -y_tolerance) & (y[1:] > y_tolerance)) | (
        (y[:-1] > y_tolerance) & (y[1:] < -y_tolerance)
    )
    for index in np.flatnonzero(strict_cross):
        fraction = -y[index] / (y[index + 1] - y[index])
        candidates.append(path[index] + fraction * (path[index + 1] - path[index]))

    if not candidates:
        return 0
    points = np.asarray(candidates, dtype=float)
    # A zero-valued integration vertex and an adjacent interpolated segment can
    # represent the same crossing; de-duplicate well below one QSL pixel.
    quantization = max(source.q_delta_code * 1.0e-4, 1.0e-8)
    keys = np.rint(points[:, [0, 2]] / quantization).astype(np.int64)
    _, unique_indices = np.unique(keys, axis=0, return_index=True)
    points = points[np.sort(unique_indices)]
    return int(points_in_source_mask(points[:, 0], points[:, 2], source).sum())


def bottom_sampling_grid(
    x: np.ndarray,
    y: np.ndarray,
    z_bottom: float,
    native_spacing: float,
    subsample: int,
    trace_limits_mm: tuple[float, float, float, float],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    if subsample <= 0:
        raise ValueError("bottom subsample must be positive")
    spacing = native_spacing / subsample
    x_all = np.arange(x[0], x[-1] + 0.5 * spacing, spacing)
    y_all = np.arange(y[0], y[-1] + 0.5 * spacing, spacing)
    x_mm = x_all * LENGTH_UNIT_MM
    y_mm = y_all * LENGTH_UNIT_MM
    keep_x = (x_mm >= trace_limits_mm[0]) & (x_mm <= trace_limits_mm[1])
    keep_y = (y_mm >= trace_limits_mm[2]) & (y_mm <= trace_limits_mm[3])
    x_selected = x_all[keep_x]
    y_selected = y_all[keep_y]
    if len(x_selected) < 3 or len(y_selected) < 3:
        raise ValueError("Bottom trace limits select too few grid points")
    yy, xx = np.meshgrid(y_selected, x_selected, indexing="ij")
    seeds = np.column_stack(
        (xx.ravel(), yy.ravel(), np.full(xx.size, z_bottom, dtype=float))
    ).astype(np.float32)
    return x_selected, y_selected, seeds, spacing


def trace_bottom_footprint(
    source: SourceMask,
    bx: np.ndarray,
    by: np.ndarray,
    bz: np.ndarray,
    x: np.ndarray,
    y: np.ndarray,
    z: np.ndarray,
    *,
    threads: int,
    subsample: int,
    trace_limits_mm: tuple[float, float, float, float],
) -> dict[str, np.ndarray | float | int]:
    native_spacing = float(x[1] - x[0])
    bottom_x, bottom_y, seeds, spacing = bottom_sampling_grid(
        x, y, float(z[0]), native_spacing, subsample, trace_limits_mm
    )
    result = trace_fastqsl(
        bx, by, bz, x, y, z, seeds, threads=threads, endpoints=True, paths=True
    )
    paths = result.get("path")
    if paths is None or len(paths) != len(seeds):
        raise ValueError("FastQSL path output does not match bottom seeds")
    crossing_count = np.fromiter(
        (
            path_crossing_count(path, source, y_tolerance=1.0e-7)
            for path in paths
        ),
        dtype=np.uint16,
        count=len(paths),
    )
    shape = (len(bottom_y), len(bottom_x))
    crossing_count = crossing_count.reshape(shape)
    rb = np.asarray(result["rboundary"], dtype=np.int16).reshape(shape)
    bottom_b = np.asarray(result["B"], dtype=float).reshape(shape + (3,))
    endpoint_status = endpoint_polarities(result, float(z[0]))
    valid_pair = np.asarray(endpoint_status["valid"], dtype=bool).reshape(shape)
    hit = crossing_count > 0
    paired = hit & valid_pair
    positive = paired & (bottom_b[..., 2] > 0.0)
    negative = paired & (bottom_b[..., 2] < 0.0)
    unpaired = hit & ~valid_pair
    ambiguous_rboundary11 = hit & (rb == 11) & ~valid_pair
    ambiguous = paired & ~(positive | negative)
    if np.any(ambiguous):
        raise ValueError(
            f"{int(ambiguous.sum())} paired bottom seeds have zero/invalid Bz at "
            f"t={source.snapshot}"
        )
    # The selected region must be safely interior to the quantitative trace
    # window.  This protects against a silently clipped footpoint component.
    for name, selected in (("all-hit", hit), ("paired", paired)):
        edge_count = int(
            selected[0].sum()
            + selected[-1].sum()
            + selected[:, 0].sum()
            + selected[:, -1].sum()
        )
        if edge_count:
            raise RuntimeError(
                f"{name} footprint touches the bottom trace-window edge at "
                f"t={source.snapshot} ({edge_count} samples); enlarge --trace-limits"
            )
    return {
        "x_code": bottom_x,
        "y_code": bottom_y,
        "spacing_code": float(spacing),
        "B_code": bottom_b,
        "rboundary": rb,
        "crossing_count": crossing_count,
        "hit": hit,
        "paired": paired,
        "positive": positive,
        "negative": negative,
        "unpaired": unpaired,
        "ambiguous_rboundary11": ambiguous_rboundary11,
    }


def endpoint_range_check(
    endpoints: dict[str, np.ndarray],
    trace_limits_mm: tuple[float, float, float, float],
    snapshot: int,
    margin_mm: float = 0.2,
) -> None:
    valid = np.asarray(endpoints["valid"], dtype=bool)
    if not np.any(valid):
        raise ValueError(f"No valid paired endpoints at t={snapshot}")
    points = np.vstack(
        (
            np.asarray(endpoints["positive_endpoint"])[valid, :2],
            np.asarray(endpoints["negative_endpoint"])[valid, :2],
        )
    ) * LENGTH_UNIT_MM
    limits = trace_limits_mm
    if (
        points[:, 0].min() <= limits[0] + margin_mm
        or points[:, 0].max() >= limits[1] - margin_mm
        or points[:, 1].min() <= limits[2] + margin_mm
        or points[:, 1].max() >= limits[3] - margin_mm
    ):
        raise RuntimeError(
            f"Forward endpoint cloud approaches --trace-limits at t={snapshot}: "
            f"x=[{points[:,0].min():.2f},{points[:,0].max():.2f}] Mm, "
            f"y=[{points[:,1].min():.2f},{points[:,1].max():.2f}] Mm"
        )


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def prefixed_flux(prefix: str, values: dict[str, float]) -> dict[str, float]:
    return {f"{prefix}_{key}": value for key, value in values.items()}


def process_snapshot(
    snapshot: int,
    *,
    magnetic_dir: Path,
    slice_dir: Path,
    mask_dir: Path,
    per_frame_dir: Path,
    threads: int,
    bottom_subsample: int,
    trace_limits_mm: tuple[float, float, float, float],
    b_unit_g: float,
) -> tuple[dict[str, Any], dict[str, np.ndarray]]:
    source = load_source_mask(mask_dir, slice_dir, snapshot)
    magnetic_path = magnetic_dir / f"nlfff_{snapshot:04d}.binary"
    if not magnetic_path.is_file():
        raise FileNotFoundError(magnetic_path)
    bx, by, bz = read_binary_data(magnetic_path)
    nz, ny, nx = bx.shape
    x = cell_centres(DOMAIN["x"], nx)
    y = cell_centres(DOMAIN["y"], ny)
    z = cell_centres(DOMAIN["z"], nz)
    native_spacing = float(x[1] - x[0])
    if not np.allclose(np.diff(y), native_spacing) or not np.allclose(
        np.diff(z), native_spacing
    ):
        raise ValueError(f"Magnetic grid is not uniform at t={snapshot}")

    forward = trace_fastqsl(
        bx,
        by,
        bz,
        x,
        y,
        z,
        source.seed_xyz_code,
        threads=threads,
        endpoints=True,
    )
    source_b_code = np.asarray(forward["B"], dtype=float)
    if source_b_code.shape != (len(source.seed_xyz_code), 3):
        raise ValueError(f"Unexpected source B shape at t={snapshot}")
    endpoints = endpoint_polarities(forward, float(z[0]))
    rb_mismatch = int(
        np.count_nonzero(endpoints["rboundary"] != source.cached_rboundary)
    )
    mismatch_limit = max(2, int(math.ceil(0.001 * len(source.seed_xyz_code))))
    if rb_mismatch > mismatch_limit:
        raise RuntimeError(
            f"Fresh tracing disagrees with {rb_mismatch} cached rboundary values "
            f"at t={snapshot}; limit={mismatch_limit}"
        )
    endpoint_range_check(endpoints, trace_limits_mm, snapshot)

    areas = source_pixel_areas_cm2(source)
    by_g = source_b_code[:, 1] * b_unit_g
    full_flux = flux_components(by_g, areas)
    closed = np.asarray(endpoints["valid"], dtype=bool)
    closed_flux = flux_components(by_g[closed], areas[closed])
    nonclosed_flux = flux_components(by_g[~closed], areas[~closed])

    bottom = trace_bottom_footprint(
        source,
        bx,
        by,
        bz,
        x,
        y,
        z,
        threads=threads,
        subsample=bottom_subsample,
        trace_limits_mm=trace_limits_mm,
    )
    bottom_bz_g = np.asarray(bottom["B_code"])[..., 2] * b_unit_g
    bottom_area_cm2 = (
        float(bottom["spacing_code"]) * LENGTH_UNIT_MM * 1.0e8
    ) ** 2
    positive_mask = np.asarray(bottom["positive"], dtype=bool)
    negative_mask = np.asarray(bottom["negative"], dtype=bool)
    phi_bz_positive = float(bottom_bz_g[positive_mask].sum() * bottom_area_cm2)
    phi_bz_negative = float(bottom_bz_g[negative_mask].sum() * bottom_area_cm2)
    phi_bz_unsigned = phi_bz_positive + abs(phi_bz_negative)
    phi_bz_pair = 0.5 * phi_bz_unsigned
    balance_denominator = phi_bz_unsigned
    balance_error = (
        abs(phi_bz_positive - abs(phi_bz_negative)) / balance_denominator
        if balance_denominator > 0.0
        else float("nan")
    )
    cross_reference = closed_flux["unsigned_Mx"]
    cross_bottom_difference = (
        (phi_bz_pair - cross_reference) / cross_reference
        if cross_reference > 0.0
        else float("nan")
    )

    rb_hist = Counter(int(value) for value in endpoints["rboundary"])
    bottom_rb_hist = Counter(
        int(value) for value in np.asarray(bottom["rboundary"])[bottom["hit"]]
    )
    crossing = np.asarray(bottom["crossing_count"])
    record: dict[str, Any] = {
        "snapshot": snapshot,
        "time_minutes": snapshot * TIME_MINUTES,
        "source_mask_pixel_count": int(source.mask.sum()),
        "source_mask_area_Mm2": float(areas.sum() / 1.0e16),
        "fresh_closed_seed_count": int(closed.sum()),
        "fresh_nonclosed_seed_count": int((~closed).sum()),
        "fresh_rboundary11_seed_count": int(
            np.count_nonzero(np.asarray(endpoints["rboundary"]) == 11)
        ),
        "fresh_rboundary11_ambiguous_polarity_count": int(
            np.count_nonzero(endpoints["both_bottom"] & ~endpoints["opposite"])
        ),
        "cached_fresh_rboundary_mismatch_count": rb_mismatch,
        "source_rboundary_histogram_json": json.dumps(
            dict(sorted(rb_hist.items())), separators=(",", ":")
        ),
        **prefixed_flux("phi_y_full", full_flux),
        **prefixed_flux("phi_y_closed", closed_flux),
        **prefixed_flux("phi_y_nonclosed", nonclosed_flux),
        "bottom_subsample_per_native_cell": bottom_subsample,
        "bottom_spacing_Mm": float(bottom["spacing_code"]) * LENGTH_UNIT_MM,
        "bottom_hit_pixel_count": int(np.asarray(bottom["hit"]).sum()),
        "bottom_paired_pixel_count": int(np.asarray(bottom["paired"]).sum()),
        "bottom_positive_pixel_count": int(positive_mask.sum()),
        "bottom_negative_pixel_count": int(negative_mask.sum()),
        "bottom_unpaired_pixel_count": int(np.asarray(bottom["unpaired"]).sum()),
        "bottom_rboundary11_ambiguous_polarity_pixel_count": int(
            np.asarray(bottom["ambiguous_rboundary11"]).sum()
        ),
        "bottom_multi_crossing_pixel_count": int(
            np.count_nonzero(np.asarray(bottom["paired"]) & (crossing > 1))
        ),
        "bottom_max_crossing_count": int(crossing.max(initial=0)),
        "bottom_hit_rboundary_histogram_json": json.dumps(
            dict(sorted(bottom_rb_hist.items())), separators=(",", ":")
        ),
        "phi_bz_positive_Mx": phi_bz_positive,
        "phi_bz_negative_signed_Mx": phi_bz_negative,
        "phi_bz_unsigned_Mx": phi_bz_unsigned,
        "phi_bz_paired_Mx": phi_bz_pair,
        "phi_bz_polarity_imbalance_fraction": balance_error,
        "phi_bz_pair_vs_phi_y_closed_unsigned_fraction": cross_bottom_difference,
        "magnetic_binary": str(magnetic_path.resolve()),
        "source_mask": str(source.mask_path.resolve()),
        "source_qsl_slice": str(source.slice_path.resolve()),
    }

    output_path = per_frame_dir / f"traced_y0_mask_footpoints_t{snapshot:04d}.npz"
    np.savez_compressed(
        output_path,
        schema_version=np.asarray(SCHEMA_VERSION, dtype=np.int16),
        snapshot=np.asarray(snapshot, dtype=np.int32),
        time_minutes=np.asarray(record["time_minutes"]),
        source_mask=source.mask.astype(np.uint8),
        source_x_Mm=source.x_mm.astype(np.float32),
        source_height_Mm=source.height_mm.astype(np.float32),
        source_log10Q=source.logq.astype(np.float32),
        source_slice_rboundary=source.slice_rboundary.astype(np.int16),
        source_seed_grid_indices_iz_ix=source.seed_indices_iz_ix,
        source_seed_xyz_amrvac=source.seed_xyz_code,
        source_seed_cell_area_cm2=areas,
        source_B_code=source_b_code.astype(np.float32),
        source_cached_rboundary=source.cached_rboundary.astype(np.int16),
        source_fresh_rboundary=np.asarray(endpoints["rboundary"], dtype=np.int16),
        source_rFs_amrvac=np.asarray(endpoints["rFs"], dtype=np.float32),
        source_rFe_amrvac=np.asarray(endpoints["rFe"], dtype=np.float32),
        source_Bs_code=np.asarray(endpoints["Bs"], dtype=np.float32),
        source_Be_code=np.asarray(endpoints["Be"], dtype=np.float32),
        source_valid_bottom_to_bottom=np.asarray(endpoints["valid"], dtype=np.uint8),
        positive_endpoint_xyz_amrvac=np.asarray(
            endpoints["positive_endpoint"], dtype=np.float32
        ),
        negative_endpoint_xyz_amrvac=np.asarray(
            endpoints["negative_endpoint"], dtype=np.float32
        ),
        positive_endpoint_xy_Mm=(
            np.asarray(endpoints["positive_endpoint"])[:, :2] * LENGTH_UNIT_MM
        ).astype(np.float32),
        negative_endpoint_xy_Mm=(
            np.asarray(endpoints["negative_endpoint"])[:, :2] * LENGTH_UNIT_MM
        ).astype(np.float32),
        bottom_x_Mm=(np.asarray(bottom["x_code"]) * LENGTH_UNIT_MM).astype(np.float32),
        bottom_y_Mm=(np.asarray(bottom["y_code"]) * LENGTH_UNIT_MM).astype(np.float32),
        bottom_Bz_code=np.asarray(bottom["B_code"], dtype=np.float32)[..., 2],
        bottom_trace_rboundary=np.asarray(bottom["rboundary"], dtype=np.int16),
        bottom_crossing_count=np.asarray(bottom["crossing_count"], dtype=np.uint16),
        bottom_hit_mask=np.asarray(bottom["hit"], dtype=np.uint8),
        bottom_paired_mask=np.asarray(bottom["paired"], dtype=np.uint8),
        positive_footpoint_mask=positive_mask.astype(np.uint8),
        negative_footpoint_mask=negative_mask.astype(np.uint8),
        unpaired_bottom_footpoint_mask=np.asarray(bottom["unpaired"], dtype=np.uint8),
        ambiguous_rboundary11_bottom_mask=np.asarray(
            bottom["ambiguous_rboundary11"], dtype=np.uint8
        ),
        B_unit_G=np.asarray(b_unit_g),
        length_unit_Mm=np.asarray(LENGTH_UNIT_MM),
        bottom_datum_z_amrvac=np.asarray(BOTTOM_DATUM_Z_CODE),
        terminal_surface_z_amrvac=np.asarray(float(z[0])),
        qsl_mask_spacing_amrvac=np.asarray(source.q_delta_code),
        bottom_sampling_spacing_amrvac=np.asarray(bottom["spacing_code"]),
        bottom_trace_limits_Mm=np.asarray(trace_limits_mm),
        magnetic_binary=np.asarray(str(magnetic_path.resolve())),
        source_mask_file=np.asarray(str(source.mask_path.resolve())),
        source_qsl_slice_file=np.asarray(str(source.slice_path.resolve())),
        quality_json=np.asarray(json.dumps(record, sort_keys=True)),
        footprint_definition=np.asarray(
            "unique bottom samples whose rboundary==11 field line crosses the solid y=0 mask at least once"
        ),
    )
    record["per_frame_npz"] = str(output_path.resolve())
    arrays = {
        "source_mask": source.mask,
        "source_x_Mm": source.x_mm,
        "source_height_Mm": source.height_mm,
        "source_logq": source.logq,
        "source_rboundary": source.slice_rboundary,
        "bottom_x_Mm": np.asarray(bottom["x_code"]) * LENGTH_UNIT_MM,
        "bottom_y_Mm": np.asarray(bottom["y_code"]) * LENGTH_UNIT_MM,
        "bottom_Bz_code": np.asarray(bottom["B_code"])[..., 2],
        "positive_mask": positive_mask,
        "negative_mask": negative_mask,
        "unpaired_mask": np.asarray(bottom["unpaired"], dtype=bool),
        "crossing_count": crossing,
    }
    return record, arrays


def outer_closed_component(rboundary: np.ndarray, source_mask: np.ndarray) -> np.ndarray:
    closed = np.asarray(rboundary) == 11
    labels, count = label(closed)
    if count == 0:
        return closed
    overlap = np.bincount(labels[source_mask].ravel(), minlength=count + 1)
    overlap[0] = 0
    component = labels == int(np.argmax(overlap))
    return binary_fill_holes(component)


def place_sign(
    axis: plt.Axes,
    mask: np.ndarray,
    x_mm: np.ndarray,
    y_mm: np.ndarray,
    sign: str,
) -> None:
    if not np.any(mask):
        return
    dy = float(np.median(np.diff(y_mm)))
    dx = float(np.median(np.diff(x_mm)))
    clearance = distance_transform_edt(mask, sampling=(dy, dx))
    iy, ix = np.unravel_index(int(np.argmax(clearance)), clearance.shape)
    axis.text(
        x_mm[ix],
        y_mm[iy],
        rf"${sign}$",
        ha="center",
        va="center",
        fontsize=25,
        fontweight="bold",
        color="#202020",
        zorder=7,
    )


def render_frame(
    record: dict[str, Any],
    arrays: dict[str, np.ndarray],
    *,
    bottom_view_limits_mm: tuple[float, float, float, float],
    y0_view_limits_mm: tuple[float, float, float, float],
    dpi: int,
) -> Image.Image:
    figure, axes = plt.subplots(1, 2, figsize=(13.2, 5.9), constrained_layout=True)
    figure.patch.set_facecolor("white")
    for axis in axes:
        axis.set_facecolor("white")

    source_mask = arrays["source_mask"]
    x0 = arrays["source_x_Mm"]
    h0 = arrays["source_height_Mm"]
    layer = np.ma.masked_where(~source_mask, np.ones(source_mask.shape))
    axes[0].pcolormesh(
        x0,
        h0,
        layer,
        shading="nearest",
        cmap=ListedColormap([MASK_COLOR]),
        vmin=0.0,
        vmax=1.0,
        rasterized=True,
        zorder=1,
    )
    axes[0].contour(
        x0,
        h0,
        source_mask.astype(float),
        levels=[0.5],
        colors=[MASK_EDGE_COLOR],
        linewidths=1.1,
        zorder=2,
    )
    outer_closed = outer_closed_component(arrays["source_rboundary"], source_mask)
    axes[0].contour(
        x0,
        h0,
        outer_closed.astype(float),
        levels=[0.5],
        colors=[CLOSED_COLOR],
        linewidths=2.0,
        zorder=4,
    )
    finite = arrays["source_logq"][np.isfinite(arrays["source_logq"])]
    if finite.size and finite.min() <= 4.0 <= finite.max():
        axes[0].contour(
            x0,
            h0,
            arrays["source_logq"],
            levels=[4.0],
            colors=["white"],
            linewidths=3.8,
            zorder=3,
        )
        axes[0].contour(
            x0,
            h0,
            arrays["source_logq"],
            levels=[4.0],
            colors=[Q4_COLOR],
            linewidths=2.0,
            zorder=3.5,
        )
    axes[0].set_xlim(y0_view_limits_mm[:2])
    axes[0].set_ylim(y0_view_limits_mm[2:])
    axes[0].set_xlabel("Solar-X [Mm]", fontsize=13)
    axes[0].set_ylabel("Height above bottom boundary [Mm]", fontsize=13)
    axes[0].set_aspect("equal")
    axes[0].text(
        0.02,
        0.98,
        r"(a) Selected $y=0$ cross-section",
        transform=axes[0].transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        zorder=8,
    )

    xb = arrays["bottom_x_Mm"]
    yb = arrays["bottom_y_Mm"]
    positive = arrays["positive_mask"]
    negative = arrays["negative_mask"]
    category = np.zeros(positive.shape, dtype=float)
    category[positive] = 1.0
    category[negative] = 0.0
    selected = positive | negative
    signed_layer = np.ma.masked_where(~selected, category)
    axes[1].pcolormesh(
        xb,
        yb,
        signed_layer,
        shading="nearest",
        cmap=ListedColormap([NEGATIVE_COLOR, POSITIVE_COLOR]),
        vmin=-0.5,
        vmax=1.5,
        alpha=0.78,
        rasterized=True,
        zorder=2,
    )
    for mask, color in ((positive, "#7A275B"), (negative, "#004A77")):
        if np.any(mask):
            axes[1].contour(
                xb,
                yb,
                mask.astype(float),
                levels=[0.5],
                colors=[color],
                linewidths=0.9,
                zorder=3,
            )
    bz = arrays["bottom_Bz_code"]
    if np.nanmin(bz) <= 0.0 <= np.nanmax(bz):
        axes[1].contour(
            xb,
            yb,
            bz,
            levels=[0.0],
            colors=[BZ0_COLOR],
            linestyles="dashed",
            linewidths=1.35,
            zorder=5,
        )
    place_sign(axes[1], positive, xb, yb, "+")
    place_sign(axes[1], negative, xb, yb, "-")
    axes[1].set_xlim(bottom_view_limits_mm[:2])
    axes[1].set_ylim(bottom_view_limits_mm[2:])
    axes[1].set_xlabel("Solar-X [Mm]", fontsize=13)
    axes[1].set_ylabel("Solar-Y [Mm]", fontsize=13)
    axes[1].set_aspect("equal")
    axes[1].text(
        0.02,
        0.98,
        "(b) Traced bottom footpoints",
        transform=axes[1].transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        zorder=8,
    )
    handles = [
        Patch(facecolor=MASK_COLOR, edgecolor=MASK_EDGE_COLOR, label="Selected cross-section"),
        Patch(facecolor=POSITIVE_COLOR, edgecolor="#7A275B", label=r"Footpoints: $B_z>0$"),
        Line2D([0], [0], color=CLOSED_COLOR, linewidth=2.0, label=r"$r_{\rm boundary}=11$"),
        Patch(facecolor=NEGATIVE_COLOR, edgecolor="#004A77", label=r"Footpoints: $B_z<0$"),
        Line2D([0], [0], color=Q4_COLOR, linewidth=2.0, label=r"$\log_{10}Q=4$"),
        Line2D([0], [0], color=BZ0_COLOR, linestyle="--", linewidth=1.5, label=r"$B_z=0$"),
    ]
    figure.legend(
        handles=handles,
        loc="outside lower center",
        ncol=3,
        frameon=False,
        fontsize=11.5,
        handlelength=2.5,
        columnspacing=1.4,
    )
    figure.suptitle(
        f"t = {record['time_minutes']:.1f} min",
        fontsize=18,
        y=0.995,
    )
    for axis in axes:
        axis.tick_params(direction="in", top=True, right=True, labelsize=11.5)

    buffer = BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi, facecolor="white")
    plt.close(figure)
    buffer.seek(0)
    with Image.open(buffer) as image:
        return image.convert("RGB")


def save_gif(
    images: list[Image.Image],
    path: Path,
    *,
    fps: float,
    endpoint_pause: float,
) -> None:
    if not images:
        raise ValueError("No images supplied for GIF")
    duration = max(1, int(round(1000.0 / fps)))
    durations = [duration] * len(images)
    pause_ms = max(duration, int(round(endpoint_pause * 1000.0)))
    durations[0] = pause_ms
    durations[-1] = pause_ms
    palette = images[0].convert("P", palette=Image.Palette.ADAPTIVE, colors=128)
    frames = [
        image.quantize(palette=palette, dither=Image.Dither.NONE) for image in images[1:]
    ]
    palette.save(
        path,
        save_all=True,
        append_images=frames,
        duration=durations,
        loop=0,
        optimize=False,
        disposal=2,
    )


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("No rows to write")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def nice_bottom_view(
    arrays_by_snapshot: list[dict[str, np.ndarray]],
    trace_limits_mm: tuple[float, float, float, float],
    padding_mm: float = 3.0,
) -> tuple[float, float, float, float]:
    xs: list[np.ndarray] = []
    ys: list[np.ndarray] = []
    for arrays in arrays_by_snapshot:
        selected = arrays["positive_mask"] | arrays["negative_mask"]
        iy, ix = np.nonzero(selected)
        if len(ix):
            xs.append(arrays["bottom_x_Mm"][ix])
            ys.append(arrays["bottom_y_Mm"][iy])
    if not xs:
        raise ValueError("No paired bottom footprint exists in the series")
    x_all = np.concatenate(xs)
    y_all = np.concatenate(ys)
    xmin = 5.0 * math.floor((float(x_all.min()) - padding_mm) / 5.0)
    xmax = 5.0 * math.ceil((float(x_all.max()) + padding_mm) / 5.0)
    ymin = 5.0 * math.floor((float(y_all.min()) - padding_mm) / 5.0)
    ymax = 5.0 * math.ceil((float(y_all.max()) + padding_mm) / 5.0)
    xmin, xmax = max(xmin, trace_limits_mm[0]), min(xmax, trace_limits_mm[1])
    ymin, ymax = max(ymin, trace_limits_mm[2]), min(ymax, trace_limits_mm[3])
    return float(xmin), float(xmax), float(ymin), float(ymax)


def padded_footprint_limits(
    arrays: dict[str, np.ndarray],
    outer_limits_mm: tuple[float, float, float, float],
    padding_mm: float = 5.0,
) -> tuple[float, float, float, float]:
    """Return a convergence-test window safely around a measured footprint."""

    selected = arrays["positive_mask"] | arrays["negative_mask"] | arrays["unpaired_mask"]
    iy, ix = np.nonzero(selected)
    if not len(ix):
        raise ValueError("Cannot construct convergence limits from an empty footprint")
    x = arrays["bottom_x_Mm"][ix]
    y = arrays["bottom_y_Mm"][iy]
    return (
        max(float(x.min()) - padding_mm, outer_limits_mm[0]),
        min(float(x.max()) + padding_mm, outer_limits_mm[1]),
        max(float(y.min()) - padding_mm, outer_limits_mm[2]),
        min(float(y.max()) + padding_mm, outer_limits_mm[3]),
    )


def make_flux_plot(rows: list[dict[str, Any]], path: Path) -> None:
    time = np.asarray([row["time_minutes"] for row in rows])
    phi_y_full = np.asarray([row["phi_y_full_unsigned_Mx"] for row in rows]) / 1.0e19
    phi_y_closed = np.asarray([row["phi_y_closed_unsigned_Mx"] for row in rows]) / 1.0e19
    phi_bz = np.asarray([row["phi_bz_paired_Mx"] for row in rows]) / 1.0e19
    imbalance = 100.0 * np.asarray(
        [row["phi_bz_polarity_imbalance_fraction"] for row in rows]
    )
    difference = 100.0 * np.asarray(
        [row["phi_bz_pair_vs_phi_y_closed_unsigned_fraction"] for row in rows]
    )
    figure, axes = plt.subplots(2, 1, figsize=(10.5, 7.4), sharex=True, constrained_layout=True)
    axes[0].plot(time, phi_y_full, color=MASK_EDGE_COLOR, linewidth=1.8, label=r"$\int_{\rm solid\ mask}|B_y|\,dA$")
    axes[0].plot(time, phi_y_closed, color=CLOSED_COLOR, linewidth=1.8, label=r"$\int_{r_b=11}|B_y|\,dA$")
    axes[0].plot(time, phi_bz, color="#202020", linewidth=2.0, label=r"$(\Phi_+ + |\Phi_-|)/2$")
    axes[0].set_ylabel(r"Magnetic flux [$10^{19}$ Mx]", fontsize=12)
    axes[0].legend(frameon=False, fontsize=10.5, ncol=3)
    axes[0].grid(alpha=0.22)
    axes[1].plot(time, imbalance, color="#D55E00", linewidth=1.6, label="Bottom polarity imbalance")
    axes[1].plot(time, difference, color="#0072B2", linewidth=1.6, label="Bottom pair vs closed y=0")
    axes[1].axhline(0.0, color="#777777", linewidth=0.8)
    axes[1].set_xlabel("Time [min]", fontsize=12)
    axes[1].set_ylabel("Relative difference [%]", fontsize=12)
    axes[1].legend(frameon=False, fontsize=10.5)
    axes[1].grid(alpha=0.22)
    for axis in axes:
        axis.tick_params(direction="in", top=True, right=True)
    figure.savefig(path, dpi=180, facecolor="white")
    plt.close(figure)


def run_convergence(
    snapshots: Iterable[int],
    factors: Iterable[int],
    *,
    magnetic_dir: Path,
    slice_dir: Path,
    mask_dir: Path,
    threads: int,
    trace_limits_mm: tuple[float, float, float, float],
    trace_limits_by_snapshot: dict[int, tuple[float, float, float, float]] | None,
    b_unit_g: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for snapshot in snapshots:
        source = load_source_mask(mask_dir, slice_dir, snapshot)
        magnetic_path = magnetic_dir / f"nlfff_{snapshot:04d}.binary"
        bx, by, bz = read_binary_data(magnetic_path)
        nz, ny, nx = bx.shape
        x = cell_centres(DOMAIN["x"], nx)
        y = cell_centres(DOMAIN["y"], ny)
        z = cell_centres(DOMAIN["z"], nz)
        current_trace_limits = (
            trace_limits_by_snapshot.get(snapshot, trace_limits_mm)
            if trace_limits_by_snapshot is not None
            else trace_limits_mm
        )
        for factor in factors:
            bottom = trace_bottom_footprint(
                source,
                bx,
                by,
                bz,
                x,
                y,
                z,
                threads=threads,
                subsample=factor,
                trace_limits_mm=current_trace_limits,
            )
            bz_g = np.asarray(bottom["B_code"])[..., 2] * b_unit_g
            area = (float(bottom["spacing_code"]) * LENGTH_UNIT_MM * 1.0e8) ** 2
            positive = np.asarray(bottom["positive"], dtype=bool)
            negative = np.asarray(bottom["negative"], dtype=bool)
            phi_pos = float(bz_g[positive].sum() * area)
            phi_neg = float(bz_g[negative].sum() * area)
            rows.append(
                {
                    "snapshot": snapshot,
                    "time_minutes": snapshot * TIME_MINUTES,
                    "subsamples_per_native_cell": factor,
                    "bottom_spacing_Mm": float(bottom["spacing_code"]) * LENGTH_UNIT_MM,
                    "trace_limits_Mm_json": json.dumps(current_trace_limits),
                    "positive_pixel_count": int(positive.sum()),
                    "negative_pixel_count": int(negative.sum()),
                    "phi_bz_positive_Mx": phi_pos,
                    "phi_bz_negative_signed_Mx": phi_neg,
                    "phi_bz_paired_Mx": 0.5 * (phi_pos + abs(phi_neg)),
                    "polarity_imbalance_fraction": abs(phi_pos - abs(phi_neg))
                    / (phi_pos + abs(phi_neg)),
                }
            )
            print(
                f"convergence t={snapshot:04d} factor={factor}: "
                f"Phi_pair={rows[-1]['phi_bz_paired_Mx']:.6e} Mx",
                flush=True,
            )
    return rows


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--magnetic-dir", type=Path, default=MAGNETIC_DIR)
    parser.add_argument("--slice-dir", type=Path, default=SLICE_DIR)
    parser.add_argument("--mask-dir", type=Path, default=MASK_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--start", type=int, default=START_SNAPSHOT)
    parser.add_argument("--stop", type=int, default=STOP_SNAPSHOT)
    parser.add_argument("--step", type=int, default=1)
    parser.add_argument("--threads", type=int, default=THREADS)
    parser.add_argument("--bottom-subsample", type=int, default=BOTTOM_SUBSAMPLE)
    parser.add_argument("--B-unit-G", type=float, default=B_UNIT_G)
    parser.add_argument(
        "--trace-limits",
        type=float,
        nargs=4,
        metavar=("XMIN", "XMAX", "YMIN", "YMAX"),
        default=BOTTOM_TRACE_LIMITS_MM,
    )
    parser.add_argument("--fps", type=float, default=8.0)
    parser.add_argument("--endpoint-pause", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=120)
    parser.add_argument(
        "--convergence-snapshots", type=int, nargs="*", default=(10, 25, 26, 95)
    )
    parser.add_argument(
        "--convergence-factors", type=int, nargs="*", default=(2, 4, 8)
    )
    parser.add_argument("--skip-convergence", action="store_true")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.step <= 0 or args.threads <= 0 or args.bottom_subsample <= 0:
        raise SystemExit("--step, --threads, and --bottom-subsample must be positive")
    if args.stop < args.start:
        raise SystemExit("--stop must not precede --start")
    if not np.isfinite(args.B_unit_G) or args.B_unit_G <= 0.0:
        raise SystemExit("--B-unit-G must be positive and finite")
    trace_limits = tuple(float(value) for value in args.trace_limits)
    if not (trace_limits[0] < trace_limits[1] and trace_limits[2] < trace_limits[3]):
        raise SystemExit("Invalid --trace-limits ordering")
    snapshots = list(range(args.start, args.stop + 1, args.step))
    if snapshots[-1] != args.stop:
        raise SystemExit("--step must land exactly on --stop")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    per_frame_dir = args.output_dir / "data"
    frame_dir = args.output_dir / "frames"
    per_frame_dir.mkdir(parents=True, exist_ok=True)
    frame_dir.mkdir(parents=True, exist_ok=True)

    rows: list[dict[str, Any]] = []
    arrays_by_snapshot: list[dict[str, np.ndarray]] = []
    for snapshot in snapshots:
        record, arrays = process_snapshot(
            snapshot,
            magnetic_dir=args.magnetic_dir,
            slice_dir=args.slice_dir,
            mask_dir=args.mask_dir,
            per_frame_dir=per_frame_dir,
            threads=args.threads,
            bottom_subsample=args.bottom_subsample,
            trace_limits_mm=trace_limits,
            b_unit_g=args.B_unit_G,
        )
        rows.append(record)
        arrays_by_snapshot.append(arrays)
        print(
            f"finished t={snapshot:04d}: n11={record['fresh_closed_seed_count']}, "
            f"Phi_y={record['phi_y_closed_signed_Mx']:.6e}, "
            f"Phi_z_pair={record['phi_bz_paired_Mx']:.6e} Mx",
            flush=True,
        )

    bottom_view = nice_bottom_view(arrays_by_snapshot, trace_limits)
    images: list[Image.Image] = []
    frame_paths: list[str] = []
    for record, arrays in zip(rows, arrays_by_snapshot):
        image = render_frame(
            record,
            arrays,
            bottom_view_limits_mm=bottom_view,
            y0_view_limits_mm=Y0_VIEW_LIMITS_MM,
            dpi=args.dpi,
        )
        frame_path = frame_dir / f"y0_to_bottom_footpoints_t{record['snapshot']:04d}.png"
        image.save(frame_path)
        images.append(image)
        frame_paths.append(str(frame_path.resolve()))
        record["frame_png"] = str(frame_path.resolve())

    tag = f"t{args.start:04d}_t{args.stop:04d}"
    gif_path = args.output_dir / f"y0_to_bottom_footpoints_{tag}.gif"
    save_gif(images, gif_path, fps=args.fps, endpoint_pause=args.endpoint_pause)
    csv_path = args.output_dir / f"y0_mask_bottom_flux_{tag}.csv"
    write_csv(csv_path, rows)
    flux_plot_path = args.output_dir / f"y0_mask_bottom_flux_{tag}.png"
    make_flux_plot(rows, flux_plot_path)

    stack_path = args.output_dir / f"y0_mask_bottom_footpoints_{tag}_stack.npz"
    np.savez_compressed(
        stack_path,
        schema_version=np.asarray(SCHEMA_VERSION, dtype=np.int16),
        snapshots=np.asarray(snapshots, dtype=np.int32),
        time_minutes=np.asarray([row["time_minutes"] for row in rows]),
        source_masks=np.asarray([arrays["source_mask"] for arrays in arrays_by_snapshot], dtype=np.uint8),
        positive_footpoint_masks=np.asarray([arrays["positive_mask"] for arrays in arrays_by_snapshot], dtype=np.uint8),
        negative_footpoint_masks=np.asarray([arrays["negative_mask"] for arrays in arrays_by_snapshot], dtype=np.uint8),
        unpaired_footpoint_masks=np.asarray([arrays["unpaired_mask"] for arrays in arrays_by_snapshot], dtype=np.uint8),
        bottom_crossing_counts=np.asarray([arrays["crossing_count"] for arrays in arrays_by_snapshot], dtype=np.uint16),
        bottom_Bz_code=np.asarray([arrays["bottom_Bz_code"] for arrays in arrays_by_snapshot], dtype=np.float32),
        source_x_Mm=arrays_by_snapshot[0]["source_x_Mm"].astype(np.float32),
        source_height_Mm=arrays_by_snapshot[0]["source_height_Mm"].astype(np.float32),
        bottom_x_Mm=arrays_by_snapshot[0]["bottom_x_Mm"].astype(np.float32),
        bottom_y_Mm=arrays_by_snapshot[0]["bottom_y_Mm"].astype(np.float32),
        phi_y_full_signed_Mx=np.asarray([row["phi_y_full_signed_Mx"] for row in rows]),
        phi_y_full_unsigned_Mx=np.asarray([row["phi_y_full_unsigned_Mx"] for row in rows]),
        phi_y_closed_signed_Mx=np.asarray([row["phi_y_closed_signed_Mx"] for row in rows]),
        phi_y_closed_unsigned_Mx=np.asarray([row["phi_y_closed_unsigned_Mx"] for row in rows]),
        phi_bz_positive_Mx=np.asarray([row["phi_bz_positive_Mx"] for row in rows]),
        phi_bz_negative_signed_Mx=np.asarray([row["phi_bz_negative_signed_Mx"] for row in rows]),
        phi_bz_paired_Mx=np.asarray([row["phi_bz_paired_Mx"] for row in rows]),
        B_unit_G=np.asarray(args.B_unit_G),
        length_unit_Mm=np.asarray(LENGTH_UNIT_MM),
        bottom_subsample_per_native_cell=np.asarray(args.bottom_subsample),
        bottom_trace_limits_Mm=np.asarray(trace_limits),
        bottom_view_limits_Mm=np.asarray(bottom_view),
        frame_pngs=np.asarray(frame_paths),
    )

    convergence_rows: list[dict[str, Any]] = []
    convergence_csv = None
    if not args.skip_convergence and args.convergence_snapshots:
        invalid = set(args.convergence_snapshots).difference(snapshots)
        if invalid:
            raise ValueError(f"Convergence snapshots are outside this run: {sorted(invalid)}")
        convergence_rows = run_convergence(
            args.convergence_snapshots,
            args.convergence_factors,
            magnetic_dir=args.magnetic_dir,
            slice_dir=args.slice_dir,
            mask_dir=args.mask_dir,
            threads=args.threads,
            trace_limits_mm=trace_limits,
            trace_limits_by_snapshot={
                snapshot: padded_footprint_limits(
                    arrays_by_snapshot[snapshots.index(snapshot)], trace_limits
                )
                for snapshot in args.convergence_snapshots
            },
            b_unit_g=args.B_unit_G,
        )
        convergence_csv = args.output_dir / f"bottom_flux_convergence_{tag}.csv"
        write_csv(convergence_csv, convergence_rows)

    metadata = {
        "schema_version": SCHEMA_VERSION,
        "method": {
            "source_cross_section": (
                "solid tracked y=0 mask; B_y integrated once per selected QSL pixel"
            ),
            "bottom_footprint": (
                "fixed bottom samples traced through the volume; select once when a "
                "path crosses the source mask; require rboundary==11 and split by local Bz"
            ),
            "paired_bottom_flux": "(Phi_Bz_positive + abs(Phi_Bz_negative))/2",
            "nonclosed_filled_pixels": (
                "included in full y=0 flux, excluded from paired bottom footpoints, reported separately"
            ),
        },
        "normalizations": {
            "B_unit_G": args.B_unit_G,
            "length_unit_Mm": LENGTH_UNIT_MM,
            "bottom_datum_z_amrvac": BOTTOM_DATUM_Z_CODE,
            "terminal_surface_note": (
                "FastQSL terminates on the first cell-centre layer, not on the datum z=0.1"
            ),
        },
        "settings": {
            "snapshots": snapshots,
            "threads": args.threads,
            "bottom_subsample_per_native_cell": args.bottom_subsample,
            "bottom_trace_limits_Mm": trace_limits,
            "bottom_view_limits_Mm": bottom_view,
            "y0_view_limits_Mm": Y0_VIEW_LIMITS_MM,
            "fps": args.fps,
            "endpoint_pause_seconds": args.endpoint_pause,
        },
        "quality": {
            "maximum_bottom_polarity_imbalance_fraction": float(
                max(row["phi_bz_polarity_imbalance_fraction"] for row in rows)
            ),
            "maximum_abs_bottom_vs_closed_y0_fraction": float(
                max(
                    abs(row["phi_bz_pair_vs_phi_y_closed_unsigned_fraction"])
                    for row in rows
                )
            ),
            "total_cached_fresh_rboundary_mismatches": int(
                sum(row["cached_fresh_rboundary_mismatch_count"] for row in rows)
            ),
            "frames_with_unpaired_bottom_hits": int(
                sum(row["bottom_unpaired_pixel_count"] > 0 for row in rows)
            ),
            "frames_with_multiple_y0_crossings": int(
                sum(row["bottom_multi_crossing_pixel_count"] > 0 for row in rows)
            ),
        },
        "outputs": {
            "gif": str(gif_path.resolve()),
            "flux_csv": str(csv_path.resolve()),
            "flux_plot": str(flux_plot_path.resolve()),
            "stack_npz": str(stack_path.resolve()),
            "per_frame_directory": str(per_frame_dir.resolve()),
            "frame_directory": str(frame_dir.resolve()),
            "convergence_csv": str(convergence_csv.resolve()) if convergence_csv else None,
        },
        "frames": rows,
        "convergence": convergence_rows,
    }
    metadata_path = args.output_dir / f"y0_mask_bottom_footpoints_{tag}_metadata.json"
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2, ensure_ascii=False)

    print(f"GIF: {gif_path}")
    print(f"Flux CSV: {csv_path}")
    print(f"Flux plot: {flux_plot_path}")
    print(f"Stack: {stack_path}")
    print(f"Metadata: {metadata_path}")


if __name__ == "__main__":
    main()

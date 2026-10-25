#!/usr/bin/env python3
"""Calculate and animate the exact-y=0 current-to-field ratio.

The magnetic inputs are the same uniform ``nlfff_XXXX.binary`` cubes used by
the FastQSL2 workflow.  They contain Bx, By, and Bz in ``(z, y, x)`` order.
For the even-sized y grid, y=0 lies halfway between the two central cell
centres.  The field is linearly interpolated to that plane, while the normal
derivatives are centred differences across the two bracketing layers.

The plotted quantity is

    |J| / |B| = |curl(B)| / |B|,

where MPI-AMRVAC's normalized convention J=curl(B) is used.  The raw ratio is
stored in inverse normalized length and is also converted to Mm^-1 using
L0=10 Mm.  One NPZ is written per snapshot together with a consolidated stack,
per-frame statistics, metadata, all rendered PNG frames, and an animated GIF.

Example
-------
/opt/anaconda3/bin/python \
  /Users/xiaoxingyun/Desktop/solar_code/workflows/Simulation_Analysis/scripts/\
make_y0_j_over_b_evolution.py
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import tempfile
from pathlib import Path
from typing import Any, Iterable

if "MPLCONFIGDIR" not in os.environ:
    mpl_cache = Path(tempfile.gettempdir()) / "solar-code-matplotlib"
    mpl_cache.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(mpl_cache)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
from matplotlib.ticker import LogFormatterMathtext
from PIL import Image


MAGNETIC_DIR = Path("/Volumes/Seagate/w2_jet_mhd/out_path1/nlfff")
OUTPUT_DIR = Path(
    "/Users/xiaoxingyun/Desktop/w2_jet_mhd/qfactor/y0_qsl/j_over_b_evolution"
)
START_SNAPSHOT = 10
STOP_SNAPSHOT = 95
TIME_MINUTES_PER_SNAPSHOT = 0.71558
FIXED_Y_CODE = 0.0
LENGTH_UNIT_MM = 10.0
BOTTOM_Z_CODE = 0.1

# Cell-centred AMRVAC normalized domain used by the QSL products.
DOMAIN = {
    "x": (-4.0, 8.0),
    "y": (-6.0, 6.0),
    "z": (0.1, 12.1),
}
DEFAULT_VIEW_MM = (-15.0, 45.0, 0.0, 60.0)
COLORMAP_NAME = "rainbow"
SCHEMA_VERSION = 1


def cell_centres(bounds: tuple[float, float], count: int) -> np.ndarray:
    spacing = (bounds[1] - bounds[0]) / count
    return bounds[0] + (np.arange(count, dtype=float) + 0.5) * spacing


def finite_stats(values: np.ndarray) -> dict[str, float | int | None]:
    array = np.asarray(values, dtype=float)
    finite = np.isfinite(array)
    selected = array[finite]
    result: dict[str, float | int | None] = {
        "finite_count": int(finite.sum()),
        "nonfinite_count": int((~finite).sum()),
    }
    for name, percentile in (
        ("min", 0.0),
        ("p01", 1.0),
        ("p05", 5.0),
        ("median", 50.0),
        ("p95", 95.0),
        ("p99", 99.0),
        ("p99_7", 99.7),
        ("max", 100.0),
    ):
        result[name] = (
            float(np.percentile(selected, percentile)) if selected.size else None
        )
    return result


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_header(path: Path) -> tuple[int, int, int]:
    with path.open("rb") as stream:
        header = np.fromfile(stream, dtype=np.int32, count=3)
    if header.size != 3:
        raise ValueError(f"Missing (nx, ny, nz) header in {path}")
    nx, ny, nz = (int(value) for value in header)
    if min(nx, ny, nz) < 3:
        raise ValueError(f"Invalid magnetic dimensions in {path}: {header.tolist()}")
    expected_size = 3 * np.dtype(np.int32).itemsize + 3 * nx * ny * nz * 8
    if path.stat().st_size != expected_size:
        raise ValueError(
            f"Unexpected size for {path}: {path.stat().st_size}, expected {expected_size}"
        )
    return nx, ny, nz


def bracketing_y_indices(y: np.ndarray, fixed_y: float) -> tuple[int, int, float]:
    if not (y[0] < fixed_y < y[-1]):
        raise ValueError(f"fixed y={fixed_y} is not inside [{y[0]}, {y[-1]}]")
    upper = int(np.searchsorted(y, fixed_y, side="right"))
    lower = upper - 1
    if lower < 0 or upper >= len(y):
        raise ValueError("Could not bracket the requested y plane")
    weight = float((fixed_y - y[lower]) / (y[upper] - y[lower]))
    if not (0.0 <= weight <= 1.0):
        raise ValueError("Invalid interpolation weight")
    return lower, upper, weight


def read_y_bracketing_slabs(
    path: Path,
    shape_xyz: tuple[int, int, int],
    lower_y: int,
    upper_y: int,
) -> np.ndarray:
    """Read only two y layers from the component-major magnetic binary."""

    nx, ny, nz = shape_xyz
    mapped = np.memmap(
        path,
        dtype=np.float64,
        mode="r",
        offset=3 * np.dtype(np.int32).itemsize,
        shape=(3, nz, ny, nx),
        order="C",
    )
    slabs = np.array(mapped[:, :, [lower_y, upper_y], :], dtype=np.float64)
    del mapped
    if slabs.shape != (3, nz, 2, nx):
        raise ValueError(f"Unexpected slab shape from {path}: {slabs.shape}")
    return slabs


def calculate_exact_y_slice(
    slabs: np.ndarray,
    *,
    x: np.ndarray,
    y_lower: float,
    y_upper: float,
    z: np.ndarray,
    interpolation_weight: float,
    b_floor_code: float,
) -> dict[str, np.ndarray]:
    """Return B, curl(B), and |curl(B)|/|B| on the requested y plane."""

    low = slabs[:, :, 0, :]
    high = slabs[:, :, 1, :]
    b_plane = (1.0 - interpolation_weight) * low + interpolation_weight * high
    db_dy = (high - low) / (y_upper - y_lower)

    bx, by, bz = b_plane
    dbx_dy, _dby_dy, dbz_dy = db_dy
    dby_dz = np.gradient(by, z, axis=0, edge_order=2)
    dbx_dz = np.gradient(bx, z, axis=0, edge_order=2)
    dbz_dx = np.gradient(bz, x, axis=1, edge_order=2)
    dby_dx = np.gradient(by, x, axis=1, edge_order=2)

    # Arrays are indexed (z, x); components remain Cartesian (x, y, z).
    curl_x = dbz_dy - dby_dz
    curl_y = dbx_dz - dbz_dx
    curl_z = dby_dx - dbx_dy
    curl_b = np.stack((curl_x, curl_y, curl_z), axis=-1)
    b_vector = np.moveaxis(b_plane, 0, -1)
    b_magnitude = np.linalg.norm(b_vector, axis=-1)
    curl_b_magnitude = np.linalg.norm(curl_b, axis=-1)

    valid = (
        np.isfinite(b_magnitude)
        & np.isfinite(curl_b_magnitude)
        & (b_magnitude > b_floor_code)
    )
    ratio_per_l0 = np.full(b_magnitude.shape, np.nan, dtype=np.float64)
    ratio_per_l0[valid] = curl_b_magnitude[valid] / b_magnitude[valid]
    ratio_per_mm = ratio_per_l0 / LENGTH_UNIT_MM
    return {
        "B_code": b_vector,
        "curl_B_per_L0": curl_b,
        "B_magnitude_code": b_magnitude,
        "curl_B_magnitude_per_L0": curl_b_magnitude,
        "j_over_b_per_L0": ratio_per_l0,
        "j_over_b_per_Mm": ratio_per_mm,
        "valid": valid,
    }


def display_selection(
    x_mm: np.ndarray,
    height_mm: np.ndarray,
    limits: tuple[float, float, float, float],
) -> tuple[np.ndarray, np.ndarray]:
    xmin, xmax, hmin, hmax = limits
    keep_x = (x_mm >= xmin) & (x_mm <= xmax)
    keep_z = (height_mm >= hmin) & (height_mm <= hmax)
    if int(keep_x.sum()) < 3 or int(keep_z.sum()) < 3:
        raise ValueError(f"Display limits select too few cells: {limits}")
    return keep_x, keep_z


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        raise ValueError("No statistics rows to write")
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def draw_frame(
    path: Path,
    *,
    ratio_per_mm: np.ndarray,
    snapshot: int,
    time_minutes: float,
    x_mm: np.ndarray,
    height_mm: np.ndarray,
    limits: tuple[float, float, float, float],
    vmin: float,
    vmax: float,
    dpi: int,
) -> None:
    cmap = plt.get_cmap(COLORMAP_NAME).copy()
    cmap.set_bad("white")
    figure, axis = plt.subplots(figsize=(9.0, 6.7), constrained_layout=True)
    image = axis.pcolormesh(
        x_mm,
        height_mm,
        ratio_per_mm,
        shading="nearest",
        cmap=cmap,
        norm=LogNorm(vmin=vmin, vmax=vmax, clip=True),
        rasterized=True,
    )
    axis.set_xlim(limits[0], limits[1])
    axis.set_ylim(limits[2], limits[3])
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("x [Mm]", fontsize=15)
    axis.set_ylabel("Height [Mm]", fontsize=15)
    axis.set_title(
        r"$y=0$: $|\nabla\times\mathbf{B}|/|\mathbf{B}|$",
        fontsize=16,
    )
    axis.text(
        0.975,
        0.965,
        f"t = {time_minutes:.1f} min\n(snapshot {snapshot})",
        transform=axis.transAxes,
        ha="right",
        va="top",
        fontsize=13,
        color="white",
        bbox={"facecolor": "black", "edgecolor": "none", "alpha": 0.58, "pad": 4.0},
    )
    axis.tick_params(direction="in", top=True, right=True, labelsize=12)
    colorbar = figure.colorbar(
        image,
        ax=axis,
        pad=0.025,
        fraction=0.050,
        format=LogFormatterMathtext(),
    )
    colorbar.set_label(r"$|\nabla\times\mathbf{B}|/|\mathbf{B}|$ [Mm$^{-1}$]", fontsize=14)
    colorbar.ax.tick_params(labelsize=11)
    figure.savefig(path, dpi=dpi, facecolor="white")
    plt.close(figure)


def save_gif(
    path: Path,
    frame_paths: Iterable[Path],
    *,
    fps: float,
    endpoint_pause_seconds: float,
) -> None:
    paths = list(frame_paths)
    if not paths:
        raise ValueError("Cannot save an empty GIF")
    adaptive = getattr(Image, "Palette", Image).ADAPTIVE
    dither_none = getattr(Image, "Dither", Image).NONE
    with Image.open(paths[0]) as source:
        first = source.convert("RGB").convert("P", palette=adaptive, colors=128)
    frames = [first]
    for frame_path in paths[1:]:
        with Image.open(frame_path) as source:
            frame = source.convert("RGB").quantize(palette=first, dither=dither_none)
        frames.append(frame)

    base_duration_ms = max(1, int(round(1000.0 / fps)))
    endpoint_duration_ms = max(
        base_duration_ms, int(round(1000.0 * endpoint_pause_seconds))
    )
    durations = [base_duration_ms] * len(frames)
    durations[0] = endpoint_duration_ms
    durations[-1] = endpoint_duration_ms
    first.save(
        path,
        save_all=True,
        append_images=frames[1:],
        duration=durations,
        loop=0,
        disposal=2,
        optimize=False,
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--magnetic-dir", type=Path, default=MAGNETIC_DIR)
    parser.add_argument("--output-dir", type=Path, default=OUTPUT_DIR)
    parser.add_argument("--start", type=int, default=START_SNAPSHOT)
    parser.add_argument("--stop", type=int, default=STOP_SNAPSHOT)
    parser.add_argument("--step", type=int, default=1)
    parser.add_argument("--fixed-y", type=float, default=FIXED_Y_CODE)
    parser.add_argument("--b-floor-code", type=float, default=1.0e-12)
    parser.add_argument(
        "--view-mm",
        type=float,
        nargs=4,
        metavar=("XMIN", "XMAX", "HMIN", "HMAX"),
        default=DEFAULT_VIEW_MM,
    )
    parser.add_argument("--lower-percentile", type=float, default=1.0)
    parser.add_argument("--upper-percentile", type=float, default=99.7)
    parser.add_argument("--vmin", type=float, default=None)
    parser.add_argument("--vmax", type=float, default=None)
    parser.add_argument("--fps", type=float, default=8.0)
    parser.add_argument("--endpoint-pause", type=float, default=1.0)
    parser.add_argument("--dpi", type=int, default=145)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.step <= 0:
        raise ValueError("--step must be positive")
    if args.stop < args.start:
        raise ValueError("--stop must be >= --start")
    if not np.isfinite(args.fixed_y) or args.fixed_y != 0.0:
        raise ValueError("This workflow is defined for exact --fixed-y 0")
    if not np.isfinite(args.b_floor_code) or args.b_floor_code < 0.0:
        raise ValueError("--b-floor-code must be finite and non-negative")
    if not (0.0 <= args.lower_percentile < args.upper_percentile <= 100.0):
        raise ValueError("Invalid display percentile range")
    if args.vmin is not None and (not np.isfinite(args.vmin) or args.vmin <= 0.0):
        raise ValueError("--vmin must be positive and finite")
    if args.vmax is not None and (not np.isfinite(args.vmax) or args.vmax <= 0.0):
        raise ValueError("--vmax must be positive and finite")
    if args.fps <= 0.0 or args.endpoint_pause < 0.0 or args.dpi <= 0:
        raise ValueError("fps and dpi must be positive; endpoint pause cannot be negative")

    snapshots = list(range(args.start, args.stop + 1, args.step))
    if snapshots[-1] != args.stop:
        raise ValueError("--step must land exactly on --stop")
    source_paths = [args.magnetic_dir / f"nlfff_{snapshot:04d}.binary" for snapshot in snapshots]
    missing = [path for path in source_paths if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"Missing {len(missing)} magnetic files; first: {missing[0]}")

    first_shape = read_header(source_paths[0])
    for path in source_paths[1:]:
        if read_header(path) != first_shape:
            raise ValueError(f"Magnetic grid shape changes at {path}")
    nx, ny, nz = first_shape
    x = cell_centres(DOMAIN["x"], nx)
    y = cell_centres(DOMAIN["y"], ny)
    z = cell_centres(DOMAIN["z"], nz)
    dx, dy, dz = float(x[1] - x[0]), float(y[1] - y[0]), float(z[1] - z[0])
    if not (np.allclose(np.diff(x), dx) and np.allclose(np.diff(y), dy) and np.allclose(np.diff(z), dz)):
        raise ValueError("Magnetic coordinate vectors are not uniform")
    if not np.allclose((dx, dy, dz), dx):
        raise ValueError(f"Magnetic cells are not cubic: dx,dy,dz={dx,dy,dz}")
    lower_y, upper_y, interpolation_weight = bracketing_y_indices(y, args.fixed_y)
    if not np.isclose(
        (1.0 - interpolation_weight) * y[lower_y] + interpolation_weight * y[upper_y],
        args.fixed_y,
        atol=1.0e-14,
    ):
        raise ValueError("Interpolated plane does not equal requested fixed y")

    x_mm = x * LENGTH_UNIT_MM
    height_code = z - BOTTOM_Z_CODE
    height_mm = height_code * LENGTH_UNIT_MM
    view_limits = tuple(float(value) for value in args.view_mm)
    keep_x, keep_z = display_selection(x_mm, height_mm, view_limits)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    data_dir = args.output_dir / "data"
    frame_dir = args.output_dir / "frames"
    data_dir.mkdir(parents=True, exist_ok=True)
    frame_dir.mkdir(parents=True, exist_ok=True)

    ratio_stack: list[np.ndarray] = []
    bvector_stack: list[np.ndarray] = []
    curlvector_stack: list[np.ndarray] = []
    bmag_stack: list[np.ndarray] = []
    curlmag_stack: list[np.ndarray] = []
    valid_stack: list[np.ndarray] = []
    statistics_rows: list[dict[str, Any]] = []
    source_records: list[dict[str, Any]] = []

    print(
        f"Calculating {len(snapshots)} exact-y=0 slices from y={y[lower_y]:+.6f}, "
        f"{y[upper_y]:+.6f}...",
        flush=True,
    )
    for index, (snapshot, source_path) in enumerate(zip(snapshots, source_paths), start=1):
        slabs = read_y_bracketing_slabs(
            source_path, first_shape, lower_y, upper_y
        )
        result = calculate_exact_y_slice(
            slabs,
            x=x,
            y_lower=float(y[lower_y]),
            y_upper=float(y[upper_y]),
            z=z,
            interpolation_weight=interpolation_weight,
            b_floor_code=args.b_floor_code,
        )
        ratio = result["j_over_b_per_Mm"]
        valid = result["valid"]
        if ratio.shape != (nz, nx) or valid.shape != (nz, nx):
            raise ValueError(f"Unexpected y=0 slice shape at snapshot {snapshot}")
        if np.any(np.isfinite(ratio) & (ratio < 0.0)):
            raise ValueError(f"Negative |J|/|B| value at snapshot {snapshot}")

        ratio_stack.append(ratio.astype(np.float32))
        bvector_stack.append(result["B_code"].astype(np.float32))
        curlvector_stack.append(result["curl_B_per_L0"].astype(np.float32))
        bmag_stack.append(result["B_magnitude_code"].astype(np.float32))
        curlmag_stack.append(result["curl_B_magnitude_per_L0"].astype(np.float32))
        valid_stack.append(valid.astype(np.uint8))
        time_minutes = snapshot * TIME_MINUTES_PER_SNAPSHOT

        full_stats = finite_stats(ratio)
        view_stats = finite_stats(ratio[np.ix_(keep_z, keep_x)])
        row: dict[str, Any] = {
            "snapshot": snapshot,
            "time_minutes": time_minutes,
            "valid_count": int(valid.sum()),
            "invalid_count": int((~valid).sum()),
        }
        row.update({f"full_{key}": value for key, value in full_stats.items()})
        row.update({f"view_{key}": value for key, value in view_stats.items()})
        statistics_rows.append(row)

        frame_npz = data_dir / f"y0_j_over_b_t{snapshot:04d}.npz"
        np.savez_compressed(
            frame_npz,
            schema_version=np.asarray(SCHEMA_VERSION, dtype=np.int16),
            snapshot=np.asarray(snapshot, dtype=np.int32),
            time_minutes=np.asarray(time_minutes, dtype=np.float64),
            fixed_y_amrvac_normalized=np.asarray(args.fixed_y, dtype=np.float64),
            y_bracketing_centres_amrvac=np.asarray(
                [y[lower_y], y[upper_y]], dtype=np.float64
            ),
            y_interpolation_weight=np.asarray(interpolation_weight, dtype=np.float64),
            x_amrvac_normalized=x.astype(np.float32),
            height_amrvac_normalized=height_code.astype(np.float32),
            z_absolute_amrvac_normalized=z.astype(np.float32),
            x_Mm=x_mm.astype(np.float32),
            height_Mm=height_mm.astype(np.float32),
            length_unit_Mm=np.asarray(LENGTH_UNIT_MM, dtype=np.float64),
            cell_spacing_amrvac_normalized=np.asarray(dx, dtype=np.float64),
            b_floor_code=np.asarray(args.b_floor_code, dtype=np.float64),
            B_code=result["B_code"].astype(np.float32),
            curl_B_per_L0=result["curl_B_per_L0"].astype(np.float32),
            B_magnitude_code=result["B_magnitude_code"].astype(np.float32),
            curl_B_magnitude_per_L0=result["curl_B_magnitude_per_L0"].astype(np.float32),
            j_over_b_per_L0=result["j_over_b_per_L0"].astype(np.float32),
            j_over_b_per_Mm=ratio.astype(np.float32),
            valid=valid.astype(np.uint8),
            source_magnetic_binary=np.asarray(str(source_path.resolve())),
        )
        source_stat = source_path.stat()
        source_records.append(
            {
                "snapshot": snapshot,
                "path": str(source_path.resolve()),
                "size_bytes": source_stat.st_size,
                "mtime_ns": source_stat.st_mtime_ns,
            }
        )
        print(
            f"  [{index:02d}/{len(snapshots)}] t={snapshot:04d}, "
            f"median={view_stats['median']:.4g}, p99.7={view_stats['p99_7']:.4g} Mm^-1",
            flush=True,
        )

    ratio_array = np.asarray(ratio_stack, dtype=np.float32)
    bvector_array = np.asarray(bvector_stack, dtype=np.float32)
    curlvector_array = np.asarray(curlvector_stack, dtype=np.float32)
    bmag_array = np.asarray(bmag_stack, dtype=np.float32)
    curlmag_array = np.asarray(curlmag_stack, dtype=np.float32)
    valid_array = np.asarray(valid_stack, dtype=np.uint8)
    times = np.asarray(snapshots, dtype=float) * TIME_MINUTES_PER_SNAPSHOT
    tag = f"t{snapshots[0]:04d}_t{snapshots[-1]:04d}"

    displayed = ratio_array[:, keep_z][:, :, keep_x]
    displayed_finite = displayed[np.isfinite(displayed) & (displayed > 0.0)]
    if not displayed_finite.size:
        raise ValueError("No positive finite J/B values in the displayed domain")
    automatic_vmin = float(np.percentile(displayed_finite, args.lower_percentile))
    automatic_vmax = float(np.percentile(displayed_finite, args.upper_percentile))
    vmin = float(args.vmin) if args.vmin is not None else automatic_vmin
    vmax = float(args.vmax) if args.vmax is not None else automatic_vmax
    if not (0.0 < vmin < vmax):
        raise ValueError(f"Invalid logarithmic color range: {vmin}, {vmax}")

    stack_path = args.output_dir / f"y0_j_over_b_{tag}_stack.npz"
    np.savez_compressed(
        stack_path,
        schema_version=np.asarray(SCHEMA_VERSION, dtype=np.int16),
        snapshots=np.asarray(snapshots, dtype=np.int32),
        time_minutes=times.astype(np.float64),
        fixed_y_amrvac_normalized=np.asarray(args.fixed_y, dtype=np.float64),
        y_bracketing_centres_amrvac=np.asarray(
            [y[lower_y], y[upper_y]], dtype=np.float64
        ),
        y_interpolation_weight=np.asarray(interpolation_weight, dtype=np.float64),
        x_amrvac_normalized=x.astype(np.float32),
        height_amrvac_normalized=height_code.astype(np.float32),
        z_absolute_amrvac_normalized=z.astype(np.float32),
        x_Mm=x_mm.astype(np.float32),
        height_Mm=height_mm.astype(np.float32),
        length_unit_Mm=np.asarray(LENGTH_UNIT_MM, dtype=np.float64),
        cell_spacing_amrvac_normalized=np.asarray(dx, dtype=np.float64),
        b_floor_code=np.asarray(args.b_floor_code, dtype=np.float64),
        B_code=bvector_array,
        curl_B_per_L0=curlvector_array,
        B_magnitude_code=bmag_array,
        curl_B_magnitude_per_L0=curlmag_array,
        j_over_b_per_L0=(ratio_array * LENGTH_UNIT_MM).astype(np.float32),
        j_over_b_per_Mm=ratio_array,
        valid=valid_array,
        source_magnetic_binaries=np.asarray([str(path.resolve()) for path in source_paths]),
    )

    csv_path = args.output_dir / f"y0_j_over_b_{tag}_statistics.csv"
    write_csv(csv_path, statistics_rows)

    print(f"Rendering {len(snapshots)} frames with fixed range {vmin:.4g}--{vmax:.4g} Mm^-1...", flush=True)
    frame_paths: list[Path] = []
    for index, (snapshot, time_minutes, ratio) in enumerate(
        zip(snapshots, times, ratio_array), start=1
    ):
        frame_path = frame_dir / f"y0_j_over_b_t{snapshot:04d}.png"
        draw_frame(
            frame_path,
            ratio_per_mm=ratio,
            snapshot=snapshot,
            time_minutes=float(time_minutes),
            x_mm=x_mm,
            height_mm=height_mm,
            limits=view_limits,
            vmin=vmin,
            vmax=vmax,
            dpi=args.dpi,
        )
        frame_paths.append(frame_path)
        if index == 1 or index % 10 == 0 or index == len(snapshots):
            print(f"  rendered [{index:02d}/{len(snapshots)}] {frame_path.name}", flush=True)

    gif_path = args.output_dir / f"y0_j_over_b_{tag}.gif"
    save_gif(
        gif_path,
        frame_paths,
        fps=args.fps,
        endpoint_pause_seconds=args.endpoint_pause,
    )

    # Round-trip checks ensure the reusable stack and animation are complete.
    with np.load(stack_path, allow_pickle=False) as archive:
        if archive["j_over_b_per_Mm"].shape != (len(snapshots), nz, nx):
            raise ValueError("Saved stack has the wrong J/B shape")
        if not np.array_equal(archive["snapshots"], np.asarray(snapshots)):
            raise ValueError("Saved stack snapshot list changed")
        stored = np.asarray(archive["j_over_b_per_Mm"], dtype=float)
        stored_b = np.asarray(archive["B_magnitude_code"], dtype=float)
        stored_curl = np.asarray(archive["curl_B_magnitude_per_L0"], dtype=float)
        stored_b_vector = np.asarray(archive["B_code"], dtype=float)
        stored_curl_vector = np.asarray(archive["curl_B_per_L0"], dtype=float)
        if not np.allclose(
            np.linalg.norm(stored_b_vector, axis=-1),
            stored_b,
            rtol=2.0e-6,
            atol=2.0e-7,
        ):
            raise ValueError("Saved B-vector stack disagrees with |B|")
        if not np.allclose(
            np.linalg.norm(stored_curl_vector, axis=-1),
            stored_curl,
            rtol=2.0e-6,
            atol=2.0e-7,
        ):
            raise ValueError("Saved curl(B)-vector stack disagrees with |curl(B)|")
        check = np.isfinite(stored) & (stored_b > args.b_floor_code)
        if not np.allclose(
            stored[check] * LENGTH_UNIT_MM * stored_b[check],
            stored_curl[check],
            rtol=2.0e-6,
            atol=2.0e-7,
        ):
            raise ValueError("Saved J/B stack fails the |curl B|/|B| identity")
    with Image.open(gif_path) as animation:
        gif_frames = int(getattr(animation, "n_frames", 1))
        gif_size = list(animation.size)
    if gif_frames != len(snapshots):
        raise ValueError(f"GIF contains {gif_frames} frames, expected {len(snapshots)}")

    metadata_path = args.output_dir / f"y0_j_over_b_{tag}_metadata.json"
    metadata = {
        "schema_version": SCHEMA_VERSION,
        "quantity": {
            "name": "curl-to-field ratio (normalized-current proxy)",
            "definition": "|J|/|B| = |curl(B)|/|B|",
            "amrvac_current_convention": "J = curl(B) in normalized code units",
            "stored_units": {
                "j_over_b_per_L0": "inverse normalized length L0^-1",
                "j_over_b_per_Mm": "Mm^-1",
                "B_magnitude_code": "normalized magnetic field",
                "curl_B_magnitude_per_L0": "normalized magnetic field per L0",
            },
            "length_unit_Mm": LENGTH_UNIT_MM,
            "b_floor_code": args.b_floor_code,
        },
        "slice": {
            "plane": "x-z at exact fixed y",
            "fixed_y_amrvac_normalized": args.fixed_y,
            "fixed_y_Mm": args.fixed_y * LENGTH_UNIT_MM,
            "bracketing_y_indices": [lower_y, upper_y],
            "bracketing_y_centres_amrvac": [float(y[lower_y]), float(y[upper_y])],
            "linear_interpolation_weight_toward_upper": interpolation_weight,
            "normal_derivative": "centered difference between bracketing y cell centres",
            "tangential_derivatives": "numpy.gradient, second-order edge treatment",
            "array_axis_order": "(time, z, x)",
            "vector_array_axis_order": "(time, z, x, component_xyz)",
        },
        "grid": {
            "input_shape_xyz": [nx, ny, nz],
            "saved_slice_shape_zx": [nz, nx],
            "cell_spacing_amrvac_normalized": dx,
            "cell_spacing_Mm": dx * LENGTH_UNIT_MM,
            "domain_amrvac_normalized": DOMAIN,
            "bottom_height_datum_amrvac_normalized": BOTTOM_Z_CODE,
        },
        "time": {
            "snapshots": snapshots,
            "time_minutes": times.tolist(),
            "minutes_per_snapshot": TIME_MINUTES_PER_SNAPSHOT,
        },
        "visualization": {
            "view_limits_Mm": list(view_limits),
            "scale": "logarithmic",
            "colormap": COLORMAP_NAME,
            "global_vmin_per_Mm": vmin,
            "global_vmax_per_Mm": vmax,
            "automatic_lower_percentile": args.lower_percentile,
            "automatic_upper_percentile": args.upper_percentile,
            "automatic_vmin_per_Mm": automatic_vmin,
            "automatic_vmax_per_Mm": automatic_vmax,
            "fps": args.fps,
            "endpoint_pause_seconds": args.endpoint_pause,
            "gif_frames": gif_frames,
            "gif_pixel_size": gif_size,
        },
        "quality": {
            "all_snapshots_present": len(source_paths) == len(snapshots),
            "all_values_nonnegative_or_nan": bool(
                np.all(~np.isfinite(ratio_array) | (ratio_array >= 0.0))
            ),
            "total_valid_count": int(valid_array.sum()),
            "total_invalid_count": int(valid_array.size - valid_array.sum()),
            "global_display_statistics_per_Mm": finite_stats(displayed),
            "stack_identity_check": "passed",
            "gif_frame_count_check": "passed",
        },
        "sources": source_records,
        "outputs": {
            "gif": str(gif_path.resolve()),
            "stack_npz": str(stack_path.resolve()),
            "statistics_csv": str(csv_path.resolve()),
            "per_snapshot_data_directory": str(data_dir.resolve()),
            "png_frame_directory": str(frame_dir.resolve()),
            "gif_sha256": sha256(gif_path),
            "stack_sha256": sha256(stack_path),
            "statistics_sha256": sha256(csv_path),
        },
    }
    with metadata_path.open("w", encoding="utf-8") as stream:
        json.dump(metadata, stream, indent=2, ensure_ascii=False)
        stream.write("\n")

    print(f"GIF: {gif_path}")
    print(f"Stack: {stack_path}")
    print(f"Per-frame data: {data_dir}")
    print(f"Statistics: {csv_path}")
    print(f"Metadata: {metadata_path}")


if __name__ == "__main__":
    main()

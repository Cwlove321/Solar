#!/usr/bin/env python3
"""Write projected AIA arrays as AMRVAC-aligned bottom-boundary VTK planes."""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import numpy as np


PARAMETERS = (
    "domain_nx1",
    "domain_nx2",
    "domain_nx3",
    "xprobmin1",
    "xprobmax1",
    "xprobmin2",
    "xprobmax2",
    "xprobmin3",
    "xprobmax3",
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--npy-dir", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--par", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument(
        "--transform",
        choices=("raw", "log10p1"),
        default="raw",
        help="Optional scalar transform before writing",
    )
    parser.add_argument(
        "--full-resolution",
        action="store_true",
        help="Use the unbinned projected NPY at twice the par x1/x2 resolution",
    )
    parser.add_argument(
        "--grid-scale",
        type=int,
        default=None,
        help="Explicit x1/x2 grid scale relative to domain_nx in the par file",
    )
    return parser.parse_args()


def parse_amrvac_par(path: Path) -> dict[str, int | float]:
    text = path.read_text(encoding="utf-8")
    values: dict[str, int | float] = {}
    for key in PARAMETERS:
        match = re.search(
            rf"\b{re.escape(key)}\s*=\s*([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[dDeE][+-]?\d+)?)",
            text,
        )
        if match is None:
            raise ValueError(f"missing {key} in {path}")
        raw = match.group(1).replace("D", "e").replace("d", "e")
        values[key] = int(float(raw)) if key.startswith("domain_nx") else float(raw)
    return values


def write_plane(
    data: np.ndarray,
    output: Path,
    params: dict[str, int | float],
    wavelength: int,
    nominal_time: str,
    date_obs: str,
    source: str,
    transform: str,
    grid_scale: int,
) -> dict[str, object]:
    nx = int(params["domain_nx1"]) * grid_scale
    ny = int(params["domain_nx2"]) * grid_scale
    if data.shape != (ny, nx):
        raise ValueError(f"{source}: array shape {data.shape} != par grid {(ny, nx)}")

    xmin = float(params["xprobmin1"])
    xmax = float(params["xprobmax1"])
    ymin = float(params["xprobmin2"])
    ymax = float(params["xprobmax2"])
    zmin = float(params["xprobmin3"])
    dz = (float(params["xprobmax3"]) - zmin) / int(params["domain_nx3"])

    dx = (xmax - xmin) / nx
    dy = (ymax - ymin) / ny
    scalar_name = (
        f"log10_AIA_{wavelength:03d}_plus_1"
        if transform == "log10p1"
        else f"AIA_{wavelength:03d}"
    )
    header = (
        "# vtk DataFile Version 3.0\n"
        f"AIA {wavelength} Angstrom bottom boundary at {nominal_time}\n"
        "BINARY\n"
        "DATASET STRUCTURED_POINTS\n"
        f"DIMENSIONS {nx + 1} {ny + 1} 1\n"
        f"ORIGIN {xmin:.17g} {ymin:.17g} {zmin:.17g}\n"
        f"SPACING {dx:.17g} {dy:.17g} {dz:.17g}\n"
        f"CELL_DATA {nx * ny}\n"
        f"SCALARS {scalar_name} float 1\n"
        "LOOKUP_TABLE default\n"
    )
    binary = np.ascontiguousarray(data, dtype=">f4").ravel(order="C").tobytes()
    with output.open("wb") as stream:
        stream.write(header.encode("ascii"))
        stream.write(binary)
        stream.write(b"\n")

    return {
        "vtk": output.name,
        "scalar": scalar_name,
        "shape_yx": [ny, nx],
        "dimensions_points": [nx + 1, ny + 1, 1],
        "number_of_cells": nx * ny,
        "bounds": [xmin, xmax, ymin, ymax, zmin, zmin],
        "origin": [xmin, ymin, zmin],
        "spacing": [dx, dy, dz],
        "wavelength_angstrom": wavelength,
        "nominal_time_utc": nominal_time,
        "date_obs_utc": date_obs,
        "source_fits": source,
        "transform": transform,
    }


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    params = parse_amrvac_par(args.par)
    metadata = json.loads(args.metadata.read_text(encoding="utf-8"))

    outputs: list[dict[str, object]] = []
    for record in metadata["records"]:
        grid_scale = args.grid_scale if args.grid_scale is not None else (2 if args.full_resolution else 1)
        if grid_scale < 1:
            raise ValueError("grid-scale must be at least 1")
        npy_key = "npy" if grid_scale > 1 else "bin2_npy"
        npy = args.npy_dir / record[npy_key]
        raw_data = np.load(npy)
        data = (
            np.log10(np.maximum(raw_data, 0.0) + 1.0).astype(np.float32)
            if args.transform == "log10p1"
            else raw_data
        )
        wavelength = int(record["wavelength_angstrom"])
        nominal_time = str(record.get("nominal_time", record["date_obs"]))
        time_tag = nominal_time.replace("-", "").replace(":", "")[0:15]
        resolution_suffix = (
            f"_{int(params['domain_nx1']) * grid_scale}x{int(params['domain_nx2']) * grid_scale}"
            if grid_scale > 1
            else ""
        )
        output = args.out / f"aia_{time_tag}_{wavelength:03d}_bottom_x3min{resolution_suffix}.vtk"
        item = write_plane(
                data,
                output,
                params,
                wavelength,
                nominal_time,
                str(record["date_obs"]),
                str(record["source"]),
                args.transform,
                grid_scale,
            )
        item["input_npy"] = npy.name
        outputs.append(item)

    manifest = {
        "amrvac_par": str(args.par),
        "parameters": params,
        "data_association": "CELL_DATA",
        "transform": args.transform,
        "transform_formula": (
            "log10(max(I, 0) + 1)" if args.transform == "log10p1" else "I"
        ),
        "grid_scale_relative_to_par": (
            args.grid_scale if args.grid_scale is not None else (2 if args.full_resolution else 1)
        ),
        "flattening": "NumPy (y,x) C-order; VTK x index varies fastest",
        "files": outputs,
    }
    (args.out / "vtk_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8"
    )


if __name__ == "__main__":
    main()

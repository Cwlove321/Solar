#!/usr/bin/env python3
"""Reproject AIA images onto a SHARP CEA pixel submap."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
import sunpy.map
from astropy.time import Time
from matplotlib.colors import AsinhNorm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sharp", type=Path, required=True)
    parser.add_argument("--aia-dir", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--x0", type=int, required=True)
    parser.add_argument("--x1", type=int, required=True)
    parser.add_argument("--y0", type=int, required=True)
    parser.add_argument("--y1", type=int, required=True)
    parser.add_argument(
        "--nominal-times",
        nargs="*",
        default=(),
        help="Requested UTC times used to group nearby wavelength records in plots",
    )
    return parser.parse_args()


def bin2_mean(data: np.ndarray) -> np.ndarray:
    ny, nx = data.shape
    if ny % 2 or nx % 2:
        raise ValueError(f"bin2 requires even dimensions, got {data.shape}")
    return np.nanmean(data.reshape(ny // 2, 2, nx // 2, 2), axis=(1, 3))


def output_stem(aia_map: sunpy.map.GenericMap, args: argparse.Namespace) -> str:
    timestamp = Time(aia_map.date).strftime("%Y%m%d_%H%M%S")
    wavelength = int(round(aia_map.wavelength.to_value(u.angstrom)))
    nx = args.x1 - args.x0
    ny = args.y1 - args.y0
    return f"aia_{timestamp}_{wavelength:03d}_sharp_cea_x{args.x0}_y{args.y0}_{nx}x{ny}"


def main() -> None:
    args = parse_args()
    args.out.mkdir(parents=True, exist_ok=True)

    sharp = sunpy.map.Map(args.sharp)
    target = sharp.submap(
        (args.x0, args.y0) * u.pix,
        top_right=(args.x1 - 1, args.y1 - 1) * u.pix,
    )
    expected_shape = (args.y1 - args.y0, args.x1 - args.x0)
    if target.data.shape != expected_shape:
        raise RuntimeError(f"target shape {target.data.shape} != {expected_shape}")

    records: list[dict[str, object]] = []
    nominal_times = [Time(value) for value in args.nominal_times]
    for aia_path in sorted(args.aia_dir.glob("*.fits")):
        aia = sunpy.map.Map(aia_path)
        projected = aia.reproject_to(
            target.wcs,
            shape_out=target.data.shape,
            algorithm="interpolation",
        )
        data = np.asarray(projected.data, dtype=np.float32)
        stem = output_stem(aia, args)
        npy_path = args.out / f"{stem}.npy"
        bin2_path = args.out / f"{stem}_bin2.npy"
        np.save(npy_path, data)
        np.save(bin2_path, bin2_mean(data).astype(np.float32))
        record = {
                "source": aia_path.name,
                "date_obs": aia.date.isot,
                "wavelength_angstrom": int(round(aia.wavelength.to_value(u.angstrom))),
                "npy": npy_path.name,
                "shape": list(data.shape),
                "bin2_npy": bin2_path.name,
                "bin2_shape": [data.shape[0] // 2, data.shape[1] // 2],
                "finite_fraction": float(np.isfinite(data).mean()),
            }
        if nominal_times:
            nearest = min(nominal_times, key=lambda value: abs((aia.date - value).sec))
            record["nominal_time"] = nearest.isot
        records.append(record)

    if not records:
        raise RuntimeError(f"no AIA FITS files found under {args.aia_dir}")

    records.sort(key=lambda item: (str(item["date_obs"]), int(item["wavelength_angstrom"])))
    metadata = {
        "sharp": str(args.sharp),
        "sharp_pixel_slice_python": {
            "x": [args.x0, args.x1],
            "y": [args.y0, args.y1],
        },
        "target_shape": list(target.data.shape),
        "bin2_shape": [target.data.shape[0] // 2, target.data.shape[1] // 2],
        "target_wcs": target.wcs.to_header().tostring(sep="\n"),
        "records": records,
    }
    (args.out / "projection_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    wavelengths = [131, 171, 211, 304]
    grouped: dict[str, dict[int, dict[str, object]]] = {}
    for record in records:
        label = str(record.get("nominal_time", record["date_obs"]))[11:16]
        grouped.setdefault(label, {})[int(record["wavelength_angstrom"])] = record
    times = sorted(grouped)
    fig, axes = plt.subplots(len(times), len(wavelengths), figsize=(14, 3.25 * len(times)), squeeze=False)
    cmaps = {131: "sdoaia131", 171: "sdoaia171", 211: "sdoaia211", 304: "sdoaia304"}
    for row, label in enumerate(times):
        for col, wavelength in enumerate(wavelengths):
            ax = axes[row, col]
            record = grouped[label].get(wavelength)
            if record is None:
                ax.set_axis_off()
                continue
            data = np.load(args.out / str(record["npy"]))
            finite = data[np.isfinite(data)]
            vmax = float(np.nanpercentile(finite, 99.7))
            ax.imshow(
                data,
                origin="lower",
                cmap=cmaps[wavelength],
                norm=AsinhNorm(linear_width=max(vmax * 0.02, 1e-6), vmin=0, vmax=vmax),
            )
            if row == 0:
                ax.set_title(f"AIA {wavelength} Å")
            if col == 0:
                ax.set_ylabel(f"{label} UT\ny pixel")
            else:
                ax.set_yticklabels([])
            if row == len(times) - 1:
                ax.set_xlabel("x pixel in SHARP submap")
            else:
                ax.set_xticklabels([])
    fig.suptitle(
        f"AIA reprojected to SHARP CEA: "
        f"x={args.x0}:{args.x1}, y={args.y0}:{args.y1}"
    )
    fig.tight_layout()
    fig.savefig(args.out / "aia_sharp_cea_mosaic.png", dpi=180, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

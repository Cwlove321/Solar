#!/usr/bin/env python3
"""Validate and render dual-view VTI output from AMRVAC thermal emission."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkIOXML import vtkXMLImageDataReader


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--top", type=Path, required=True)
    parser.add_argument("--side", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--wave", type=int, default=193)
    parser.add_argument("--frame", required=True)
    return parser.parse_args()


def aia_cmap(wave: int):
    try:
        import sunpy.visualization.colormaps  # noqa: F401

        return plt.get_cmap(f"sdoaia{wave}")
    except (ImportError, ValueError):
        return plt.get_cmap("inferno")


def read_vti(path: Path, wave: int) -> tuple[np.ndarray, list[float], str]:
    if not path.is_file() or path.stat().st_size == 0:
        raise FileNotFoundError(f"Missing or empty VTI: {path}")

    reader = vtkXMLImageDataReader()
    reader.SetFileName(str(path))
    reader.Update()
    image = reader.GetOutput()
    cell_data = image.GetCellData()
    names = [cell_data.GetArrayName(i) for i in range(cell_data.GetNumberOfArrays())]
    preferred = f"AIA{wave}"
    name = preferred if preferred in names else next(
        (item for item in names if item and item.upper().startswith("AIA")), None
    )
    if name is None:
        raise ValueError(f"No AIA CellData array in {path}; arrays={names}")

    values = vtk_to_numpy(cell_data.GetArray(name)).astype(np.float64, copy=False)
    nxp, nyp, nzp = image.GetDimensions()
    nx = max(nxp - 1, 1)
    ny = max(nyp - 1, 1)
    if values.size != nx * ny:
        raise ValueError(
            f"Unexpected VTI shape for {path}: dimensions={image.GetDimensions()}, "
            f"array_size={values.size}"
        )
    values = values.reshape(ny, nx)
    if not np.isfinite(values).all():
        raise ValueError(f"Non-finite AIA intensity in {path}")
    if not np.any(values > 0):
        raise ValueError(f"AIA intensity is entirely non-positive in {path}")
    return values, list(image.GetBounds()[:4]), name


def crop_positive(
    values: np.ndarray, bounds: list[float], margin_fraction: float = 0.08
) -> tuple[np.ndarray, list[float]]:
    rows, cols = np.where(values > 0)
    row_margin = max(4, int((rows.max() - rows.min() + 1) * margin_fraction))
    col_margin = max(4, int((cols.max() - cols.min() + 1) * margin_fraction))
    r0 = max(0, rows.min() - row_margin)
    r1 = min(values.shape[0], rows.max() + row_margin + 1)
    c0 = max(0, cols.min() - col_margin)
    c1 = min(values.shape[1], cols.max() + col_margin + 1)

    xmin, xmax, ymin, ymax = bounds
    x_edges = np.linspace(xmin, xmax, values.shape[1] + 1)
    y_edges = np.linspace(ymin, ymax, values.shape[0] + 1)
    return values[r0:r1, c0:c1], [x_edges[c0], x_edges[c1], y_edges[r0], y_edges[r1]]


def stats(values: np.ndarray) -> dict[str, float | int]:
    positive = values[values > 0]
    return {
        "shape_yx": list(values.shape),
        "positive_pixels": int(positive.size),
        "linear_min_positive": float(positive.min()),
        "linear_max": float(positive.max()),
        "log10_p01": float(np.percentile(np.log10(positive), 1)),
        "log10_p50": float(np.percentile(np.log10(positive), 50)),
        "log10_p997": float(np.percentile(np.log10(positive), 99.7)),
    }


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    inputs = {"top": args.top, "side": args.side}
    raw: dict[str, np.ndarray] = {}
    cropped: dict[str, tuple[np.ndarray, list[float]]] = {}
    arrays: dict[str, str] = {}
    all_positive = []
    metadata: dict[str, object] = {
        "wave_A": args.wave,
        "frame": args.frame,
        "radiation_model": "MPI-AMRVAC native optically thin thermal emission",
        "views": {},
    }

    for view, path in inputs.items():
        values, bounds, array_name = read_vti(path, args.wave)
        raw[view] = values
        cropped[view] = crop_positive(values, bounds)
        arrays[view] = array_name
        all_positive.append(values[values > 0])
        metadata["views"][view] = {
            "vti": str(path.resolve()),
            "array": array_name,
            "full_bounds_arcsec": bounds,
            **stats(values),
        }

    log_values = np.log10(np.concatenate(all_positive))
    # Exclude very faint PSF tails, which otherwise flatten all structure.
    vmin = float(np.percentile(log_values, 5))
    vmax = float(np.percentile(log_values, 99.5))
    if vmax <= vmin:
        vmax = vmin + 1.0
    metadata["shared_log10_limits"] = [vmin, vmax]

    cmap = aia_cmap(args.wave).copy()
    cmap.set_bad("black")
    labels = {"top": "Radial top view", "side": "South-to-north tangent view"}

    for view, (values, extent) in cropped.items():
        log_image = np.ma.masked_where(values <= 0, np.log10(np.maximum(values, 1e-300)))
        fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
        im = ax.imshow(
            log_image,
            origin="lower",
            extent=extent,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_title(f"AMRVAC native AIA {args.wave} - {labels[view]}")
        ax.set_xlabel("Detector X [arcsec]")
        ax.set_ylabel("Detector Y [arcsec]")
        fig.colorbar(im, ax=ax, label=r"log$_{10}$ intensity [DN s$^{-1}$ pixel$^{-1}$]")
        fig.savefig(
            args.output_dir / f"aia{args.wave}_{args.frame}_{view}.png",
            dpi=180,
            facecolor="white",
        )
        plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(13.2, 5.8), constrained_layout=True)
    image = None
    for ax, view in zip(axes, ("top", "side")):
        values, extent = cropped[view]
        log_image = np.ma.masked_where(values <= 0, np.log10(np.maximum(values, 1e-300)))
        image = ax.imshow(
            log_image,
            origin="lower",
            extent=extent,
            cmap=cmap,
            vmin=vmin,
            vmax=vmax,
            interpolation="nearest",
            aspect="equal",
        )
        ax.set_title(labels[view])
        ax.set_xlabel("Detector X [arcsec]")
        ax.set_ylabel("Detector Y [arcsec]")
    assert image is not None
    fig.colorbar(
        image,
        ax=axes,
        shrink=0.88,
        label=r"log$_{10}$ intensity [DN s$^{-1}$ pixel$^{-1}$]",
    )
    fig.suptitle(f"MPI-AMRVAC native optically thin AIA {args.wave} - {args.frame}")
    combined = args.output_dir / f"aia{args.wave}_{args.frame}_dualview.png"
    fig.savefig(combined, dpi=180, facecolor="white")
    plt.close(fig)

    metadata["combined_png"] = str(combined.resolve())
    metadata_path = args.output_dir / f"aia{args.wave}_{args.frame}_native.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    print(f"Validated and rendered native AIA output: {combined}")


if __name__ == "__main__":
    main()

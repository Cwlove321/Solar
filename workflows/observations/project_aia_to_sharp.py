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
from sunpy.coordinates import propagate_with_solar_surface


ROTATION_MODELS = ("howard", "snodgrass", "allen", "rigid")


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sharp", type=Path, required=True)
    sources = parser.add_mutually_exclusive_group(required=True)
    sources.add_argument("--aia-dir", type=Path)
    sources.add_argument("--aia", type=Path, nargs="+", help="Selected AIA FITS files.")
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--x0", type=int)
    parser.add_argument("--x1", type=int)
    parser.add_argument("--y0", type=int)
    parser.add_argument("--y1", type=int)
    parser.add_argument("--bin2", action="store_true", help="Additionally save 2x2 averages; requires even dimensions.")
    parser.add_argument("--raw-dn", action="store_true", help="Preserve source intensity units instead of exposure-normalized DN/s.")
    parser.add_argument(
        "--rotation-model",
        choices=ROTATION_MODELS,
        default="howard",
        help=(
            "Solar differential-rotation model used to propagate every AIA map "
            "to the SHARP target time (default: howard)"
        ),
    )
    parser.add_argument(
        "--no-differential-rotation",
        action="store_true",
        help="Reproject spatially without propagating AIA data to the SHARP target time",
    )
    parser.add_argument(
        "--nominal-times",
        nargs="*",
        default=(),
        help="Requested UTC times used to group nearby wavelength records in plots",
    )
    args = parser.parse_args(argv)
    crop = (args.x0, args.x1, args.y0, args.y1)
    if any(v is not None for v in crop) and not all(v is not None for v in crop):
        parser.error("Supply all of --x0/--x1/--y0/--y1, or omit them to use the full target.")
    return args


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


def reproject_to_sharp(
    aia_map: sunpy.map.GenericMap,
    target: sunpy.map.GenericMap,
    *,
    differential_rotation: bool = True,
    rotation_model: str = "howard",
) -> sunpy.map.GenericMap:
    """Reproject an AIA map onto a SHARP target, optionally aligning its time.

    The target WCS supplies both the SHARP CEA geometry and its observation
    time.  ``propagate_with_solar_surface`` makes the WCS transformation follow
    the rotating solar surface, so AIA maps observed at different times are all
    sampled at the same SHARP reference time.
    """

    reproject_kwargs = {
        "shape_out": target.data.shape,
        "algorithm": "interpolation",
    }
    if not differential_rotation:
        return aia_map.reproject_to(target.wcs, **reproject_kwargs)

    with propagate_with_solar_surface(rotation_model=rotation_model):
        return aia_map.reproject_to(target.wcs, **reproject_kwargs)


def project_aia_product(aia, target, *, differential_rotation=True, rotation_model="howard", normalize=True):
    """Project onto the target's exact WCS/time and preserve its finite ROI.

    Unlike the low-level compatibility API above, this returns a reloadable AIA
    map with target reference time, original source time, units and ROI mask.
    """
    source_time = aia.coordinate_frame.obstime
    reference = target.coordinate_frame.obstime
    exposure = aia.exposure_time.to_value(u.s)
    unit = str(aia.meta.get("bunit", "DN"))
    source = aia
    if normalize:
        rate_unit = aia.unit is not None and aia.unit.is_equivalent(u.DN/u.s)
        if not rate_unit:
            if aia.unit is not None and not aia.unit.is_equivalent(u.DN):
                raise ValueError(f"Expected DN or DN/s input; got {aia.unit}.")
            if not np.isfinite(exposure) or exposure <= 0:
                raise ValueError("AIA exposure must be positive to convert DN to DN/s.")
            meta = aia.meta.copy(); meta['bunit'] = 'DN / s'
            source = sunpy.map.Map(aia.data.astype(np.float32)/exposure, meta)
        unit = 'DN / s'
    projected = reproject_to_sharp(source, target, differential_rotation=differential_rotation, rotation_model=rotation_model)
    data = np.asarray(projected.data, dtype=np.float32).copy()
    data[~np.isfinite(target.data)] = np.nan
    header = projected.wcs.to_header()
    header['BUNIT'] = unit
    header['WAVELNTH'] = aia.wavelength.to_value(u.angstrom); header['WAVEUNIT'] = 'angstrom'
    header['TELESCOP'] = 'SDO/AIA'; header['INSTRUME'] = 'AIA'
    header['DATE-OBS'] = reference.utc.isot; header['T_OBS'] = reference.utc.isot
    header['SRC_TIME'] = source_time.utc.isot; header['REF_TIME'] = reference.utc.isot
    header['SRC_EXPT'] = float(exposure)
    header['EXPTIME'] = 1. if normalize else float(exposure)
    header['ROTMODEL'] = rotation_model if differential_rotation else 'none'
    header['CONTENT'] = 'AIA on target CEA grid; DATE-OBS is reference, SRC_TIME is source'
    header.add_history('Target CEA WCS and finite ROI preserved. Source time in SRC_TIME.')
    return sunpy.map.Map(data, header)


def main(argv=None) -> None:
    args = parse_args(argv)
    args.out.mkdir(parents=True, exist_ok=True)

    sharp = sunpy.map.Map(args.sharp)
    if args.x0 is None:
        args.x0, args.y0 = 0, 0
        args.y1, args.x1 = sharp.data.shape
    if not (0 <= args.x0 < args.x1 <= sharp.data.shape[1] and 0 <= args.y0 < args.y1 <= sharp.data.shape[0]):
        raise ValueError("Requested pixel crop falls outside the CEA target.")
    target = sharp.submap(
        (args.x0, args.y0) * u.pix,
        top_right=(args.x1 - 1, args.y1 - 1) * u.pix,
    )
    expected_shape = (args.y1 - args.y0, args.x1 - args.x0)
    if target.data.shape != expected_shape:
        raise RuntimeError(f"target shape {target.data.shape} != {expected_shape}")
    if args.bin2 and any(n % 2 for n in target.data.shape):
        raise ValueError("--bin2 requires even target dimensions; omit it to preserve the full grid.")

    differential_rotation = not args.no_differential_rotation
    reference_time = Time(target.coordinate_frame.obstime)

    records: list[dict[str, object]] = []
    nominal_times = [Time(value) for value in args.nominal_times]
    aia_paths = args.aia if args.aia is not None else sorted(p for p in args.aia_dir.glob("*.fits") if not p.name.startswith('._'))
    for aia_path in aia_paths:
        aia = sunpy.map.Map(aia_path)
        projected = project_aia_product(
            aia,
            target,
            differential_rotation=differential_rotation,
            rotation_model=args.rotation_model,
            normalize=not args.raw_dn,
        )
        data = np.asarray(projected.data, dtype=np.float32)
        stem = output_stem(aia, args)
        npy_path = args.out / f"{stem}.npy"
        bin2_path = args.out / f"{stem}_bin2.npy"
        np.save(npy_path, data)
        fits_path = args.out / f"{stem}.fits"
        projected.save(fits_path, overwrite=True)
        if args.bin2:
            np.save(bin2_path, bin2_mean(data).astype(np.float32))
        record = {
            "source": str(aia_path.resolve()),
            "date_obs": aia.date.isot,
            "source_reference_time_utc": aia.coordinate_frame.obstime.utc.isot,
            "fits": fits_path.name,
            "unit": str(projected.unit),
            "source_exposure_seconds": float(aia.exposure_time.to_value(u.s)),
            "time_offset_to_sharp_seconds": float(
                (reference_time - Time(aia.coordinate_frame.obstime)).to_value(u.s)
            ),
            "wavelength_angstrom": int(round(aia.wavelength.to_value(u.angstrom))),
            "npy": npy_path.name,
            "shape": list(data.shape),
            "bin2_npy": bin2_path.name if args.bin2 else None,
            "bin2_shape": [data.shape[0] // 2, data.shape[1] // 2] if args.bin2 else None,
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
        "bin2_shape": [target.data.shape[0] // 2, target.data.shape[1] // 2] if args.bin2 else None,
        "exposure_normalized": not args.raw_dn,
        "target_roi_mask_preserved": True,
        "target_wcs": target.wcs.to_header().tostring(sep="\n"),
        "temporal_alignment": {
            "reference": "SHARP target WCS observation time",
            "reference_time_utc": reference_time.utc.isot,
            "differential_rotation_enabled": differential_rotation,
            "rotation_model": args.rotation_model if differential_rotation else None,
            "implementation": "sunpy.coordinates.propagate_with_solar_surface",
        },
        "records": records,
    }
    (args.out / "projection_metadata.json").write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    wavelengths = sorted({int(record['wavelength_angstrom']) for record in records})
    grouped: dict[str, dict[int, dict[str, object]]] = {}
    for record in records:
        label = str(record.get("nominal_time", record["date_obs"]))[:16]
        grouped.setdefault(label, {})[int(record["wavelength_angstrom"])] = record
    times = sorted(grouped)
    fig, axes = plt.subplots(len(times), len(wavelengths), figsize=(14, 3.25 * len(times)), squeeze=False)
    cmaps = {wave: f"sdoaia{wave}" for wave in wavelengths}
    for row, label in enumerate(times):
        for col, wavelength in enumerate(wavelengths):
            ax = axes[row, col]
            record = grouped[label].get(wavelength)
            if record is None:
                ax.set_axis_off()
                continue
            data = np.load(args.out / str(record["npy"]))
            finite = data[np.isfinite(data)]
            if not finite.size:
                ax.text(.5,.5,'No valid overlap',ha='center',transform=ax.transAxes)
                continue
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

#!/usr/bin/env python
"""Batch synthetic radiation images from Cartesian AMRVAC ``.dat`` snapshots."""

from __future__ import annotations

import argparse
import gc
from collections import namedtuple
from glob import glob
from pathlib import Path
import sys
import traceback

import numpy as np

REPO_SRC = Path(__file__).resolve().parents[3] / "src"
if REPO_SRC.exists() and str(REPO_SRC) not in sys.path:
    sys.path.insert(0, str(REPO_SRC))

from solarutils.sim.radiation import (
    RadiationPlane,
    SlabBounds,
    density_pressure_temperature_from_amrvac,
    log_intensity_image,
    plot_radiation_image,
    synthesize_radiation_image,
)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.simesh_path:
        sys.path.insert(0, str(Path(args.simesh_path).expanduser().resolve()))

    try:
        from simesh import amr_loader
    except ImportError:
        from simesh.legacy.frontends.amrvac.io import amr_loader

    dat_files = sorted(glob(args.dat_glob))
    if not dat_files:
        raise FileNotFoundError(f"No AMRVAC .dat files matched: {args.dat_glob}")

    output_dir = Path(args.output_dir).expanduser()
    output_dir.mkdir(parents=True, exist_ok=True)
    plane = RadiationPlane.from_notebook_defaults(normal=args.normal, resolution=args.resolution)
    slab = SlabBounds(lower=args.slab_lower, upper=args.slab_upper)
    wavelengths = [parse_wavelength(item) for item in args.wavelength]
    indices = select_indices(args, len(dat_files))

    for file_index in indices:
        dat_path = Path(dat_files[file_index])
        try:
            print(f"[{file_index:04d}] loading {dat_path}")
            number_density, pressure, temperature = load_uniform_physical_cube(
                dat_path,
                amr_loader=amr_loader,
                slab=slab,
                refinement=args.uniform_refinement,
            )
            if args.save_physical:
                np.savez_compressed(
                    output_dir / f"data_{file_index:04d}.npz",
                    numberdensity=number_density,
                    pressure=pressure,
                    temperature=temperature,
                    source=str(dat_path),
                )

            for wavelength in wavelengths:
                image = synthesize_radiation_image(
                    number_density,
                    pressure,
                    temperature,
                    plane,
                    wavelength=wavelength,
                    absorption=args.absorption,
                    dda_root=args.dda_root,
                )
                image_log = log_intensity_image(image)
                wave_dir = output_dir / f"{'abs' if args.absorption else 'thin'}_{wavelength}"
                wave_dir.mkdir(parents=True, exist_ok=True)
                np.save(wave_dir / f"image_{file_index:04d}.npy", image_log)
                if args.plot:
                    plot_radiation_image(
                        image_log,
                        wavelength=wavelength,
                        frame=file_index,
                        output=wave_dir / f"aia{wavelength}_{file_index:04d}.png",
                        extent=args.extent,
                        vmin=args.vmin,
                        vmax=args.vmax,
                        crop_pixels=args.crop_pixels,
                        image_is_log=True,
                    )
            del number_density, pressure, temperature
            gc.collect()
        except Exception as exc:
            print(f"[{file_index:04d}] failed: {exc}")
            if args.stop_on_error:
                raise
            traceback.print_exc()
            continue

    return 0


def load_uniform_physical_cube(dat_path, *, amr_loader, slab: SlabBounds, refinement: int):
    """Load one AMRVAC file and export notebook-compatible physical cubes."""

    ds = amr_loader(str(dat_path))
    try:
        ds.update()
        nx, ny, nz = refinement * np.asarray(ds.header["domain_nx"], dtype=int)
        names = ds.header["w_names"]
        variables = namedtuple("W", names)
        slab_all = ds.mesh.export_slab_uniform(slab.lower, slab.upper, int(nx), int(ny), int(nz))
        w = variables(*np.transpose(slab_all))
        if not hasattr(w, "rho") or not hasattr(w, "e"):
            raise AttributeError("AMRVAC variables must contain rho and e")
        return density_pressure_temperature_from_amrvac(w.rho, w.e)
    finally:
        del ds
        gc.collect()


def parse_wavelength(value: str) -> int | str:
    if value in {"radio", "white_light"}:
        return value
    return int(value)


def select_indices(args, n_files: int) -> list[int]:
    if args.indices:
        indices = []
        for item in args.indices:
            if ":" in item:
                parts = [int(part) if part else None for part in item.split(":")]
                indices.extend(range(n_files)[slice(*parts)])
            else:
                indices.append(int(item))
    else:
        stop = args.stop if args.stop is not None else n_files
        indices = list(range(args.start, stop, args.step))
    bad = [index for index in indices if index < 0 or index >= n_files]
    if bad:
        raise IndexError(f"Frame index out of range for {n_files} files: {bad[0]}")
    return indices


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dat-glob", required=True, help="Glob for AMRVAC .dat files.")
    parser.add_argument("--output-dir", required=True, help="Directory for radsyn products.")
    parser.add_argument("--dda-root", help="Directory containing RST_cart and AIA_response.")
    parser.add_argument("--simesh-path", help="Optional path to prepend for importing simesh.")
    parser.add_argument("--wavelength", nargs="+", default=["171"], help="AIA wavelength(s), radio, or white_light.")
    parser.add_argument("--start", type=int, default=0, help="First file index when --indices is not used.")
    parser.add_argument("--stop", type=int, help="Stop file index when --indices is not used.")
    parser.add_argument("--step", type=int, default=1, help="Step between file indices.")
    parser.add_argument("--indices", nargs="*", help="Explicit file indices or Python slices, e.g. 1 4 10:20:3.")
    parser.add_argument("--uniform-refinement", type=int, default=8, help="Multiplier for domain_nx.")
    parser.add_argument("--resolution", type=int, default=512, help="Detector plane pixels per side.")
    parser.add_argument("--normal", type=float, nargs=3, default=(0.0, 1.0, 0.0), help="LOS plane normal.")
    parser.add_argument("--slab-lower", type=float, nargs=3, default=(-1.5, -3.0, 0.1), help="Slab lower corner.")
    parser.add_argument("--slab-upper", type=float, nargs=3, default=(4.5, 3.0, 12.1), help="Slab upper corner.")
    parser.add_argument("--extent", type=float, nargs=4, default=(-1.5, 4.5, 0.0, 12.0), help="Plot extent in Mm.")
    parser.add_argument("--crop-pixels", type=int, default=128, help="Side crop for quicklook plots.")
    parser.add_argument("--vmin", type=float, default=0.0, help="Quicklook lower color limit.")
    parser.add_argument("--vmax", type=float, help="Quicklook upper color limit.")
    parser.add_argument("--absorption", action="store_true", help="Include EUV absorption.")
    parser.add_argument("--save-physical", action="store_true", help="Save number density, pressure, and temperature cubes.")
    parser.add_argument("--plot", action="store_true", help="Save PNG quicklooks in addition to NPY arrays.")
    parser.add_argument("--stop-on-error", action="store_true", help="Abort on the first failed frame.")
    return parser.parse_args(argv)


if __name__ == "__main__":
    raise SystemExit(main())

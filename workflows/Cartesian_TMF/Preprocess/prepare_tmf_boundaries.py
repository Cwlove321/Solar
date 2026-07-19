#!/usr/bin/env python3
"""Generate data-driven TMF, potential, and NLFFF boundaries from SHARP FITS.

This command-line entry point runs the project-local ``driven_boundary``
DAVE4VM workflow and adds input validation, portable arguments, macOS metadata
filtering, and an output manifest.
"""

from __future__ import annotations

import argparse
import json
from glob import glob
from pathlib import Path

import numpy as np


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", required=True, type=Path, help="Directory containing chronological *.Br/.Bt/.Bp.fits files.")
    parser.add_argument("--output-dir", required=True, type=Path, help="Event preprocessing output directory.")
    parser.add_argument("--x1", required=True, type=int, help="Inclusive Python x start pixel of the SHARP crop.")
    parser.add_argument("--x2", required=True, type=int, help="Exclusive Python x end pixel of the SHARP crop.")
    parser.add_argument("--y1", required=True, type=int, help="Inclusive Python y start pixel of the SHARP crop.")
    parser.add_argument("--y2", required=True, type=int, help="Exclusive Python y end pixel of the SHARP crop.")
    parser.add_argument("--level", type=int, default=2, help="TMF multigrid level; level=2 is 2x2 binning.")
    parser.add_argument(
        "--prep", type=int, default=0,
        help="TMF B-boundary preprocessing; default 0 preserves the observed time evolution.",
    )
    parser.add_argument("--driven-fast", type=float, default=12.0, help="Temporal acceleration factor used by TMF.")
    parser.add_argument("--window-size", type=int, default=23, help="DAVE4VM window half-width/control size.")
    parser.add_argument("--threshold", type=float, default=100.0, help="Field-strength threshold in G for velocity cleanup.")
    parser.add_argument("--nlfff-level", type=int, default=2, help="Multigrid level for potential/NLFFF boundaries.")
    parser.add_argument(
        "--nlfff-prep", type=int, default=0,
        help="Initial NLFFF preprocessing; default 0 preserves the first observed field for TMF.",
    )
    parser.add_argument("--diagnostic-plots", action="store_true", help="Keep per-frame DAVE4VM velocity PNGs; disabled by default.")
    return parser.parse_args(argv)


def clean_glob(pattern: str) -> list[str]:
    """Return FITS paths while ignoring macOS AppleDouble ``._*`` files."""

    return sorted(path for path in glob(pattern) if not Path(path).name.startswith("._"))


def validate_inputs(args: argparse.Namespace) -> dict[str, list[str]]:
    if not args.input_dir.is_dir():
        raise FileNotFoundError(f"input directory does not exist: {args.input_dir}")
    if not (args.x2 > args.x1 and args.y2 > args.y1):
        raise ValueError("crop bounds must satisfy x2>x1 and y2>y1")
    if args.level < 1 or args.nlfff_level < 1:
        raise ValueError("--level and --nlfff-level must be positive")
    components = {
        "Br": clean_glob(str(args.input_dir / "*.Br.fits")),
        "Bt": clean_glob(str(args.input_dir / "*.Bt.fits")),
        "Bp": clean_glob(str(args.input_dir / "*.Bp.fits")),
    }
    counts = {name: len(paths) for name, paths in components.items()}
    if min(counts.values()) < 2:
        raise ValueError(f"at least two frames of each component are required; found {counts}")
    if len(set(counts.values())) != 1:
        raise ValueError(f"Br/Bt/Bp frame counts must agree; found {counts}")
    return components


def read_boundary_header(path: Path, components: int) -> dict[str, object]:
    with path.open("rb") as stream:
        nx, ny = np.fromfile(stream, dtype=np.int32, count=2)
        xc, yc, dx, dy = np.fromfile(stream, dtype=np.float64, count=4)
    return {
        "path": str(path),
        "nx": int(nx),
        "ny": int(ny),
        "components": components,
        "xc_cm": float(xc),
        "yc_cm": float(yc),
        "dx_cm": float(dx),
        "dy_cm": float(dy),
    }


def read_driven_header(path: Path) -> dict[str, object]:
    with path.open("rb") as stream:
        nx, ny, nt = np.fromfile(stream, dtype=np.int32, count=3)
        dx, dy, dz = np.fromfile(stream, dtype=np.float64, count=3)
    return {
        "path": str(path),
        "nx": int(nx),
        "ny": int(ny),
        "nt": int(nt),
        "dx_km": float(dx),
        "dy_km": float(dy),
        "dz_km": float(dz),
    }


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    components = validate_inputs(args)
    output = args.output_dir
    boundary_dir = output / "tmf" / "bd"
    figure_dir = output / "fig" / "v_fig"
    boundary_dir.mkdir(parents=True, exist_ok=True)
    figure_dir.mkdir(parents=True, exist_ok=True)

    try:
        from .prepare_nlfff_boundary import write_nlfff_from_cartesian  # noqa: PLC0415
        from . import driven_boundary  # noqa: PLC0415
    except ImportError:  # Direct execution: python Preprocess/prepare_tmf_boundaries.py
        from prepare_nlfff_boundary import write_nlfff_from_cartesian  # noqa: PLC0415
        import driven_boundary  # noqa: PLC0415

    # Keep the time-series/DAVE calculation in driven_boundary, but use the project
    # single-frame processor for the initial potential and NLFFF files.  This
    # prevents a second, subtly different NLFFF writer from drifting away from
    # the standalone workflow.
    def write_initial_nlfff(outpath, bx, by, bz, nlfff_level, nlfff_prep=0, xc=0, yc=0):
        arcsec_to_cm = 7.32825e7
        return write_nlfff_from_cartesian(
            Path(outpath) / "tmf" / "bd",
            bx,
            by,
            bz,
            level=nlfff_level,
            preprocess=bool(nlfff_prep),
            # Native SHARP CEA spacing; the shared NLFFF writer applies the
            # requested multigrid factor exactly once in the file header.
            dx_cm=0.50 * arcsec_to_cm,
            xc_cm=xc,
            yc_cm=yc,
        )

    # Replace the module-level glob with the metadata-safe version.
    driven_boundary.glob = clean_glob
    driven_boundary.write_nlfff_boundary = write_initial_nlfff
    if not args.diagnostic_plots:
        driven_boundary.vector_plot = lambda *unused_args, **unused_kwargs: None

    driven_boundary.write_driven_boundary(
        str(args.input_dir) + "/",
        str(output) + "/",
        x1=args.x1,
        x2=args.x2,
        y1=args.y1,
        y2=args.y2,
        level=args.level,
        prep=args.prep,
        driven_fast=args.driven_fast,
        window_size=args.window_size,
        thr=args.threshold,
        nlfff_level=args.nlfff_level,
        nlfff_prep=args.nlfff_prep,
    )

    utc_path = boundary_dir / "UTC.npy"
    utc = np.load(utc_path, allow_pickle=True)
    metadata = {
        "inputs": {name: {"count": len(paths), "first": paths[0], "last": paths[-1]} for name, paths in components.items()},
        "crop_original_pixels": {"x1": args.x1, "x2": args.x2, "y1": args.y1, "y2": args.y2},
        "tmf": {
            "level": args.level,
            "prep": args.prep,
            "driven_fast": args.driven_fast,
            "window_size": args.window_size,
            "threshold_gauss": args.threshold,
            "boundary": read_driven_header(boundary_dir / "B_boundary.dat"),
            "utc_first": str(utc[0]),
            "utc_last": str(utc[-1]),
        },
        "nlfff": {
            "level": args.nlfff_level,
            "prep": args.nlfff_prep,
            "potential_boundary": read_boundary_header(boundary_dir / "potential_boundary.dat", components=1),
            "nlfff_boundary": read_boundary_header(boundary_dir / "nlfff_boundary.dat", components=3),
            "note": "nlfff_boundary includes two zero-padded ghost pixels on each side.",
        },
    }
    (output / "tmf_preprocess_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    print(f"Wrote TMF products to {boundary_dir}")
    print(f"Wrote metadata to {output / 'tmf_preprocess_metadata.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

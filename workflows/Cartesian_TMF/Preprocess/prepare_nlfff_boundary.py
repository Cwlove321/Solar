#!/usr/bin/env python
"""Prepare Guo-style Cartesian potential/NLFFF bottom-boundary files from CEA FITS.

This bundled preprocessing entry point prepares the initial potential/NLFFF
boundary used by the portable ``Cartesian_TMF`` project.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

try:
    from .nlfff_core.cartesian import (
        SOLAR_RADIUS_CM,
        boundary_metadata,
        cartesian_preprocess_diagnostics,
        crop_components_with_padding,
        format_guo_amrvac_parameters,
        local_cartesian_from_heliographic,
        multigrid_reduce_components,
        multigrid_reducer,
        preprocess_cartesian_field,
        write_guo_allboundaries,
        write_guo_grid_ini,
        write_guo_nd_ini,
        write_nlfff_boundary,
        write_potential_boundary,
    )
except ImportError:  # Direct execution: python Preprocess/prepare_nlfff_boundary.py
    from nlfff_core.cartesian import (
        SOLAR_RADIUS_CM,
        boundary_metadata,
        cartesian_preprocess_diagnostics,
        crop_components_with_padding,
        format_guo_amrvac_parameters,
        local_cartesian_from_heliographic,
        multigrid_reduce_components,
        multigrid_reducer,
        preprocess_cartesian_field,
        write_guo_allboundaries,
        write_guo_grid_ini,
        write_guo_nd_ini,
        write_nlfff_boundary,
        write_potential_boundary,
    )


DEFAULT_NLFFF_PREPROCESS_MU3 = 0.1
DEFAULT_NLFFF_PREPROCESS_MU4 = 0.1
DEFAULT_NLFFF_PREPROCESS_MAX_ITER = 5000
DEFAULT_NLFFF_PREPROCESS_TOL = 1.0e-4


def validate_physical_grid(nx: int, ny: int) -> dict[str, object]:
    """Enforce this workflow's 2/3/5-factor rule on the unpadded binned grid."""

    factors = {}
    invalid = []
    for name, value in (("nx", nx), ("ny", ny)):
        if isinstance(value, bool) or not isinstance(value, (int, np.integer)) or value <= 0:
            raise ValueError(f"physical {name} must be a positive integer; got {value!r}")
        remaining = int(value)
        powers = {}
        for prime in (2, 3, 5):
            exponent = 0
            while remaining % prime == 0:
                remaining //= prime
                exponent += 1
            if exponent:
                powers[str(prime)] = exponent
        factors[name] = powers
        if remaining != 1:
            invalid.append(f"{name}={value} (remaining factor {remaining})")
    if invalid:
        raise ValueError(
            "physical nx/ny after binning, excluding file ghosts, must contain only "
            "prime factors 2, 3, 5: " + "; ".join(invalid)
            + ". Adjust the crop before writing boundary files; e.g. 360 x 216."
        )
    return {
        "passed": True,
        "physical_shape_xy": [int(nx), int(ny)],
        "allowed_prime_factors": [2, 3, 5],
        "factorization": factors,
        "excludes_file_ghosts": True,
    }


def write_nlfff_from_cartesian(
    output_dir: Path,
    bx: np.ndarray,
    by: np.ndarray,
    bz: np.ndarray,
    *,
    level: int,
    preprocess: bool = False,
    preprocess_mu3: float = DEFAULT_NLFFF_PREPROCESS_MU3,
    preprocess_mu4: float = DEFAULT_NLFFF_PREPROCESS_MU4,
    preprocess_max_iter: int = DEFAULT_NLFFF_PREPROCESS_MAX_ITER,
    preprocess_tol: float = DEFAULT_NLFFF_PREPROCESS_TOL,
    nghost: int = 2,
    dx_cm: float,
    dy_cm: float | None = None,
    xc_cm: float = 0.0,
    yc_cm: float = 0.0,
    z_min_10mm: float = 0.1,
) -> dict[str, object]:
    """Write potential/NLFFF products from Cartesian physical pixels.

    This is the shared non-CLI part of the existing NLFFF workflow.  TMF uses
    it for its initial frame, while this file's command-line interface reads
    and crops a single CEA FITS map before applying the same conventions.
    """

    reducer = multigrid_reducer(level)
    if bx.shape != by.shape or bx.shape != bz.shape:
        raise ValueError("Bx, By, and Bz must have identical shapes")
    if bx.ndim != 2:
        raise ValueError("Bx, By, and Bz must be two-dimensional")
    if nghost < 0:
        raise ValueError("nghost must be non-negative")
    if bx.shape[0] % reducer or bx.shape[1] % reducer:
        raise ValueError(f"array shape {bx.shape} is not divisible by reducer={reducer}")
    grid_validation = validate_physical_grid(bx.shape[1] // reducer, bx.shape[0] // reducer)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    bx, by, bz = multigrid_reduce_components(bx, by, bz, level=level)
    if preprocess:
        bx, by, bz, _ = preprocess_cartesian_field(
            bx,
            by,
            bz,
            mu3=preprocess_mu3,
            mu4=preprocess_mu4,
            max_iter=preprocess_max_iter,
            tol=preprocess_tol,
        )
    # The Guo writers expect a ghost-padded array: potential output strips the
    # file ghosts, while NLFFF output retains them.
    bx = np.pad(bx, ((nghost, nghost), (nghost, nghost)), mode="constant")
    by = np.pad(by, ((nghost, nghost), (nghost, nghost)), mode="constant")
    bz = np.pad(bz, ((nghost, nghost), (nghost, nghost)), mode="constant")
    dy_cm = dx_cm if dy_cm is None else dy_cm
    meta = boundary_metadata(
        nx=bx.shape[1], ny=bx.shape[0], xc_cm=xc_cm, yc_cm=yc_cm,
        dx_cm=dx_cm * reducer, dy_cm=dy_cm * reducer, nghost=nghost,
        z_min_10mm=z_min_10mm,
    )
    write_potential_boundary(output_dir / "potential_boundary.dat", bz, meta)
    write_nlfff_boundary(output_dir / "nlfff_boundary.dat", bx, by, bz, meta)
    return {"meta": meta, "grid_validation": grid_validation}


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    reducer = multigrid_reducer(args.level)
    grid_validation = validate_physical_window(args, reducer)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    print(f"Validated physical grid {args.nx // reducer} x {args.ny // reducer}: prime factors 2/3/5.")

    br_path, bt_path, bp_path = resolve_component_paths(args)
    print("Loading CEA vector FITS...")
    br, br_header = read_fits_data(br_path)
    bt, _ = read_fits_data(bt_path)
    bp, _ = read_fits_data(bp_path)
    bx, by, bz = local_cartesian_from_heliographic(br, bt, bp)

    ghost_pixels = args.nghost * reducer

    print("Cropping potential window plus ghost pixels...")
    bx_crop, by_crop, bz_crop, padding = crop_components_with_padding(
        bx,
        by,
        bz,
        x0=args.x0,
        y0=args.y0,
        nx=args.nx,
        ny=args.ny,
        pad_x=ghost_pixels,
        pad_y=ghost_pixels,
        fill_value=args.blank_value,
    )

    print(f"Applying multigrid level={args.level} reducer={reducer}...")
    bx_mg, by_mg, bz_mg = multigrid_reduce_components(bx_crop, by_crop, bz_crop, level=args.level)

    preprocess_metrics = None
    diagnostics_before = cartesian_preprocess_diagnostics(bx_mg, by_mg, bz_mg)
    diagnostics_after = diagnostics_before
    if args.preprocess:
        print("Preprocessing vector boundary with Wiegelmann-style optimization...")
        bx_mg, by_mg, bz_mg, preprocess_metrics = preprocess_cartesian_field(
            bx_mg,
            by_mg,
            bz_mg,
            mu3=args.preprocess_mu3,
            mu4=args.preprocess_mu4,
            max_iter=args.preprocess_max_iter,
            tol=args.preprocess_tol,
            report_interval=args.preprocess_report_interval,
            progress=print_preprocess_progress,
        )
        diagnostics_after = cartesian_preprocess_diagnostics(bx_mg, by_mg, bz_mg)
    else:
        print("Skipping preprocessing.")

    dx_cm, dy_cm = infer_spacing_cm(br_header, reducer, args)
    xc_cm, yc_cm = infer_window_center_cm(br_header, args)
    meta = boundary_metadata(
        nx=bx_mg.shape[1],
        ny=bx_mg.shape[0],
        xc_cm=xc_cm,
        yc_cm=yc_cm,
        dx_cm=dx_cm,
        dy_cm=dy_cm,
        nghost=args.nghost,
        z_min_10mm=args.z_min_10mm,
    )

    print("Writing Guo/AMRVAC boundary products...")
    write_guo_allboundaries(output_dir / "allboundaries.dat", bx_mg, by_mg, bz_mg)
    write_guo_grid_ini(output_dir / "grid.ini", nx=meta.nx, ny=meta.ny)
    write_guo_nd_ini(output_dir / "nd.ini", nxmax=bx_crop.shape[1], nymax=bx_crop.shape[0])
    write_potential_boundary(output_dir / "potential_boundary.dat", bz_mg, meta)
    write_nlfff_boundary(output_dir / "nlfff_boundary.dat", bx_mg, by_mg, bz_mg, meta)

    parameters_txt = format_guo_amrvac_parameters(
        meta,
        domain_nx3=args.domain_nx3,
        block_nx1=args.block_nx1,
        block_nx2=args.block_nx2,
        block_nx3=args.block_nx3,
    )
    (output_dir / "amrvac_parameters.txt").write_text(parameters_txt, encoding="utf-8")

    metadata = {
        "inputs": {
            "Br": str(br_path),
            "Bt": str(bt_path),
            "Bp": str(bp_path),
        },
        "selection": {
            "x0": args.x0,
            "y0": args.y0,
            "nx": args.nx,
            "ny": args.ny,
            "pixel_origin": "0-based Python pixels",
            "meaning": "potential/physical region before multigrid; ghost pixels are added automatically",
        },
        "ghost_padding": {
            "nghost_output_pixels_each_side": args.nghost,
            "input_ghost_pixels_each_side": ghost_pixels,
            "blank_value": args.blank_value,
            "padding_pixels_before_multigrid": padding,
            "boundary_input_shape_yx": list(bz_crop.shape),
            "boundary_output_shape_yx": list(bz_mg.shape),
        },
        "multigrid": {
            "level": args.level,
            "reducer": reducer,
            "shape_after_multigrid_yx": list(bz_mg.shape),
            "physical_shape_after_multigrid_yx": [args.ny // reducer, args.nx // reducer],
        },
        "preprocess": {
            "enabled": args.preprocess,
            "mu3": args.preprocess_mu3,
            "mu4": args.preprocess_mu4,
            "max_iter": args.preprocess_max_iter,
            "tol": args.preprocess_tol,
            "metrics": preprocess_metrics,
            "diagnostics_before": diagnostics_before,
            "diagnostics_after": diagnostics_after,
            "diagnostics_ratio_after_over_before": diagnostics_ratio(diagnostics_before, diagnostics_after),
        },
        "coordinate": {
            "solar_radius_cm": args.solar_radius_cm,
            "dx_cm": dx_cm,
            "dy_cm": dy_cm,
            "xc_cm": xc_cm,
            "yc_cm": yc_cm,
            "note": "x/y are local CEA-plane coordinates relative to the input patch reference pixel unless overridden.",
        },
        "amrvac": meta.as_dict,
        "grid_validation": grid_validation,
        "outputs": {
            "potential_boundary": str(output_dir / "potential_boundary.dat"),
            "nlfff_boundary": str(output_dir / "nlfff_boundary.dat"),
            "allboundaries": str(output_dir / "allboundaries.dat"),
            "grid_ini": str(output_dir / "grid.ini"),
            "nd_ini": str(output_dir / "nd.ini"),
            "parameters_txt": str(output_dir / "amrvac_parameters.txt"),
        },
    }
    (output_dir / "parameters.json").write_text(
        json.dumps(json_safe(metadata), indent=2, allow_nan=False),
        encoding="utf-8",
    )

    if args.quicklook:
        plot_quicklook(output_dir / "boundary_quicklook.png", bx_mg, by_mg, bz_mg, vmax=args.vmax)

    print(parameters_txt)
    print(f"Saved products in {output_dir}")
    return 0


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cea-dir", type=Path, help="Directory containing Br.fits, Bt.fits, and Bp.fits.")
    parser.add_argument("--br", type=Path, help="Br FITS path.")
    parser.add_argument("--bt", type=Path, help="Bt FITS path.")
    parser.add_argument("--bp", type=Path, help="Bp FITS path.")
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--x0", required=True, type=int, help="0-based starting x pixel of the potential/physical region.")
    parser.add_argument("--y0", required=True, type=int, help="0-based starting y pixel of the potential/physical region.")
    parser.add_argument("--nx", required=True, type=int, help="Potential/physical-region x size before multigrid.")
    parser.add_argument("--ny", required=True, type=int, help="Potential/physical-region y size before multigrid.")
    parser.add_argument("--level", type=int, default=1, help="Guo multigrid level: 1 keeps native pixels, 3 means 4x4 bin.")
    parser.add_argument("--nghost", type=int, default=2, help="Ghost pixels on each side used by AMRVAC.")
    parser.add_argument("--blank-value", type=float, default=0.0, help="Fill value for ghost pixels outside the FITS frame.")
    parser.add_argument("--solar-radius-cm", type=float, default=SOLAR_RADIUS_CM)
    parser.add_argument("--dx-cm", type=float, help="Override output pixel dx in cm after multigrid.")
    parser.add_argument("--dy-cm", type=float, help="Override output pixel dy in cm after multigrid.")
    parser.add_argument("--xc-cm", type=float, help="Override selected-window center x in cm.")
    parser.add_argument("--yc-cm", type=float, help="Override selected-window center y in cm.")
    parser.add_argument("--z-min-10mm", type=float, default=0.1, help="xprobmin3 in 10 Mm units.")
    parser.add_argument("--domain-nx3", type=int, help="Override vertical domain_nx3; default ny-4.")
    parser.add_argument("--block-nx1", type=int)
    parser.add_argument("--block-nx2", type=int)
    parser.add_argument("--block-nx3", type=int)
    parser.add_argument("--preprocess", dest="preprocess", action="store_true", default=True)
    parser.add_argument("--no-preprocess", dest="preprocess", action="store_false")
    parser.add_argument("--preprocess-mu3", type=float, default=DEFAULT_NLFFF_PREPROCESS_MU3)
    parser.add_argument("--preprocess-mu4", type=float, default=DEFAULT_NLFFF_PREPROCESS_MU4)
    parser.add_argument("--preprocess-max-iter", type=int, default=DEFAULT_NLFFF_PREPROCESS_MAX_ITER)
    parser.add_argument("--preprocess-report-interval", type=int, default=50)
    parser.add_argument("--preprocess-tol", type=float, default=DEFAULT_NLFFF_PREPROCESS_TOL)
    parser.add_argument("--quicklook", dest="quicklook", action="store_true", default=True)
    parser.add_argument("--no-quicklook", dest="quicklook", action="store_false")
    parser.add_argument("--vmax", type=float, default=500.0)
    return parser.parse_args(argv)


def validate_physical_window(args, reducer: int):
    if args.nx <= 0 or args.ny <= 0:
        raise ValueError("physical nx, ny before binning must be positive")
    if args.nx % reducer or args.ny % reducer:
        raise ValueError(
            f"physical nx, ny must be divisible by reducer={reducer}; got nx={args.nx}, ny={args.ny}"
        )
    if args.nghost < 0:
        raise ValueError("--nghost must be non-negative")
    nx, ny = args.nx // reducer, args.ny // reducer
    report = validate_physical_grid(nx, ny)
    nz = ny if args.domain_nx3 is None else args.domain_nx3
    if nz <= 0:
        raise ValueError("--domain-nx3 must be positive")
    for axis, size in enumerate((nx, ny, nz), 1):
        block = getattr(args, f"block_nx{axis}")
        if block is not None and (block <= 0 or block % 2 or size % block):
            raise ValueError(
                f"block_nx{axis}={block} must be positive, even and divide domain_nx{axis}={size}"
            )
    report["domain_shape_xyz"] = [nx, ny, nz]
    report["block_shape_xyz"] = [args.block_nx1, args.block_nx2, args.block_nx3]
    return report


def resolve_component_paths(args):
    if args.cea_dir is not None:
        br = args.br or args.cea_dir / "Br.fits"
        bt = args.bt or args.cea_dir / "Bt.fits"
        bp = args.bp or args.cea_dir / "Bp.fits"
    else:
        br, bt, bp = args.br, args.bt, args.bp
    missing = [name for name, path in (("Br", br), ("Bt", bt), ("Bp", bp)) if path is None or not Path(path).exists()]
    if missing:
        raise FileNotFoundError(f"missing component FITS: {', '.join(missing)}")
    return Path(br), Path(bt), Path(bp)


def read_fits_data(path: Path):
    with fits.open(path) as hdul:
        # JSOC SHARP exports use a blank primary HDU followed by a compressed
        # image extension, whereas locally produced FITS files often store the
        # image in the primary HDU.  Accept either layout.
        for hdu in hdul:
            if hdu.data is not None and np.ndim(hdu.data) == 2:
                return np.asarray(hdu.data, dtype=float), hdu.header.copy()
    raise ValueError(f"no two-dimensional image data found in {path}")


def infer_spacing_cm(header, reducer: int, args):
    if args.dx_cm is not None:
        dx_cm = float(args.dx_cm)
    else:
        cdelt1 = abs(float(header.get("CDELT1")))
        dx_cm = float(args.solar_radius_cm) * np.radians(cdelt1) * reducer
    if args.dy_cm is not None:
        dy_cm = float(args.dy_cm)
    else:
        cdelt2 = abs(float(header.get("CDELT2", header.get("CDELT1"))))
        dy_cm = float(args.solar_radius_cm) * np.radians(cdelt2) * reducer
    return dx_cm, dy_cm


def infer_window_center_cm(header, args):
    if args.xc_cm is not None:
        xc_cm = float(args.xc_cm)
    else:
        x_center_original = args.x0 + 0.5 * (args.nx - 1)
        crpix1 = float(header.get("CRPIX1", (int(header["NAXIS1"]) + 1.0) / 2.0))
        cdelt1 = float(header.get("CDELT1"))
        x_degree = (x_center_original + 1.0 - crpix1) * cdelt1
        xc_cm = float(args.solar_radius_cm) * np.radians(x_degree)

    if args.yc_cm is not None:
        yc_cm = float(args.yc_cm)
    else:
        y_center_original = args.y0 + 0.5 * (args.ny - 1)
        crpix2 = float(header.get("CRPIX2", (int(header["NAXIS2"]) + 1.0) / 2.0))
        cdelt2 = float(header.get("CDELT2", header.get("CDELT1")))
        y_degree = (y_center_original + 1.0 - crpix2) * cdelt2
        yc_cm = float(args.solar_radius_cm) * np.radians(y_degree)
    return xc_cm, yc_cm


def print_preprocess_progress(values):
    print(
        "  "
        f"it={values['iteration']:5d} "
        f"dL={values['dL']:.3e} "
        f"L12={values['L1'] + values['L2']:.3e} "
        f"L3={values['L3']:.3e} "
        f"L4={values['L4']:.3e}"
    )


def diagnostics_ratio(before, after):
    keys = ("eps_force", "eps_torque", "eps_smooth", "L1", "L2", "L4")
    ratios = {}
    for key in keys:
        denom = before.get(key)
        numer = after.get(key)
        if denom in (None, 0.0):
            ratios[key] = None
        else:
            ratios[key] = numer / denom
    return ratios


def plot_quicklook(path: Path, bx, by, bz, *, vmax: float):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4), dpi=160, constrained_layout=True)
    for ax, data, title in zip(axes, (bx, by, bz), ("Bx = Bp", "By = -Bt", "Bz = Br"), strict=True):
        im = ax.imshow(data, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax)
        ax.set_title(title)
        ax.set_xlabel("x pixel")
        ax.set_ylabel("y pixel")
        fig.colorbar(im, ax=ax, shrink=0.8)
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def json_safe(value):
    if isinstance(value, dict):
        return {key: json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(item) for item in value]
    if isinstance(value, np.generic):
        return json_safe(value.item())
    if isinstance(value, float) and not np.isfinite(value):
        return None
    return value


if __name__ == "__main__":
    raise SystemExit(main())

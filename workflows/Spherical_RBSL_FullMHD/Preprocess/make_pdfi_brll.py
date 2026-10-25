#!/usr/bin/env python3
"""Create the PDFI_SS Br(theta, phi) FITS file from workflow products."""

from __future__ import annotations

import argparse
import json
import struct
from pathlib import Path

import numpy as np
from astropy.io import fits


def read_cutout(path: Path) -> np.ndarray:
    with path.open("rb") as handle:
        ntheta, nphi = struct.unpack("<ii", handle.read(8))
        values = np.fromfile(handle, dtype="<f8")
    expected = 3 * ntheta * nphi
    if values.size != expected:
        raise ValueError(f"{path}: expected {expected} values, found {values.size}")
    return values.reshape(3, nphi, ntheta)


def read_pfss_boundary(path: Path) -> np.ndarray:
    with path.open("rb") as handle:
        nphi, ntheta = struct.unpack("<ii", handle.read(8))
        np.fromfile(handle, dtype="<f8", count=ntheta)
        np.fromfile(handle, dtype="<f8", count=nphi)
        values = np.fromfile(handle, dtype="<f8", count=ntheta * nphi)
    if values.size != ntheta * nphi:
        raise ValueError(f"{path}: incomplete radial-field array")
    return values.reshape(ntheta, nphi)


def build_brll(product_dir: Path, output: Path) -> None:
    parameters = json.loads((product_dir / "parameters.json").read_text(encoding="utf-8"))
    grid = parameters["grid"]
    nghost = int(grid["effective_nghost"])
    nlat = int(grid["nlat"])
    nlon = int(grid["nlon"])
    x0 = int(grid["x0"])
    y0 = int(grid["y0"])

    cutout = read_cutout(product_dir / "cutout.dat")
    br_phi_theta = cutout[0]
    if nghost:
        brll = br_phi_theta[nghost:-nghost, nghost:-nghost]
    else:
        brll = br_phi_theta
    if brll.shape != (nlon, nlat):
        raise ValueError(f"physical Br shape {brll.shape}, expected {(nlon, nlat)}")

    global_br = read_pfss_boundary(product_dir / "pfss_boundary.dat")
    global_patch = global_br[y0 : y0 + nlat, x0 : x0 + nlon]
    mismatch = global_patch - brll.T[::-1, :]
    max_abs = float(np.max(np.abs(mismatch)))
    if max_abs != 0.0:
        raise ValueError(
            "cutout Br and the patch inserted into pfss_boundary.dat differ: "
            f"max_abs={max_abs}"
        )

    header = fits.Header()
    header["BUNIT"] = "G"
    header["CONTENT"] = "Br(theta,phi) for PDFI_SS"
    header["NTHETA"] = nlat
    header["NPHI"] = nlon
    header["NGHOST"] = 0
    header["LONMIN"] = float(grid["lon_range"][0])
    header["LONMAX"] = float(grid["lon_range"][1])
    header["LATMIN"] = float(grid["lat_range"][0])
    header["LATMAX"] = float(grid["lat_range"][1])
    output.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(brll.astype(np.float32), header=header).writeto(output, overwrite=True)

    print(f"wrote: {output}")
    print(f"Python FITS shape (phi,theta): {brll.shape}")
    print(f"Fortran Brll(theta,phi): ({nlat},{nlon})")
    print(f"Br range: {brll.min():.6f} .. {brll.max():.6f} G")
    print("ghost cells: 0")
    print("global PFSS patch comparison: exact")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "product_dir",
        type=Path,
        help="Directory containing cutout.dat, pfss_boundary.dat, and parameters.json",
    )
    parser.add_argument(
        "--output",
        type=Path,
        help="Output FITS path; default is PRODUCT_DIR/Brll0000.fits",
    )
    args = parser.parse_args()
    product_dir = args.product_dir.expanduser().resolve()
    output = args.output.expanduser().resolve() if args.output else product_dir / "Brll0000.fits"
    build_brll(product_dir, output)


if __name__ == "__main__":
    main()

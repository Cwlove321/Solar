#!/usr/bin/env python3
"""Compare simesh stretched-DAT coordinates with an AMRVAC VTU export.

Example (on the host where rbsl14 is mounted):
  PYTHONPATH=src python scripts/validate_stretched_dat_vtu.py \
    --dat data0000.dat --vtu data0000.vtu --par amrvac.par --spherical
"""
from __future__ import annotations

import argparse
import re
import struct
from functools import lru_cache
import numpy as np

from simesh.amrvac import open_dataset


def spherical(points):
    r = np.linalg.norm(points, axis=1)
    theta = np.arccos(np.clip(points[:, 2] / r, -1.0, 1.0))
    phi = np.mod(np.arctan2(points[:, 1], points[:, 0]), 2.0 * np.pi)
    return np.column_stack((r, theta, phi))


@lru_cache(maxsize=4)
def _raw_vtu_layout(path):
    with open(path, "rb") as fh:
        prefix = fh.read(32 * 1024 * 1024)
    marker = b"<AppendedData"
    header_end = prefix.find(marker)
    if header_end < 0:
        raise ValueError("VTU has no AppendedData section in its XML header")
    header = prefix[:header_end].decode("ascii")
    raw_start = header_end + prefix[header_end:].find(b"_") + 1
    pieces = re.findall(r"<Piece\b[^>]*>(.*?)</Piece>", header, flags=re.S)
    return raw_start, pieces


def read_amrvac_raw_vtu_piece(path, piece_index):
    """Read one AMRVAC raw-appended VTU piece without loading its multi-GB file.

    MPI-AMRVAC emits raw XML appended data; several VTK builds reject that XML
    form even though the blocks are valid.  The offsets are byte offsets from
    the underscore directly after ``<AppendedData ...>``.
    """
    raw_start, pieces = _raw_vtu_layout(path)
    if piece_index >= len(pieces):
        raise IndexError(f"VTU has {len(pieces)} pieces, requested {piece_index}")
    piece = pieces[piece_index]
    point_match = re.search(r"<Points>\s*<DataArray\b[^>]*offset=\"(\d+)\"", piece, flags=re.S)
    connectivity_match = re.search(r"Name=\"connectivity\"[^>]*offset=\"(\d+)\"", piece)
    if not point_match or not connectivity_match:
        raise ValueError("VTU piece lacks Points or connectivity appended offsets")

    def block(offset, dtype):
        with open(path, "rb") as fh:
            fh.seek(raw_start + int(offset))
            nbytes = struct.unpack("<I", fh.read(4))[0]
            return np.frombuffer(fh.read(nbytes), dtype=dtype).copy()

    points = block(point_match.group(1), "<f4").reshape(-1, 3).astype(float)
    connectivity = block(connectivity_match.group(1), "<i4").reshape(-1, 8)
    return points, connectivity, len(pieces)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dat", required=True)
    parser.add_argument("--vtu", required=True)
    parser.add_argument("--par", required=True, help="MPI-AMRVAC par file containing meshlist stretching")
    parser.add_argument("--samples", type=int, default=100)
    parser.add_argument("--seed", type=int, default=20241009)
    parser.add_argument("--spherical", action="store_true", help="VTU points are Cartesian images of (r, theta, phi)")
    args = parser.parse_args()

    ds = open_dataset(args.dat, grid_config=args.par)
    if not ds.uses_stretched_coordinates:
        raise SystemExit("The supplied par file does not describe a stretched grid.")
    rng = np.random.default_rng(args.seed)
    leaves = rng.choice(int(ds.nleafs), size=min(args.samples, int(ds.nleafs)), replace=False)
    nodes = ds.block_node_coordinates
    centers = ds.block_center_coordinates
    radial_errors = []
    corner_errors = []
    center_errors = []

    # VTU stores Cartesian corner points, not AMR block ids. Select its points
    # by each reconstructed leaf's physical interval and compare the sets. The
    # AMRVAC VTU writer emits pieces in the same SFC leaf order as DAT blocks.
    for leaf in leaves:
        cart_points, connectivity, npieces = read_amrvac_raw_vtu_piece(args.vtu, int(leaf))
        points = spherical(cart_points) if args.spherical else cart_points[:, :int(ds.ndim)]
        # A Cartesian average of spherical-cell vertices is not the midpoint in
        # (r, theta, phi). Transform each vertex first, then average the native
        # coordinate components just as AMRVAC defines a cell center.
        vtu_vertices = cart_points[connectivity]
        vtu_cell_centers = spherical(vtu_vertices.reshape(-1, 3)).reshape(vtu_vertices.shape) if args.spherical else vtu_vertices[:, :, :int(ds.ndim)]
        vtu_cell_centers = vtu_cell_centers.mean(axis=1)
        lo = np.array([axis[leaf, 0] for axis in nodes])
        hi = np.array([axis[leaf, -1] for axis in nodes])
        actual = points
        expected = np.array(np.meshgrid(*[axis[leaf] for axis in nodes], indexing="ij")).reshape(int(ds.ndim), -1).T
        # VTU can duplicate shared corners. Compare each expected corner with
        # its nearest VTU point; this is insensitive to VTK cell ordering.
        distances = np.sqrt(((expected[:, None, :] - actual[None, :, :]) ** 2).sum(axis=2))
        corner_errors.extend(distances.min(axis=1))
        expected_centers = np.array(np.meshgrid(*[axis[leaf] for axis in centers], indexing="ij")).reshape(int(ds.ndim), -1).T
        # Cell centers are not necessarily VTU points, so use the same VTU
        # cells' geometric centers when available.
        if len(vtu_cell_centers):
            vtu_centers = vtu_cell_centers
            distances = np.sqrt(((expected_centers[:, None, :] - vtu_centers[None, :, :]) ** 2).sum(axis=2))
            center_errors.extend(distances.min(axis=1))
            radial_errors.extend(np.abs(expected_centers[:, 0, None] - vtu_centers[None, :, 0]).min(axis=1))

    print(f"sampled leaf blocks: {len(leaves)}")
    print(f"VTU pieces: {npieces}")
    print(f"corner max error: {np.max(corner_errors):.16e}")
    print(f"center max error: {np.max(center_errors):.16e}")
    print(f"radial max error: {np.max(radial_errors):.16e}")
    print(f"radial RMS error: {np.sqrt(np.mean(np.square(radial_errors))):.16e}")
    print(f"bottom-level minimum dr: {min(np.diff(axis).min() for axis in nodes[0]):.16e}")


if __name__ == "__main__":
    main()

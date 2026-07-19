#!/usr/bin/env python3
from __future__ import annotations

import argparse
import os
from collections import namedtuple
from glob import glob
from pathlib import Path

import numpy as np
from pyevtk.hl import gridToVTK
from simesh.amrvac import open_dataset, read_uniform


def init_boundary(filename, B1, B2, B3):
    # Arrays are stored as z,y,x for compatibility with the old script.
    ax, bx = -B1[:, :, 0], B1[:, :, -1]
    ay, by = -B2[:, 0, :], B2[:, -1, :]
    az, bz = -B3[0, :, :], B3[-1, :, :]
    total = 0.0
    with open(filename, "wb") as f:
        for arr in (ax, bx, ay, by, az, bz):
            arr.tofile(f)
            total += float(np.sum(arr))
    return total


def write_binary_data(file, dx, B1, B2, B3):
    nz, ny, nx = B1.shape
    with open(file, "wb") as f:
        np.array(dx, dtype=np.double).tofile(f)
        np.array(nx, dtype=np.int32).tofile(f)
        np.array(ny, dtype=np.int32).tofile(f)
        np.array(nz, dtype=np.int32).tofile(f)
        B1.tofile(f)
        B2.tofile(f)
        B3.tofile(f)


def curl(b1, b2, b3, dx=1, dy=1, dz=1, dim="zyx"):
    if dim == "zyx":
        curl_x = np.gradient(b3, axis=1) / dy - np.gradient(b2, axis=0) / dz
        curl_y = np.gradient(b1, axis=0) / dz - np.gradient(b3, axis=2) / dx
        curl_z = np.gradient(b2, axis=2) / dx - np.gradient(b1, axis=1) / dy
    else:
        curl_x = np.gradient(b3, axis=1) / dy - np.gradient(b2, axis=2) / dz
        curl_y = np.gradient(b1, axis=2) / dz - np.gradient(b3, axis=0) / dx
        curl_z = np.gradient(b2, axis=0) / dx - np.gradient(b1, axis=1) / dy
    return curl_x, curl_y, curl_z


def process_data(path, out_path, start, end, level=2, save_nlfff=False, save_bd=False, save_vtk=False, save_npz=False):
    out_path = Path(out_path)
    for name in ["npz", "vtu", "A", "potential", "bd", "nlfff"]:
        (out_path / name).mkdir(parents=True, exist_ok=True)

    file_dat = sorted(glob(os.path.join(path, "*.dat")))
    if end > len(file_dat):
        end = len(file_dat)

    for i in range(start, end):
        dat = file_dat[i]
        ds = open_dataset(dat)
        domain_nx = np.asarray(ds.domain_nx, dtype=int)
        nx, ny, nz = (level * domain_nx).tolist()
        xmin, ymin, zmin = np.asarray(ds.physical_domain[0], dtype=float).tolist()
        xmax, ymax, zmax = np.asarray(ds.physical_domain[1], dtype=float).tolist()
        dx = (xmax - xmin) / nx

        names = list(ds.wnames)
        W = namedtuple("W", names)
        print(f"[{i}] {dat} resolution=({nx},{ny},{nz})", flush=True)

        udata = read_uniform(
            dat,
            resolution=(nx, ny, nz),
            bounds=(np.array([xmin, ymin, zmin]), np.array([xmax, ymax, zmax])),
            field_indices=list(range(len(names))),
            interpolation="zero",
        )
        # New simesh returns nx,ny,nz,nvar. Keep old script convention nvar,nz,ny,nx.
        slab_all = np.transpose(udata, (3, 2, 1, 0))
        del udata
        w = W(*slab_all)
        num = f"{i:04d}"

        if save_nlfff:
            file = out_path / "nlfff" / f"nlfff_{num}.binary"
            write_binary_data(file, dx, w.b1, w.b2, w.b3)

        if save_bd:
            file = out_path / "bd" / f"bd_{num}.binary"
            init_boundary(file, w.b1, w.b2, w.b3)

        j1, j2, j3 = curl(w.b1, w.b2, w.b3, dx=dx, dy=dx, dz=dx, dim="zyx")
        j = np.sqrt(j1**2 + j2**2 + j3**2)
        v1 = w.m1 / w.rho * 116.45
        v2 = w.m2 / w.rho * 116.45
        v3 = w.m3 / w.rho * 116.45

        if save_vtk:
            x = np.linspace(xmin, xmax, nx)
            y = np.linspace(ymin, ymax, ny)
            z = np.linspace(zmin, zmax, nz)
            gridToVTK(
                str(out_path / "vtu" / f"eru{num}"),
                x,
                y,
                z,
                pointData={
                    "b1": np.ascontiguousarray(w.b1.T),
                    "b2": np.ascontiguousarray(w.b2.T),
                    "b3": np.ascontiguousarray(w.b3.T),
                    "j": np.ascontiguousarray(j.T),
                    "e": np.ascontiguousarray(w.e.T),
                    "rho": np.ascontiguousarray(w.rho.T),
                    "v1": np.ascontiguousarray(v1.T),
                    "v2": np.ascontiguousarray(v2.T),
                    "v3": np.ascontiguousarray(v3.T),
                },
            )

        if save_npz:
            data_dict = w._asdict()
            exclude = {"b1", "b2", "b3"}
            filtered = {k: v for k, v in data_dict.items() if k not in exclude}
            np.savez_compressed(out_path / "npz" / f"eru{num}.npz", **filtered)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process AMRVAC dat files with new simesh.amrvac API.")
    parser.add_argument("path", type=str, help="Path to input data files.")
    parser.add_argument("out_path", type=str, help="Path to output directory.")
    parser.add_argument("start", type=int, help="Start index.")
    parser.add_argument("end", type=int, help="End index.")
    parser.add_argument("--level", type=int, default=2, help="Resolution multiplier (nx = level * domain_nx).")
    parser.add_argument("--save_nlfff", action="store_true", help="Save NLFFF binary files.")
    parser.add_argument("--save_bd", action="store_true", help="Save boundary binary files.")
    parser.add_argument("--save_vtk", action="store_true", help="Save VTK files.")
    parser.add_argument("--save_npz", action="store_true", help="Save other variables as NPZ files.")
    args = parser.parse_args()
    process_data(args.path, args.out_path, args.start, args.end, args.level, args.save_nlfff, args.save_bd, args.save_vtk, args.save_npz)

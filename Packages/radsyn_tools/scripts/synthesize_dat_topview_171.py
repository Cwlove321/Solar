#!/usr/bin/env python3
"""Create an AIA 171 Å top view from a stretched MPI-AMRVAC DAT file.

The DAT is sampled onto an axis-aligned *uniform Cartesian* array before the
pure-Python Method-1 DDA radiative transfer is applied.  This deliberately
keeps AMR/stretch-coordinate handling inside simesh.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "DDA_method1"))
from RST_cart.DDA import VolumeContainer, integrate_plane_points
from RST_cart.generate_plane_points import generate_plane_points
from RST_cart.physical_module import cal_euv
from simesh.amrvac import open_dataset

RSUN_CM = 6.955e10
NUMBER_DENSITY_UNIT = 1.0e9
TEMPERATURE_UNIT = 1.0e6


def sphere_to_cartesian(spherical: np.ndarray) -> np.ndarray:
    r, theta, phi = spherical.T
    return np.column_stack((
        r * np.sin(theta) * np.cos(phi),
        r * np.sin(theta) * np.sin(phi),
        r * np.cos(theta),
    ))


def cartesian_to_sphere(cartesian: np.ndarray) -> np.ndarray:
    r = np.linalg.norm(cartesian, axis=1)
    theta = np.arccos(np.clip(cartesian[:, 2] / r, -1.0, 1.0))
    phi = np.mod(np.arctan2(cartesian[:, 1], cartesian[:, 0]), 2.0 * np.pi)
    return np.column_stack((r, theta, phi))


def cartesian_bounds(xmin: np.ndarray, xmax: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    corners = np.array(np.meshgrid(*zip(xmin, xmax), indexing="ij")).reshape(3, -1).T
    xyz = sphere_to_cartesian(corners)
    return xyz.min(axis=0), xyz.max(axis=0)


def parse_resolution(text: str) -> tuple[int, int, int]:
    value = tuple(int(item) for item in text.split(","))
    if len(value) != 3 or min(value) < 2:
        raise argparse.ArgumentTypeError("resolution must be nx,ny,nz with each value >= 2")
    return value


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dat", type=Path)
    parser.add_argument("par", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--resolution", type=parse_resolution, default=(96, 72, 48))
    parser.add_argument("--image-size", type=int, default=192)
    parser.add_argument("--wave", type=int, default=171, choices=(94, 131, 171, 193, 211, 304, 335))
    parser.add_argument(
        "--view", choices=("bottom-inward", "south-to-north-tangent"),
        default="bottom-inward",
    )
    parser.add_argument(
        "--axis-path", type=Path,
        help="path.dat whose maximum-radius point defines the rope centre for the tangent view",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    dataset = open_dataset(str(args.dat), grid_config=str(args.par))
    dataset.load_data()
    names = list(dataset.metadata["w_names"])
    original_indices = [names.index(name) for name in ("rho", "e")]
    param_names = list(dataset.metadata["param_names"])
    gamma = float(dataset.metadata["params"][param_names.index("gamma")])
    spherical_left, spherical_right = (np.asarray(edge, dtype=float) for edge in dataset.physical_domain)
    cart_left, cart_right = cartesian_bounds(spherical_left, spherical_right)
    padding = 1.0e-10 * np.maximum(cart_right - cart_left, 1.0)
    cart_left -= padding
    cart_right += padding

    grid_shape = args.resolution
    axes = [
        np.linspace(cart_left[i], cart_right[i], grid_shape[i], endpoint=False)
        + (cart_right[i] - cart_left[i]) / (2.0 * grid_shape[i])
        for i in range(3)
    ]
    xx, yy, zz = np.meshgrid(*axes, indexing="ij")
    sample_xyz = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    sample_spherical = cartesian_to_sphere(sample_xyz)
    valid = np.all(sample_spherical >= spherical_left, axis=1) & np.all(sample_spherical <= spherical_right, axis=1)
    values = np.zeros((sample_xyz.shape[0], 2), dtype=float)
    # This is the only DAT-to-uniform interpolation step.  simesh chooses the
    # AMR leaf from true stretched physical bounds, then linearly samples it.
    values[valid] = dataset.sample_points(
        sample_spherical[valid], field_indices=original_indices, interpolation="linear"
    )
    rho_code = values[:, 0].reshape(grid_shape)
    internal_energy = values[:, 1].reshape(grid_shape)
    temperature = (gamma - 1.0) * internal_energy / np.maximum(rho_code, 1.0e-30) * TEMPERATURE_UNIT
    temperature[~np.isfinite(temperature) | (rho_code <= 0.0)] = 1.0
    electron_density = rho_code * NUMBER_DENSITY_UNIT
    np.save(args.output / "rho_cm3_uniform_cartesian.npy", electron_density)
    np.save(args.output / "temperature_K_uniform_cartesian.npy", temperature)

    # Bottom-boundary centre -> sphere centre is the default inward view.
    bottom_spherical = np.array([[spherical_left[0], *(0.5 * (spherical_left[1:] + spherical_right[1:]))]])
    bottom_center = sphere_to_cartesian(bottom_spherical)[0]
    view_centre = bottom_center
    view_description = "bottom-boundary centre to sphere centre"
    rope_spherical = None
    if args.view == "bottom-inward":
        view_direction = -bottom_center / np.linalg.norm(bottom_center)
    else:
        if args.axis_path is None:
            parser.error("--axis-path is required for --view south-to-north-tangent")
        axis_path = np.loadtxt(args.axis_path)
        rope = axis_path[int(np.nanargmax(axis_path[:, 2])), :3]
        # path.dat columns are (longitude deg, latitude deg, radius R_sun).
        rope_spherical = np.array([[rope[2], np.radians(90.0 - rope[1]), np.radians(rope[0])]])
        view_centre = sphere_to_cartesian(rope_spherical)[0]
        theta, phi = rope_spherical[0, 1:]
        e_theta = np.array((np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta)))
        # -e_theta raises latitude: an observer/ray travelling from south to north.
        view_direction = -e_theta
        view_description = "south to north along the tangent plane at the rope centre"
    field_width = 1.10 * max(cart_right - cart_left) * RSUN_CM
    plane_points = generate_plane_points(
        view_centre * RSUN_CM, view_direction, field_width, args.image_size
    )
    kappa, emissivity = cal_euv(
        electron_density, electron_density * temperature, temperature, args.wave, absorption=False
    )
    volume = VolumeContainer(kappa, emissivity, cart_left * RSUN_CM, cart_right * RSUN_CM)
    image = integrate_plane_points(plane_points, view_direction, volume).reshape(args.image_size, args.image_size)
    np.save(args.output / f"aia_{args.wave}_topview_linear.npy", image)
    positive = image[image > 0]
    log_image = np.where(image > 0, np.log10(image), np.nan)
    vmin, vmax = np.percentile(np.log10(positive), (1, 99.7)) if len(positive) else (0.0, 1.0)
    fig, ax = plt.subplots(figsize=(7.2, 6.2), constrained_layout=True)
    # Matplotlib does not ship the SunPy-specific ``sdoaia171`` map, so use a
    # portable AIA-like warm palette for a directly runnable script.
    shown = ax.imshow(log_image.T, origin="lower", cmap="afmhot", vmin=vmin, vmax=vmax)
    ax.set_title(f"Synthetic AIA {args.wave} Å, {view_description}")
    ax.set_xlabel("Image-plane pixel")
    ax.set_ylabel("Image-plane pixel")
    fig.colorbar(shown, ax=ax, label=r"$\log_{10}$ intensity (arbitrary units)")
    fig.savefig(args.output / f"aia_{args.wave}_topview.png", dpi=220)
    plt.close(fig)
    (args.output / "metadata.json").write_text(json.dumps({
        "dat": str(args.dat), "par": str(args.par), "wave_A": args.wave,
        "uniform_cartesian_shape": grid_shape,
        "coordinate_source": "simesh true stretched AMR coordinates",
        "view": args.view,
        "view_description": view_description,
        "view_direction": view_direction.tolist(),
        "bottom_boundary_center_spherical": bottom_spherical[0].tolist(),
        "rope_center_spherical": None if rope_spherical is None else rope_spherical[0].tolist(),
        "cartesian_domain_code_units": [cart_left.tolist(), cart_right.tolist()],
    }, indent=2) + "\n")
    print(args.output / f"aia_{args.wave}_topview.png")


if __name__ == "__main__":
    main()

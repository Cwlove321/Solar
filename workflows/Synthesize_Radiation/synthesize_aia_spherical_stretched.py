#!/usr/bin/env python3
"""Synthesize optically thin and absorbing AIA 304 views from spherical AMRVAC.

The input adapter reconstructs a stretched ``(r, theta, phi)`` AMR mesh with
``simesh`` and resamples it to a uniform Cartesian volume. Parallel physical
lines of sight are then integrated with centimetre-unit ray-box intersection
and DDA traversal. The absorbing result is a continuum dark-cloud experiment,
not an NLTE He II 304 line calculation.
"""

from __future__ import annotations

import argparse
import gc
import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

RSUN_CM = 6.955e10
NUMBER_DENSITY_UNIT = 1.0e9
TEMPERATURE_UNIT = 1.0e6
VIEW_NAMES = ("bottom-inward", "south-to-north-tangent")
DISK_TOP_RADIUS_RSUN = 1.005
HYDROGEN_EDGE_A = 911.753
HELIUM_I_EDGE_A = 504.259
SIGMA_H_I_EDGE_CM2 = 6.30e-18
SIGMA_HE_I_EDGE_CM2 = 7.83e-18
HELIUM_ABUNDANCE = 0.10
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DDA_ROOT = REPO_ROOT / "Packages" / "radsyn_tools" / "DDA_method1"


def sphere_to_cartesian(spherical: np.ndarray) -> np.ndarray:
    r, theta, phi = spherical.T
    return np.column_stack(
        (
            r * np.sin(theta) * np.cos(phi),
            r * np.sin(theta) * np.sin(phi),
            r * np.cos(theta),
        )
    )


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


def view_geometry(
    view: str,
    bottom_center: np.ndarray,
    rope_spherical: np.ndarray,
    outer_radius: float,
    outer_margin: float,
) -> tuple[np.ndarray, np.ndarray, str]:
    if view == "bottom-inward":
        radial = bottom_center / np.linalg.norm(bottom_center)
        centre = radial * (outer_radius + outer_margin)
        return centre, -radial, "Radial top-down"

    centre = sphere_to_cartesian(rope_spherical)[0]
    theta, phi = rope_spherical[0, 1:]
    e_theta = np.array(
        (np.cos(theta) * np.cos(phi), np.cos(theta) * np.sin(phi), -np.sin(theta))
    )
    return centre, -e_theta, "South-to-north tangent"


def finite_log_image(image: np.ndarray, floor: float) -> np.ndarray:
    result = np.full(image.shape, floor, dtype=float)
    valid = np.isfinite(image) & (image > 0.0)
    result[valid] = np.log10(image[valid])
    return result


def image_stats(image: np.ndarray) -> dict[str, float | int]:
    positive = image[np.isfinite(image) & (image > 0.0)]
    if positive.size == 0:
        return {"positive_pixels": 0}
    log_values = np.log10(positive)
    return {
        "positive_pixels": int(positive.size),
        "linear_min": float(positive.min()),
        "linear_max": float(positive.max()),
        "log10_p01": float(np.percentile(log_values, 1.0)),
        "log10_p50": float(np.percentile(log_values, 50.0)),
        "log10_p997": float(np.percentile(log_values, 99.7)),
    }


def transfer_sample(vc, _v_pos, _v_dir, enter_t, exit_t, cur_ind, ray) -> None:
    """Apply the exact constant-coefficient transfer solution in one cell."""
    ds = exit_t - enter_t
    if ds <= 0.0:
        return
    ix, iy, iz = cur_ind
    if not (0 <= ix < vc.dims[0] and 0 <= iy < vc.dims[1] and 0 <= iz < vc.dims[2]):
        return

    opacity = max(float(vc.kappa_view[ix, iy, iz]), 0.0)
    emissivity = max(float(vc.j_view[ix, iy, iz]), 0.0)
    cell_tau = opacity * ds
    if cell_tau < 1.0e-8:
        emitted = emissivity * ds
    else:
        emitted = emissivity * (-np.expm1(-cell_tau)) / opacity
    ray.sum_val += np.exp(-min(ray.tau, 745.0)) * emitted
    ray.tau += cell_tau


def ray_box_interval(start, direction, left_edge, right_edge):
    """Return the forward ray interval inside an axis-aligned box in cm."""
    enter = 0.0
    leave = np.inf
    for axis in range(3):
        component = float(direction[axis])
        if abs(component) < 1.0e-15:
            if start[axis] < left_edge[axis] or start[axis] > right_edge[axis]:
                return None
            continue
        first = (left_edge[axis] - start[axis]) / component
        second = (right_edge[axis] - start[axis]) / component
        enter = max(enter, min(first, second))
        leave = min(leave, max(first, second))
        if leave <= enter:
            return None
    return enter, leave


def walk_volume_cm(volume, start, direction, ray) -> bool:
    """Traverse one ray through the Cartesian volume using centimetre distances."""
    interval = ray_box_interval(
        start, direction, volume.left_edge, volume.right_edge
    )
    if interval is None:
        return False
    enter_t, leave_t = interval
    if not np.isfinite(leave_t) or leave_t <= enter_t:
        return False

    probe_offset = min(
        1.0e-9 * float(np.min(volume.dds)), 0.25 * (leave_t - enter_t)
    )
    probe = start + (enter_t + probe_offset) * direction
    cell = np.floor((probe - volume.left_edge) * volume.idds).astype(int)
    cell = np.clip(cell, 0, np.asarray(volume.dims) - 1)
    step = np.sign(direction).astype(int)
    next_t = np.full(3, np.inf, dtype=float)
    delta_t = np.full(3, np.inf, dtype=float)

    for axis in range(3):
        if step[axis] == 0:
            continue
        boundary_index = cell[axis] + (1 if step[axis] > 0 else 0)
        boundary = volume.left_edge[axis] + boundary_index * volume.dds[axis]
        next_t[axis] = (boundary - start[axis]) / direction[axis]
        delta_t[axis] = volume.dds[axis] / abs(direction[axis])

    current_t = enter_t
    while current_t < leave_t and np.all((cell >= 0) & (cell < volume.dims)):
        segment_end = min(leave_t, float(np.min(next_t)))
        if segment_end > current_t:
            transfer_sample(
                volume, start, direction, current_t, segment_end, cell, ray
            )
        if segment_end >= leave_t:
            break

        tolerance = 1.0e-12 * max(abs(segment_end), 1.0)
        crossed = np.flatnonzero(next_t <= segment_end + tolerance)
        for axis in crossed:
            cell[axis] += step[axis]
            next_t[axis] += delta_t[axis]
        current_t = segment_end
    return True


def integrate_with_tau(dda, plane_points, direction, volume):
    """Return emergent intensity and total optical depth for each ray."""
    unit_direction = np.asarray(direction, dtype=float)
    unit_direction /= np.linalg.norm(unit_direction)
    intensity = np.zeros(plane_points.shape[0], dtype=float)
    optical_depth = np.zeros_like(intensity)
    for index, point in enumerate(plane_points):
        ray = dda.RayData()
        walk_volume_cm(volume, point, unit_direction, ray)
        intensity[index] = ray.sum_val
        optical_depth[index] = ray.tau
    return intensity, optical_depth


def ratio_stats(absorbed: np.ndarray, thin: np.ndarray, tau: np.ndarray) -> dict[str, float | int]:
    valid = np.isfinite(thin) & np.isfinite(absorbed) & (thin > 0.0)
    thick = valid & np.isfinite(tau) & (tau >= 1.0)
    result: dict[str, float | int] = {
        "tau_max": float(np.nanmax(tau)),
        "tau_p99": float(np.nanpercentile(tau, 99.0)),
        "pixels_tau_ge_1": int(np.count_nonzero(thick)),
    }
    if np.any(valid):
        ratio = absorbed[valid] / thin[valid]
        result["absorbed_to_thin_p50"] = float(np.nanpercentile(ratio, 50.0))
    if np.any(thick):
        ratio = absorbed[thick] / thin[thick]
        result["absorbed_to_thin_tau_ge_1_p50"] = float(np.nanpercentile(ratio, 50.0))
    return result


def continuum_opacity_304(hydrogen_density: np.ndarray, temperature: np.ndarray) -> np.ndarray:
    """Approximate cool-plasma H I and He I bound-free opacity at 304 A.

    The MHD density unit is 1.4 m_H times 1e9 hydrogen nuclei cm^-3, so the
    sampled code density directly gives total hydrogen nuclei after multiplying
    by NUMBER_DENSITY_UNIT.  The temperature taper is an explicit upper-bound
    cloud assumption because the MHD state has no ion fractions.
    """
    wavelength = 304.0
    sigma_h_i = SIGMA_H_I_EDGE_CM2 * (wavelength / HYDROGEN_EDGE_A) ** 3
    sigma_he_i = SIGMA_HE_I_EDGE_CM2 * (wavelength / HELIUM_I_EDGE_A) ** 2
    neutral_fraction = np.clip((1.0e5 - temperature) / (1.0e5 - 2.0e4), 0.0, 1.0)
    return hydrogen_density * neutral_fraction * (
        sigma_h_i + HELIUM_ABUNDANCE * sigma_he_i
    )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("dat", type=Path)
    parser.add_argument("par", type=Path)
    parser.add_argument("axis_path", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--dda-root", type=Path, default=DEFAULT_DDA_ROOT)
    parser.add_argument("--frame", required=True)
    parser.add_argument("--resolution", type=parse_resolution, default=(160, 120, 80))
    parser.add_argument("--image-size", type=int, default=320)
    parser.add_argument("--outer-margin-rsun", type=float, default=0.05)
    parser.add_argument("--wave", type=int, default=304, choices=(304,))
    parser.add_argument(
        "--volume-output",
        type=Path,
        help="Save the reusable uniform Cartesian radiation volume as NPY files",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    sys.path.insert(0, str(args.dda_root))
    import RST_cart.DDA as dda
    from RST_cart.generate_plane_points import generate_plane_points
    from RST_cart.physical_module import cal_euv
    from simesh.amrvac import open_dataset

    dataset = open_dataset(str(args.dat), grid_config=str(args.par))
    dataset.load_data()
    names = list(dataset.metadata["w_names"])
    original_indices = [names.index(name) for name in ("rho", "e")]
    param_names = list(dataset.metadata["param_names"])
    gamma = float(dataset.metadata["params"][param_names.index("gamma")])
    time_code = float(dataset.metadata["time"])
    spherical_left, spherical_right = (
        np.asarray(edge, dtype=float) for edge in dataset.physical_domain
    )
    cart_left, cart_right = cartesian_bounds(spherical_left, spherical_right)
    padding = 1.0e-10 * np.maximum(cart_right - cart_left, 1.0)
    cart_left -= padding
    cart_right += padding

    axes = [
        np.linspace(cart_left[i], cart_right[i], args.resolution[i], endpoint=False)
        + (cart_right[i] - cart_left[i]) / (2.0 * args.resolution[i])
        for i in range(3)
    ]
    xx, yy, zz = np.meshgrid(*axes, indexing="ij")
    sample_xyz = np.column_stack((xx.ravel(), yy.ravel(), zz.ravel()))
    sample_spherical = cartesian_to_sphere(sample_xyz)
    valid = np.all(sample_spherical >= spherical_left, axis=1) & np.all(
        sample_spherical <= spherical_right, axis=1
    )
    sample_radius = np.linalg.norm(sample_xyz, axis=1).reshape(args.resolution)
    values = np.zeros((sample_xyz.shape[0], 2), dtype=float)
    values[valid] = dataset.sample_points(
        sample_spherical[valid], field_indices=original_indices, interpolation="linear"
    )
    del dataset, sample_xyz, sample_spherical, xx, yy, zz
    gc.collect()

    rho_code = values[:, 0].reshape(args.resolution)
    internal_energy = values[:, 1].reshape(args.resolution)
    temperature = (
        (gamma - 1.0)
        * internal_energy
        / np.maximum(rho_code, 1.0e-30)
        * TEMPERATURE_UNIT
    )
    temperature[~np.isfinite(temperature) | (rho_code <= 0.0)] = 1.0
    hydrogen_density = np.maximum(rho_code, 0.0) * NUMBER_DENSITY_UNIT
    electron_density = hydrogen_density

    _, emissivity = cal_euv(
        electron_density,
        electron_density * temperature,
        temperature,
        args.wave,
        absorption=False,
    )
    kappa = continuum_opacity_304(hydrogen_density, temperature)
    cloud_mask = (sample_radius >= DISK_TOP_RADIUS_RSUN) & (kappa > 0.0)
    kappa = np.where(cloud_mask, kappa, 0.0)
    kappa = np.nan_to_num(kappa, nan=0.0, posinf=0.0, neginf=0.0)
    kappa = np.maximum(kappa, 0.0)
    emissivity = np.nan_to_num(emissivity, nan=0.0, posinf=0.0, neginf=0.0)
    emissivity = np.maximum(emissivity, 0.0)
    cloud_emissivity = np.where(cloud_mask, 0.0, emissivity)

    volume_files: dict[str, str] = {}
    if args.volume_output is not None:
        args.volume_output.mkdir(parents=True, exist_ok=True)
        volume_arrays = {
            "x_centres_rsun": axes[0],
            "y_centres_rsun": axes[1],
            "z_centres_rsun": axes[2],
            "hydrogen_density_cm3": hydrogen_density,
            "temperature_k": temperature,
            f"aia{args.wave}_emissivity": emissivity,
            f"opacity_{args.wave}_cm-1": kappa,
            "valid_domain_mask": valid.reshape(args.resolution),
            "cool_cloud_mask": cloud_mask,
        }
        for name, array in volume_arrays.items():
            path = args.volume_output / f"{name}.npy"
            if array.dtype == np.bool_:
                saved = np.asarray(array, dtype=np.bool_)
            else:
                saved = np.asarray(array, dtype=np.float32)
            np.save(path, saved)
            volume_files[name] = str(path.resolve())
            del saved
        volume_metadata = {
            "dat": str(args.dat.resolve()),
            "par": str(args.par.resolve()),
            "frame": args.frame,
            "time_code": time_code,
            "gamma": gamma,
            "wave_A": args.wave,
            "shape": list(args.resolution),
            "floating_dtype": "float32",
            "coordinate_system": "uniform Cartesian cell centres in Rsun",
            "cartesian_domain_rsun": [cart_left.tolist(), cart_right.tolist()],
            "coordinate_source": "simesh true stretched AMR coordinates",
            "valid_voxels": int(np.count_nonzero(valid)),
            "total_voxels": int(valid.size),
            "files": volume_files,
        }
        volume_metadata_path = args.volume_output / "volume_metadata.json"
        volume_metadata_path.write_text(json.dumps(volume_metadata, indent=2) + "\n")
    absorbing_volume = dda.VolumeContainer(
        kappa, cloud_emissivity, cart_left * RSUN_CM, cart_right * RSUN_CM
    )
    thin_volume = dda.VolumeContainer(
        np.zeros_like(kappa), emissivity, cart_left * RSUN_CM, cart_right * RSUN_CM
    )
    bottom_spherical = np.array(
        [[spherical_left[0], *(0.5 * (spherical_left[1:] + spherical_right[1:]))]]
    )
    bottom_center = sphere_to_cartesian(bottom_spherical)[0]
    axis_path = np.loadtxt(args.axis_path)
    rope = axis_path[int(np.nanargmax(axis_path[:, 2])), :3]
    rope_spherical = np.array(
        [[rope[2], np.radians(90.0 - rope[1]), np.radians(rope[0])]]
    )
    field_width = 1.10 * max(cart_right - cart_left) * RSUN_CM

    thin_images: dict[str, np.ndarray] = {}
    absorbed_images: dict[str, np.ndarray] = {}
    tau_images: dict[str, np.ndarray] = {}
    view_meta: dict[str, dict[str, object]] = {}
    for view in VIEW_NAMES:
        centre, direction, description = view_geometry(
            view,
            bottom_center,
            rope_spherical,
            float(spherical_right[0]),
            args.outer_margin_rsun,
        )
        plane_points = generate_plane_points(
            centre * RSUN_CM, direction, field_width, args.image_size
        )
        thin, _ = integrate_with_tau(dda, plane_points, direction, thin_volume)
        absorbed, tau = integrate_with_tau(dda, plane_points, direction, absorbing_volume)
        thin = thin.reshape(args.image_size, args.image_size)
        absorbed = absorbed.reshape(args.image_size, args.image_size)
        tau = tau.reshape(args.image_size, args.image_size)
        disk_boundary_intensity = 0.0
        if view == "bottom-inward":
            positive_thin = thin[np.isfinite(thin) & (thin > 0.0)]
            disk_boundary_intensity = float(np.percentile(positive_thin, 75.0))
            disk_mask = np.isfinite(thin) & (thin > 0.0)
            absorbed += (
                disk_boundary_intensity
                * np.exp(-np.minimum(tau, 745.0))
                * disk_mask
            )
        thin_images[view] = thin
        absorbed_images[view] = absorbed
        tau_images[view] = tau
        np.save(args.output / f"aia{args.wave}_{args.frame}_{view}_thin_linear.npy", thin)
        np.save(
            args.output / f"aia{args.wave}_{args.frame}_{view}_absorbed_linear.npy",
            absorbed,
        )
        np.save(args.output / f"aia{args.wave}_{args.frame}_{view}_tau.npy", tau)
        view_meta[view] = {
            "description": description,
            "direction_cartesian": direction.tolist(),
            "centre_cartesian_rsun": centre.tolist(),
            "thin_stats": image_stats(thin),
            "absorbed_stats": image_stats(absorbed),
            "transfer_stats": ratio_stats(absorbed, thin, tau),
            "imposed_disk_boundary_intensity": disk_boundary_intensity,
        }

    all_positive_logs = np.concatenate(
        [
            np.log10(image[np.isfinite(image) & (image > 0.0)])
            for image in thin_images.values()
        ]
    )
    vmin, vmax = np.percentile(all_positive_logs, (1.0, 99.7))
    fig, axes_plot = plt.subplots(2, 2, figsize=(11.5, 10.0), constrained_layout=True)
    shown = None
    for column, view in enumerate(VIEW_NAMES):
        for row, (model, image) in enumerate(
            (
                ("Optically thin", thin_images[view]),
                ("Dark-cloud transfer", absorbed_images[view]),
            )
        ):
            ax = axes_plot[row, column]
            shown = ax.imshow(
                finite_log_image(image, vmin).T,
                origin="lower",
                cmap="afmhot",
                vmin=vmin,
                vmax=vmax,
            )
            ax.set_title(f"{view_meta[view]['description']} | {model}")
            ax.set_xlabel("Image-plane pixel")
            ax.set_ylabel("Image-plane pixel")
            ax.set_aspect("equal")
    fig.suptitle(
        f"Synthetic AIA {args.wave} A continuum-absorption test | "
        f"frame {args.frame} | t={time_code:.3f}"
    )
    fig.colorbar(
        shown,
        ax=axes_plot.ravel().tolist(),
        location="right",
        shrink=0.86,
        label="log10 intensity (arbitrary units)",
    )
    combined_png = args.output / f"aia{args.wave}_{args.frame}_thin_vs_absorbed.png"
    fig.savefig(combined_png, dpi=220)
    plt.close(fig)

    absorbed_positive_logs = np.concatenate(
        [
            np.log10(image[np.isfinite(image) & (image > 0.0)])
            for image in absorbed_images.values()
        ]
    )
    absorbed_vmin, absorbed_vmax = np.percentile(absorbed_positive_logs, (1.0, 99.7))
    fig, absorbed_axes = plt.subplots(1, 2, figsize=(11.5, 5.2), constrained_layout=True)
    absorbed_shown = None
    for ax, view in zip(absorbed_axes, VIEW_NAMES):
        absorbed_shown = ax.imshow(
            finite_log_image(absorbed_images[view], absorbed_vmin).T,
            origin="lower",
            cmap="afmhot",
            vmin=absorbed_vmin,
            vmax=absorbed_vmax,
        )
        ax.set_title(view_meta[view]["description"])
        ax.set_xlabel("Image-plane pixel")
        ax.set_ylabel("Image-plane pixel")
        ax.set_aspect("equal")
    fig.suptitle(
        f"Synthetic absorbing AIA {args.wave} A | frame {args.frame} | t={time_code:.3f}"
    )
    fig.colorbar(
        absorbed_shown,
        ax=absorbed_axes,
        location="right",
        shrink=0.86,
        label="log10 intensity (arbitrary units)",
    )
    absorbed_png = args.output / f"aia{args.wave}_{args.frame}_absorbed_dualview.png"
    fig.savefig(absorbed_png, dpi=220)
    plt.close(fig)

    fig, tau_axes = plt.subplots(1, 2, figsize=(11.5, 5.2), constrained_layout=True)
    tau_shown = None
    positive_tau = np.concatenate(
        [tau[np.isfinite(tau) & (tau > 0.0)] for tau in tau_images.values()]
    )
    tau_max = float(np.percentile(positive_tau, 99.7)) if positive_tau.size else 1.0
    for ax, view in zip(tau_axes, VIEW_NAMES):
        tau_shown = ax.imshow(
            np.log10(1.0 + tau_images[view]).T,
            origin="lower",
            cmap="viridis",
            vmin=0.0,
            vmax=np.log10(1.0 + tau_max),
        )
        ax.set_title(view_meta[view]["description"])
        ax.set_xlabel("Image-plane pixel")
        ax.set_ylabel("Image-plane pixel")
        ax.set_aspect("equal")
    fig.suptitle(f"AIA {args.wave} A optical depth | frame {args.frame} | t={time_code:.3f}")
    fig.colorbar(
        tau_shown,
        ax=tau_axes,
        location="right",
        shrink=0.86,
        label="log10(1 + tau)",
    )
    tau_png = args.output / f"aia{args.wave}_{args.frame}_optical_depth.png"
    fig.savefig(tau_png, dpi=220)
    plt.close(fig)

    metadata = {
        "dat": str(args.dat.resolve()),
        "par": str(args.par.resolve()),
        "axis_path": str(args.axis_path.resolve()),
        "frame": args.frame,
        "time_code": time_code,
        "wave_A": args.wave,
        "radiation_model": (
            "AIA response emissivity plus a pure-absorbing cool cloud with approximate "
            "H I and He I bound-free opacity; exact constant-cell formal transfer solution"
        ),
        "model_scope": (
            "Dark-cloud limiting experiment only; cool-cloud source function is set to "
            "zero, no NLTE He II 304 resonance-line transfer"
        ),
        "opacity_assumption": (
            "MHD rho interpreted as total hydrogen nuclei density using 1e9 cm^-3 per "
            "code unit; H I and He I neutral fraction is 1 below 2e4 K, tapers linearly "
            "to 0 at 1e5 K; He/H=0.1; He II continuum opacity is zero because 304 A "
            "is longward of the 228 A He II ionization edge; opacity below r=1.005 Rsun "
            "is excluded as the chromospheric boundary"
        ),
        "disk_boundary_assumption": (
            "Radial top-down view adds a uniform lower-boundary intensity equal to the "
            "75th percentile of its optically thin positive pixels; tangent view has no "
            "disk boundary"
        ),
        "uniform_cartesian_shape": list(args.resolution),
        "image_shape": [args.image_size, args.image_size],
        "top_view_detector_radius_rsun": float(spherical_right[0] + args.outer_margin_rsun),
        "coordinate_source": "simesh true stretched AMR coordinates",
        "ray_traversal": (
            "axis-aligned slab entry plus Cartesian DDA; all distances in cm"
        ),
        "shared_log10_limits": [float(vmin), float(vmax)],
        "absorbed_log10_limits": [float(absorbed_vmin), float(absorbed_vmax)],
        "bottom_boundary_center_spherical": bottom_spherical[0].tolist(),
        "rope_apex_spherical": rope_spherical[0].tolist(),
        "cartesian_domain_rsun": [cart_left.tolist(), cart_right.tolist()],
        "reusable_volume_output": (
            str(args.volume_output.resolve()) if args.volume_output is not None else None
        ),
        "reusable_volume_files": volume_files,
        "views": view_meta,
        "output_png": str(combined_png.resolve()),
        "absorbed_dualview_png": str(absorbed_png.resolve()),
        "optical_depth_png": str(tau_png.resolve()),
    }
    metadata_path = args.output / f"aia{args.wave}_{args.frame}_absorption.json"
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    print(combined_png)


if __name__ == "__main__":
    main()

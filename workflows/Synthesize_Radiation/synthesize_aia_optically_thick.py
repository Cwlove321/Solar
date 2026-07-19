#!/usr/bin/env python3
"""Reusable AIA optically-thick synthesis from AMRVAC .dat files.

The script uses the current simesh API (`simesh.amrvac.read_uniform`) and the
local RST_cart radiation helpers. It does not use simesh.legacy.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import sunpy.map  # noqa: F401; registers sdoaia* colormaps


MASS_H = 1.674e-24
K_B = 1.3806e-16
AXIS_EPS = 1e-12
REPO_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_DDA_ROOT = REPO_ROOT / "Packages" / "radsyn_tools" / "DDA_method1"
RADSYN_ROOT = DEFAULT_DDA_ROOT
DDA = None
physical_module = None


def json_default(obj):
    if isinstance(obj, np.integer):
        return int(obj)
    if isinstance(obj, np.floating):
        return float(obj)
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj).__name__} is not JSON serializable")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Synthesize optically-thick AIA images from an AMRVAC .dat file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("dat", help="Path to the AMRVAC .dat file.")
    parser.add_argument("--output-dir", default=None, help="Directory for PNG/NPY/metadata outputs.")
    parser.add_argument(
        "--dda-root",
        default=str(DEFAULT_DDA_ROOT),
        help="Directory containing RST_cart/ and AIA_response/.",
    )
    parser.add_argument("--wavelength", type=int, default=171, help="AIA wavelength, e.g. 171, 211, 304.")
    parser.add_argument("--normal", type=float, nargs=3, default=[0.0, 1.0, 0.0], help="Line-of-sight direction.")
    parser.add_argument("--resolution", type=int, nargs=3, default=None, help="Uniform grid resolution nx ny nz.")
    parser.add_argument(
        "--native-resolution",
        action="store_true",
        help="Use domain_nx * 2**(levmax-1). If omitted, this is still the default when --resolution is not set.",
    )
    parser.add_argument("--z-min", type=float, default=None, help="Discard data below this z height in Mm before synthesis.")
    parser.add_argument("--i0", type=float, default=0.0, help="Background intensity I0.")
    parser.add_argument("--kappa-scale", type=float, default=1.0, help="Multiply absorption coefficient by this factor.")
    parser.add_argument("--emission-scale", type=float, default=1.0, help="Multiply emissivity by this factor.")
    parser.add_argument("--clip-negative-kappa", action=argparse.BooleanOptionalAction, default=True)
    parser.add_argument("--vmin", type=float, default=None, help="Plot vmin in log10 intensity.")
    parser.add_argument("--vmax", type=float, default=None, help="Plot vmax in log10 intensity.")
    parser.add_argument("--percentile-vmin", type=float, default=5.0, help="Auto vmin percentile of positive log data.")
    parser.add_argument("--percentile-vmax", type=float, default=99.0, help="Auto vmax percentile of positive log data.")
    parser.add_argument("--label", default=None, help="Extra label used in output filenames.")
    parser.add_argument("--dpi", type=int, default=220, help="PNG DPI.")
    return parser.parse_args()


def configure_radsyn_tools(dda_root: str | Path) -> Path:
    global RADSYN_ROOT, DDA, physical_module
    RADSYN_ROOT = Path(dda_root).expanduser().resolve()
    if not (RADSYN_ROOT / "RST_cart").exists():
        raise FileNotFoundError(f"RST_cart not found under --dda-root: {RADSYN_ROOT}")
    if str(RADSYN_ROOT) not in sys.path:
        sys.path.insert(0, str(RADSYN_ROOT))
    import RST_cart.DDA as dda_module
    import RST_cart.physical_module as physical_module_module

    DDA = dda_module
    physical_module = physical_module_module
    return RADSYN_ROOT


def load_uniform_physical(dat_file: Path, resolution: tuple[int, int, int], box_min: np.ndarray, box_max: np.ndarray):
    from simesh.amrvac import read_uniform

    udata = read_uniform(
        str(dat_file),
        resolution=resolution,
        bounds=(box_min, box_max),
        field_indices=[0, 4],
        interpolation="zero",
    )
    rho = udata[..., 0]
    energy = udata[..., 1]
    numberdensity = rho / (1.4 * MASS_H) * 2.341e-15
    pressure = 2.0 / 3.0 * energy * 0.317538
    temperature = pressure / numberdensity / (2.3 * K_B)
    return numberdensity, pressure, temperature


def calculate_euv_terms(numberdensity: np.ndarray, pressure: np.ndarray, temperature: np.ndarray, wavelength: int):
    if physical_module is None:
        raise RuntimeError("RST_cart.physical_module is not configured")
    kappa, emissivity = physical_module.cal_euv(
        np.nan_to_num(numberdensity, nan=0.0),
        np.nan_to_num(pressure, nan=0.0),
        np.nan_to_num(temperature, nan=0.0),
        wavelength,
        True,
    )
    return kappa, emissivity


def normalized_normal(normal: np.ndarray) -> np.ndarray:
    norm = np.linalg.norm(normal)
    if norm < AXIS_EPS:
        raise ValueError("normal vector length is zero")
    return normal / norm


def axis_index_and_sign(normal_unit: np.ndarray) -> tuple[int, int] | None:
    for axis in range(3):
        basis = np.zeros(3)
        basis[axis] = 1.0
        if np.linalg.norm(normal_unit - basis) < AXIS_EPS:
            return axis, 1
        if np.linalg.norm(normal_unit + basis) < AXIS_EPS:
            return axis, -1
    return None


def integrate_axis(
    kappa: np.ndarray,
    emissivity: np.ndarray,
    normal_unit: np.ndarray,
    box_min: np.ndarray,
    box_max: np.ndarray,
    resolution: tuple[int, int, int],
    i0: float,
    emission_scale: float,
):
    axis_info = axis_index_and_sign(normal_unit)
    if axis_info is None:
        raise ValueError("normal is not axis aligned")
    axis, sign = axis_info
    ds_cm = (box_max[axis] - box_min[axis]) * 1e9 / resolution[axis]
    dtau = kappa * ds_cm
    j = emissivity * emission_scale
    if sign < 0:
        dtau_work = np.flip(dtau, axis=axis)
        j_work = np.flip(j, axis=axis)
    else:
        dtau_work = dtau
        j_work = j
    tau_before = np.cumsum(dtau_work, axis=axis) - dtau_work
    image = np.sum(np.exp(-tau_before) * j_work * ds_cm, axis=axis)
    tau = np.sum(dtau, axis=axis)
    background = i0 * np.exp(-tau)
    image = image + background
    return orient_axis_image(image, tau, axis, box_min, box_max), ds_cm


def orient_axis_image(image: np.ndarray, tau: np.ndarray, axis: int, box_min: np.ndarray, box_max: np.ndarray):
    if axis == 0:
        # image axes are y,z -> plot as z,y
        return image.T, tau.T, [float(box_min[1]), float(box_max[1]), float(box_min[2]), float(box_max[2])], "Y [Mm]", "Z [Mm]"
    if axis == 1:
        # image axes are x,z -> plot as z,x
        return image.T, tau.T, [float(box_min[0]), float(box_max[0]), float(box_min[2]), float(box_max[2])], "X [Mm]", "Z [Mm]"
    # axis == 2: image axes are x,y -> plot as y,x
    return image.T, tau.T, [float(box_min[0]), float(box_max[0]), float(box_min[1]), float(box_max[1])], "X [Mm]", "Y [Mm]"


def detector_basis(normal: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    n = normalized_normal(normal)
    vertical = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(n, vertical)) > 0.98:
        vertical = np.array([1.0, 0.0, 0.0])
    u = np.cross(n, vertical)
    u /= np.linalg.norm(u)
    v = np.cross(u, n)
    v /= np.linalg.norm(v)
    return n, u, v


def projected_box_detector(
    left_edge: np.ndarray,
    right_edge: np.ndarray,
    center: np.ndarray,
    normal: np.ndarray,
    pixel_cm: float,
) -> tuple[np.ndarray, list[float], tuple[int, int], tuple[np.ndarray, np.ndarray, np.ndarray]]:
    n, u, v = detector_basis(normal)
    corners = np.array(
        [
            [x, y, z]
            for x in (left_edge[0], right_edge[0])
            for y in (left_edge[1], right_edge[1])
            for z in (left_edge[2], right_edge[2])
        ],
        dtype=float,
    )
    rel = corners - center
    u_proj = rel @ u
    v_proj = rel @ v
    u_min, u_max = float(u_proj.min()), float(u_proj.max())
    v_min, v_max = float(v_proj.min()), float(v_proj.max())
    n_u = int(np.ceil((u_max - u_min) / pixel_cm))
    n_v = int(np.ceil((v_max - v_min) / pixel_cm))
    u_coords = u_min + (np.arange(n_u) + 0.5) * (u_max - u_min) / n_u
    v_coords = v_min + (np.arange(n_v) + 0.5) * (v_max - v_min) / n_v
    uu, vv = np.meshgrid(u_coords, v_coords)
    points = center + np.outer(uu.ravel(), u) + np.outer(vv.ravel(), v)
    extent_m_m = [u_min / 1e9, u_max / 1e9, v_min / 1e9, v_max / 1e9]
    return points, extent_m_m, (n_v, n_u), (n, u, v)


def integrate_dda(
    kappa: np.ndarray,
    emissivity: np.ndarray,
    normal: np.ndarray,
    box_min: np.ndarray,
    box_max: np.ndarray,
    resolution: tuple[int, int, int],
    i0: float,
    emission_scale: float,
):
    if DDA is None:
        raise RuntimeError("RST_cart.DDA is not configured")
    left_edge_cm = box_min * 1e9
    right_edge_cm = box_max * 1e9
    center_cm = (left_edge_cm + right_edge_cm) / 2.0
    cell_sizes_cm = (box_max - box_min) * 1e9 / np.array(resolution)
    pixel_cm = float(max(cell_sizes_cm[0], cell_sizes_cm[1]))
    points, extent, shape, basis = projected_box_detector(left_edge_cm, right_edge_cm, center_cm, normal, pixel_cm)
    vc = DDA.VolumeContainer(kappa, emissivity * emission_scale, left_edge_cm, right_edge_cm)
    n = normalized_normal(normal)
    max_t = float(np.linalg.norm(right_edge_cm - left_edge_cm) * 2.0)
    intensity = np.zeros(points.shape[0], dtype=np.float64)
    tau = np.zeros(points.shape[0], dtype=np.float64)
    for i in range(points.shape[0]):
        rd = DDA.RayData()
        DDA.walk_volume(vc, points[i].astype(np.float64), n, rd, max_t=max_t)
        intensity[i] = rd.sum_val
        tau[i] = rd.tau
    image = intensity.reshape(shape) + i0 * np.exp(-tau.reshape(shape))
    tau_image = tau.reshape(shape)
    return image, tau_image, extent, "Detector U [Mm]", "Detector V [Mm]", basis


def log_image(image: np.ndarray) -> np.ndarray:
    image_log = np.zeros_like(image)
    mask = image > 0
    image_log[mask] = np.log10(image[mask])
    return np.nan_to_num(image_log, nan=0.0, neginf=0.0, posinf=0.0)


def auto_limits(image: np.ndarray, image_log: np.ndarray, args: argparse.Namespace) -> tuple[float, float]:
    positive = image_log[np.isfinite(image_log) & (image > 0)]
    finite = image_log[np.isfinite(image_log)]
    if args.vmin is not None:
        vmin = args.vmin
    elif positive.size:
        vmin = float(np.nanpercentile(positive, args.percentile_vmin))
    else:
        vmin = float(np.nanmin(finite)) if finite.size else 0.0
    if args.vmax is not None:
        vmax = args.vmax
    elif positive.size:
        vmax = float(np.nanpercentile(positive, args.percentile_vmax))
    else:
        vmax = float(np.nanmax(finite)) if finite.size else 1.0
    return vmin, vmax


def label_from_args(dat_file: Path, args: argparse.Namespace) -> str:
    if args.label:
        return args.label
    normal_label = "normal" + "_".join(f"{x:g}".replace("-", "m").replace(".", "p") for x in args.normal)
    return f"aia{args.wavelength}_{normal_label}_{dat_file.stem}"


def main() -> None:
    args = parse_args()
    from simesh.amrvac.api import open_dataset

    dda_root = configure_radsyn_tools(args.dda_root)
    dat_file = Path(args.dat).expanduser().resolve()
    output_dir = Path(args.output_dir).expanduser().resolve() if args.output_dir else dat_file.parent / "radiation"
    output_dir.mkdir(parents=True, exist_ok=True)

    ds = open_dataset(str(dat_file))
    ds.load_metadata()
    original_box_min = ds.metadata["xmin"].astype(float)
    box_max = ds.metadata["xmax"].astype(float)
    domain_nx = ds.metadata["domain_nx"].astype(int)
    levmax = int(ds.metadata.get("levmax", 1))
    resolution = tuple(args.resolution) if args.resolution else tuple((domain_nx * (2 ** (levmax - 1))).astype(int))
    box_min = original_box_min.copy()
    if args.z_min is not None:
        box_min[2] = max(box_min[2], args.z_min)

    numberdensity, pressure, temperature = load_uniform_physical(dat_file, resolution, box_min, box_max)
    kappa, emissivity = calculate_euv_terms(numberdensity, pressure, temperature, args.wavelength)
    if args.clip_negative_kappa:
        kappa = np.maximum(kappa, 0.0)
    kappa = kappa * args.kappa_scale

    normal = np.array(args.normal, dtype=float)
    normal_unit = normalized_normal(normal)
    axis_info = axis_index_and_sign(normal_unit)
    if axis_info is not None:
        (image, tau_image, extent, xlabel, ylabel), ds_cm = integrate_axis(
            kappa, emissivity, normal_unit, box_min, box_max, resolution, args.i0, args.emission_scale
        )
        integration = "axis"
        detector_basis_info = None
    else:
        image, tau_image, extent, xlabel, ylabel, detector_basis_info = integrate_dda(
            kappa, emissivity, normal, box_min, box_max, resolution, args.i0, args.emission_scale
        )
        ds_cm = None
        integration = "dda"

    image_log = log_image(image)
    vmin, vmax = auto_limits(image, image_log, args)
    cmap_name = f"sdoaia{args.wavelength}"
    cmap = cmap_name if cmap_name in plt.colormaps() else "hot"

    label = label_from_args(dat_file, args)
    prefix = f"{label}_thick_I0_{args.i0:g}".replace(".", "p").replace("+", "")
    raw_path = output_dir / f"{prefix}_raw.npy"
    log_path = output_dir / f"{prefix}_log.npy"
    tau_path = output_dir / f"{prefix}_tau.npy"
    png_path = output_dir / f"{prefix}.png"
    metadata_path = output_dir / f"{prefix}_metadata.json"

    np.save(raw_path, image)
    np.save(log_path, image_log)
    np.save(tau_path, tau_image)

    fig, ax = plt.subplots(figsize=(9.0, 7.6))
    im = ax.imshow(
        image_log,
        origin="lower",
        extent=extent,
        cmap=cmap,
        vmin=vmin,
        vmax=vmax,
        aspect="equal",
        interpolation="nearest",
    )
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"AIA {args.wavelength}, normal={args.normal}")
    cbar = fig.colorbar(im, ax=ax, shrink=0.86, pad=0.03)
    cbar.set_label("log10(Intensity)")
    fig.tight_layout()
    fig.savefig(png_path, dpi=args.dpi)
    plt.close(fig)

    metadata = {
        "source_dat": str(dat_file),
        "dda_root": str(dda_root),
        "simesh_api": "simesh.amrvac.read_uniform",
        "uses_simesh_legacy": False,
        "uses_cached_npz": False,
        "wavelength": args.wavelength,
        "absorption": True,
        "I0": args.i0,
        "normal": args.normal,
        "normal_unit": normal_unit.tolist(),
        "integration": integration,
        "original_box_min_Mm": original_box_min.tolist(),
        "box_min_Mm": box_min.tolist(),
        "box_max_Mm": box_max.tolist(),
        "domain_nx": domain_nx.tolist(),
        "levmax": levmax,
        "uniform_resolution_xyz": list(resolution),
        "z_min_cut_Mm": args.z_min,
        "kappa_scale": args.kappa_scale,
        "emission_scale": args.emission_scale,
        "clip_negative_kappa": args.clip_negative_kappa,
        "plot_extent": extent,
        "plot_labels": [xlabel, ylabel],
        "output_shape": list(image.shape),
        "axis_ds_cm": ds_cm,
        "detector_basis": None
        if detector_basis_info is None
        else {
            "normal_unit": detector_basis_info[0].tolist(),
            "u_unit": detector_basis_info[1].tolist(),
            "v_unit": detector_basis_info[2].tolist(),
        },
        "colormap": cmap,
        "plot_vmin_vmax": [vmin, vmax],
        "raw_min_max": [float(np.nanmin(image)), float(np.nanmax(image))],
        "log_min_max": [float(np.nanmin(image_log)), float(np.nanmax(image_log))],
        "tau_min_max": [float(np.nanmin(tau_image)), float(np.nanmax(tau_image))],
        "numberdensity_min_max": [float(np.nanmin(numberdensity)), float(np.nanmax(numberdensity))],
        "temperature_min_max": [float(np.nanmin(temperature)), float(np.nanmax(temperature))],
        "output_raw_npy": str(raw_path),
        "output_log_npy": str(log_path),
        "output_tau_npy": str(tau_path),
        "output_png": str(png_path),
    }
    metadata_path.write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2, default=json_default),
        encoding="utf-8",
    )

    print(f"raw: {raw_path}")
    print(f"log: {log_path}")
    print(f"tau: {tau_path}")
    print(f"png: {png_path}")
    print(f"metadata: {metadata_path}")


if __name__ == "__main__":
    main()

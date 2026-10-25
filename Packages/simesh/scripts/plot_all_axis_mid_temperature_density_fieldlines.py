#!/usr/bin/env python3
"""Plot one fixed-phi temperature/density frame for every available DAT."""

from __future__ import annotations

import gc
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle
from simesh.amrvac import open_dataset


ROOT = Path("/share/home/guoyang1/zhouchang/20241009/rbsl14")
RUN = ROOT / (
    "fullmhd_stretched_100km_block12_cubicfix_dynamicamr_"
    "relaxed0p30_rbsl_flux1p8_rho20_eta200_tc_trac_norad"
)
INSERT = ROOT / (
    "rbsl_insert_stretched_100km_block12_cubicfix_box_l4_200Mm_"
    "relaxed0p30_flux1p8_rho20_ctfix"
)
PREFIX = "r14_f18_"

RSUN_MM = 695.5
UNIT_NUMBER_DENSITY = 1.0e9
UNIT_TEMPERATURE = 1.0e6
UNIT_TIME_S = 5972.477161555993
HEIGHT_RANGE_MM = (0.0, 120.0)
LATITUDE_RANGE_DEG = (8.0, 30.0)
TEMPERATURE_NORM = LogNorm(8.0e3, 4.0e6)
DENSITY_NORM = LogNorm(2.0e8, 1.5e14)


def axis_midpoint() -> tuple[float, float, float]:
    path = np.loadtxt(INSERT / "path.dat")
    index = int(np.nanargmax(path[:, 2]))
    return tuple(float(value) for value in path[index, :3])


def extract_fixed_phi(
    path: Path, phi_deg: float, height: np.ndarray, latitude: np.ndarray
) -> dict[str, np.ndarray | float]:
    """Sample a fixed-phi plane using simesh's true AMRVAC coordinates.

    ``sample_points`` receives spherical (r, theta, phi) positions, locates the
    owning AMR leaf with its physical node bounds, and interpolates its values.
    Passing ``amrvac.par`` is deliberate: a DAT header has no complete stretch
    definition, so this prevents a stretched radial mesh being treated as uniform.
    """
    dataset = open_dataset(str(path), grid_config=str(RUN / "amrvac.par"))
    dataset.load_data()
    names = list(dataset.metadata["w_names"])
    param_names = list(dataset.metadata["param_names"])
    gamma = float(dataset.metadata["params"][param_names.index("gamma")])
    wanted = ("rho", "e", "b1", "b2")
    hh, ll = np.meshgrid(height, latitude, indexing="ij")
    points = np.column_stack((
        1.0 + hh.ravel() / RSUN_MM,
        np.radians(90.0 - ll.ravel()),
        np.full(hh.size, np.radians(phi_deg)),
    ))
    indices = [names.index(name) for name in wanted]
    samples = dataset.sample_points(points, field_indices=indices, interpolation="linear")
    fields = {
        name: samples[:, index].reshape(hh.shape)
        for index, name in enumerate(wanted)
    }
    fields["height"] = height
    fields["latitude"] = latitude
    fields["gamma"] = gamma
    fields["time"] = float(dataset.metadata["time"])
    fields["coordinate_mode"] = "simesh true stretched AMR coordinates"
    # Leaf blocks that physically intersect this phi plane.  These bounds are
    # from the reconstructed AMRVAC coordinate object, not logical AMR indices.
    node_r, node_theta, node_phi = dataset.block_node_coordinates
    phi = np.radians(phi_deg)
    epsilon = 32.0 * np.finfo(float).eps
    intersects_phi = (node_phi[:, 0] - epsilon <= phi) & (phi <= node_phi[:, -1] + epsilon)
    fields["amr_bounds"] = np.column_stack((
        node_r[intersects_phi, 0], node_r[intersects_phi, -1],
        node_theta[intersects_phi, 0], node_theta[intersects_phi, -1],
        np.asarray(dataset.tree_info[0], dtype=int)[intersects_phi],
    ))
    del dataset
    gc.collect()
    return fields


def plot_frame(
    path: Path,
    temperature_ax,
    density_ax,
    phi_deg: float,
    axis_latitude_deg: float,
    axis_height_mm: float,
) -> dict[str, float | str]:
    height = np.linspace(*HEIGHT_RANGE_MM, 421)
    latitude = np.linspace(*LATITUDE_RANGE_DEG, 521)
    data = extract_fixed_phi(path, phi_deg, height, latitude)
    rho = np.asarray(data["rho"])
    eint = np.asarray(data["e"])
    br = np.asarray(data["b1"])
    btheta = np.asarray(data["b2"])
    gamma = float(data["gamma"])
    time_code = float(data["time"])
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        temperature = (gamma - 1.0) * eint / rho * UNIT_TEMPERATURE
        density = rho * UNIT_NUMBER_DENSITY
    temperature[~np.isfinite(temperature) | (temperature <= 0.0)] = np.nan
    density[~np.isfinite(density) | (density <= 0.0)] = np.nan

    distance = np.radians(latitude - axis_latitude_deg) * RSUN_MM
    temperature_ax.pcolormesh(
        distance, height, temperature, shading="auto", cmap="inferno",
        norm=TEMPERATURE_NORM, rasterized=True,
    )
    density_ax.pcolormesh(
        distance, height, density, shading="auto", cmap="viridis",
        norm=DENSITY_NORM, rasterized=True,
    )

    rr = 1.0 + height[:, None] / RSUN_MM
    horizontal_b = -btheta / rr
    magnitude = np.hypot(horizontal_b, br)
    horizontal_b = np.nan_to_num(horizontal_b / np.maximum(magnitude, 1.0e-30))
    vertical_b = np.nan_to_num(br / np.maximum(magnitude, 1.0e-30))
    step = (slice(None, None, 4), slice(None, None, 4))
    for ax in (temperature_ax, density_ax):
        ax.streamplot(
            distance[step[1]], height[step[0]], horizontal_b[step], vertical_b[step],
            color="white", density=0.85, linewidth=0.45, arrowsize=0.50,
            minlength=0.15, maxlength=4.0,
        )
        ax.scatter(
            [0.0], [axis_height_mm], marker="*", s=46, c="#37d5ff",
            edgecolors="black", linewidths=0.4, zorder=8,
        )
        ax.axvline(0.0, color="#37d5ff", linewidth=0.6, alpha=0.55)
        ax.set_aspect("equal", adjustable="box")
        ax.set_xlim(float(distance[0]), float(distance[-1]))
        ax.set_ylim(*HEIGHT_RANGE_MM)
        ax.grid(color="white", alpha=0.08, linewidth=0.35)

    title = f"{path.stem[-4:]}  t={time_code:.3f} ({time_code * UNIT_TIME_S:.1f} s)"
    temperature_ax.set_title(title, fontsize=10)
    finite_temperature = temperature[np.isfinite(temperature)]
    finite_density = density[np.isfinite(density)]
    del data, rho, eint, br, btheta, temperature, density
    gc.collect()
    return {
        "file": path.name,
        "time_code": time_code,
        "time_s": time_code * UNIT_TIME_S,
        "temperature_min_K": float(np.min(finite_temperature)),
        "temperature_max_K": float(np.max(finite_temperature)),
        "density_min_cm-3": float(np.min(finite_density)),
        "density_max_cm-3": float(np.max(finite_density)),
    }


def render_frame(
    path: Path,
    output: Path,
    phi_deg: float,
    axis_latitude_deg: float,
    axis_height_mm: float,
) -> dict[str, float | str]:
    fig, axes = plt.subplots(
        1, 2, figsize=(12.0, 6.3), constrained_layout=True,
        sharex=True, sharey=True,
    )
    metrics = plot_frame(
        path,
        axes[0],
        axes[1],
        phi_deg,
        axis_latitude_deg,
        axis_height_mm,
    )

    axes[0].set_title("Temperature", fontsize=11)
    axes[1].set_title("Number density", fontsize=11)
    for ax in axes:
        ax.set_xlabel("Surface distance from axis latitude [Mm]")
    axes[0].set_ylabel("Height [Mm]")
    fig.colorbar(
        ScalarMappable(norm=TEMPERATURE_NORM, cmap="inferno"),
        ax=axes[0],
        location="right",
        shrink=0.94,
        pad=0.02,
        label="Temperature [K]",
    )
    fig.colorbar(
        ScalarMappable(norm=DENSITY_NORM, cmap="viridis"),
        ax=axes[1],
        location="right",
        shrink=0.94,
        pad=0.02,
        label=r"Number density [cm$^{-3}$]",
    )
    fig.suptitle(
        rf"{path.stem[-4:]}  $t={metrics['time_code']:.3f}$ "
        rf"({metrics['time_s']:.1f} s); fixed $\phi={phi_deg:.3f}^\circ$; "
        rf"axis: latitude $={axis_latitude_deg:.3f}^\circ$, "
        rf"height $={axis_height_mm:.1f}$ Mm",
        fontsize=12,
    )
    metrics_output = output.with_suffix(".json")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)
    metrics_output.write_text(
        json.dumps(
            {
                "phi_deg": phi_deg,
                "axis_latitude_deg": axis_latitude_deg,
                "axis_height_Mm": axis_height_mm,
                "height_range_Mm": HEIGHT_RANGE_MM,
                "latitude_range_deg": LATITUDE_RANGE_DEG,
                "coordinate_mode": "simesh true stretched AMR coordinates",
                "frame": metrics,
            },
            indent=2,
        )
        + "\n"
    )
    return metrics


def render_amr_grid_frame(
    path: Path, output: Path, phi_deg: float, axis_latitude_deg: float,
) -> None:
    """Render temperature with physical AMR leaf-block boundaries overlaid."""
    height = np.linspace(*HEIGHT_RANGE_MM, 421)
    latitude = np.linspace(*LATITUDE_RANGE_DEG, 521)
    data = extract_fixed_phi(path, phi_deg, height, latitude)
    rho = np.asarray(data["rho"])
    eint = np.asarray(data["e"])
    gamma = float(data["gamma"])
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        temperature = (gamma - 1.0) * eint / rho * UNIT_TEMPERATURE
    temperature[~np.isfinite(temperature) | (temperature <= 0.0)] = np.nan

    distance = np.radians(latitude - axis_latitude_deg) * RSUN_MM
    fig, ax = plt.subplots(figsize=(12.0, 6.3), constrained_layout=True)
    image = ax.pcolormesh(
        distance, height, temperature, shading="auto", cmap="inferno",
        norm=TEMPERATURE_NORM, rasterized=True,
    )
    bounds = np.asarray(data["amr_bounds"], dtype=float)
    levels = bounds[:, 4].astype(int)
    level_min, level_max = int(levels.min()), int(levels.max())
    colors = plt.get_cmap("cool")(np.linspace(0.25, 0.95, level_max - level_min + 1))
    for r0, r1, theta0, theta1, level in bounds:
        latitude0, latitude1 = 90.0 - np.degrees((theta0, theta1))
        x0 = np.radians(min(latitude0, latitude1) - axis_latitude_deg) * RSUN_MM
        x1 = np.radians(max(latitude0, latitude1) - axis_latitude_deg) * RSUN_MM
        y0, y1 = (r0 - 1.0) * RSUN_MM, (r1 - 1.0) * RSUN_MM
        ax.add_patch(Rectangle(
            (x0, y0), x1 - x0, y1 - y0, fill=False,
            edgecolor=colors[int(level) - level_min], linewidth=0.35, alpha=0.85,
        ))
    ax.set_xlim(float(distance[0]), float(distance[-1]))
    ax.set_ylim(*HEIGHT_RANGE_MM)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Surface distance from axis latitude [Mm]")
    ax.set_ylabel("Height [Mm]")
    ax.set_title(
        f"{path.stem[-4:]}: temperature with AMR leaf-block boundaries "
        f"at $\\phi={phi_deg:.3f}^\\circ$ (levels {level_min}–{level_max})"
    )
    fig.colorbar(image, ax=ax, shrink=0.94, pad=0.02, label="Temperature [K]")
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    phi_deg, axis_latitude_deg, axis_radius = axis_midpoint()
    axis_height_mm = (axis_radius - 1.0) * RSUN_MM
    frames_dir = RUN / "python" / "fixed_phi_axis_mid_temperature_density_fieldlines_frames"
    frames_dir.mkdir(parents=True, exist_ok=True)
    processed: set[str] = set()
    all_metrics: dict[str, dict[str, float | str]] = {}

    while True:
        files = sorted((RUN / "data").glob(f"{PREFIX}[0-9][0-9][0-9][0-9].dat"))
        if not files:
            raise FileNotFoundError(f"No DAT files found below {RUN / 'data'}")
        pending = [path for path in files if path.name not in processed]
        if not pending:
            break

        for path in pending:
            frame_id = path.stem[-4:]
            output = frames_dir / f"fixed_phi_axis_mid_temperature_density_fieldlines_{frame_id}.png"
            metrics_output = output.with_suffix(".json")
            if (
                output.is_file()
                and metrics_output.is_file()
                and json.loads(metrics_output.read_text()).get("coordinate_mode")
                == "simesh true stretched AMR coordinates"
            ):
                saved = json.loads(metrics_output.read_text())
                metrics = saved["frame"]
                print(f"SKIP {frame_id}: {output}", flush=True)
            else:
                metrics = render_frame(
                    path,
                    output,
                    phi_deg,
                    axis_latitude_deg,
                    axis_height_mm,
                )
                print(f"DONE {frame_id}: {output}", flush=True)
            grid_output = frames_dir / f"fixed_phi_axis_mid_temperature_amr_grid_{frame_id}.png"
            if not grid_output.is_file():
                render_amr_grid_frame(path, grid_output, phi_deg, axis_latitude_deg)
                print(f"DONE AMR GRID {frame_id}: {grid_output}", flush=True)
            processed.add(path.name)
            all_metrics[frame_id] = metrics

    summary = frames_dir / "fixed_phi_axis_mid_temperature_density_fieldlines_metrics.json"
    summary.write_text(
        json.dumps(
            {
                "phi_deg": phi_deg,
                "axis_latitude_deg": axis_latitude_deg,
                "axis_height_Mm": axis_height_mm,
                "height_range_Mm": HEIGHT_RANGE_MM,
                "latitude_range_deg": LATITUDE_RANGE_DEG,
                "coordinate_mode": "simesh true stretched AMR coordinates",
                "frames": [all_metrics[key] for key in sorted(all_metrics)],
            },
            indent=2,
        )
        + "\n"
    )
    print(f"SUMMARY {summary}", flush=True)


if __name__ == "__main__":
    main()

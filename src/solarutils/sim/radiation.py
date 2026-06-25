"""Synthetic radiation helpers for Cartesian MHD cubes."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
import os
import sys

import matplotlib.pyplot as plt
import numpy as np


MASS_H_CGS = 1.674e-24
K_B_CGS = 1.3806e-16


@dataclass(frozen=True)
class RadiationPlane:
    """Square detector plane used for DDA line-of-sight integration."""

    left_edge: np.ndarray
    right_edge: np.ndarray
    center: np.ndarray
    normal: np.ndarray
    width: float | tuple[float, float]
    resolution: int | tuple[int, int]

    @classmethod
    def from_notebook_defaults(cls, *, normal=(0.0, 1.0, 0.0), resolution: int = 512):
        """Return the plane geometry used by ``Desktop/jet_mhd/辐射合成.ipynb``."""

        return cls(
            left_edge=np.array([-3.0, -3.0, 0.1]) * 1e9,
            right_edge=np.array([3.0, 3.0, 12.1]) * 1e9,
            center=np.array([0.0, 0.0, 6.1]) * 1e9,
            normal=np.asarray(normal, dtype=float),
            width=12.0e9,
            resolution=resolution,
        )


@dataclass(frozen=True)
class SlabBounds:
    """Uniform slab bounds in simulation coordinates."""

    lower: tuple[float, float, float] = (-1.5, -3.0, 0.1)
    upper: tuple[float, float, float] = (4.5, 3.0, 12.1)


def density_pressure_temperature_from_amrvac(
    rho_code,
    e_code,
    *,
    rho_scale: float = 2.341e-15,
    pressure_scale: float = 0.317538,
    mean_particle_factor: float = 2.3,
    helium_factor: float = 1.4,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Convert AMRVAC ``rho`` and ``e`` arrays to number density, pressure, and temperature.

    This preserves the conversion used in the reference notebook:
    ``n = rho / (1.4 m_H) * 2.341e-15``,
    ``p = 2/3 e * 0.317538``, and ``T = p / n / (2.3 k_B)``.
    """

    rho_code = np.asarray(rho_code, dtype=float)
    e_code = np.asarray(e_code, dtype=float)
    if rho_code.shape != e_code.shape:
        raise ValueError("rho_code and e_code must have the same shape")

    number_density = rho_code / (helium_factor * MASS_H_CGS) * rho_scale
    pressure = (2.0 / 3.0) * e_code * pressure_scale
    temperature = np.divide(
        pressure,
        number_density * mean_particle_factor * K_B_CGS,
        out=np.zeros_like(pressure, dtype=float),
        where=number_density != 0.0,
    )
    return number_density, pressure, temperature


def clean_cube_for_dda(cube, *, transpose: bool = True) -> np.ndarray:
    """Return a finite 3D array in the orientation expected by the DDA routines."""

    arr = np.asarray(cube, dtype=float)
    if arr.ndim != 3:
        raise ValueError("cube must be 3D")
    if transpose:
        arr = arr.T
    return np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)


def synthesize_radiation_image(
    number_density,
    pressure,
    temperature,
    plane: RadiationPlane,
    *,
    wavelength: int | str = 171,
    radius=None,
    absorption: bool = False,
    dda_root: str | Path | None = None,
    transpose: bool = True,
) -> np.ndarray:
    """Integrate a synthetic radiation image through a Cartesian cube.

    Parameters are cgs arrays matching the reference DDA tools. The returned
    image has shape ``(plane.resolution, plane.resolution)``.
    """

    rho = clean_cube_for_dda(number_density, transpose=transpose)
    pres = clean_cube_for_dda(pressure, transpose=transpose)
    temp = clean_cube_for_dda(temperature, transpose=transpose)
    if rho.shape != pres.shape or rho.shape != temp.shape:
        raise ValueError("number_density, pressure, and temperature must have the same shape")

    if radius is None:
        radius = np.zeros(rho.shape, dtype=float)
    else:
        radius = clean_cube_for_dda(radius, transpose=transpose)

    with import_dda_tools(dda_root) as tools:
        physical_module = tools["physical_module"]
        points = tools["points"]
        dda = tools["dda"]

        if wavelength == "radio":
            kappa_data, j_data = physical_module.cal_radio(rho, pres, temp)
        elif wavelength == "white_light":
            kappa_data, j_data = physical_module.cal_white_light(rho, radius)
        else:
            kappa_data, j_data = physical_module.cal_euv(rho, pres, temp, int(wavelength), absorption)

        vc = dda.VolumeContainer(
            np.asarray(kappa_data, dtype=float),
            np.asarray(j_data, dtype=float),
            np.asarray(plane.left_edge, dtype=float),
            np.asarray(plane.right_edge, dtype=float),
        )
        plane_points, image_shape = generate_detector_points(plane, points_module=points)
        image = dda.integrate_plane_points(
            plane_points,
            np.asarray(plane.normal, dtype=float),
            vc,
            max_t=float(np.max(np.asarray(plane.right_edge) - np.asarray(plane.left_edge))),
        )
    return image.reshape(image_shape)


def generate_detector_points(plane: RadiationPlane, *, points_module=None) -> tuple[np.ndarray, tuple[int, int]]:
    """Return detector sample points and output image shape for a plane.

    Scalar ``width``/``resolution`` keeps the original square-plane behavior.
    Two-element values create a rectangular detector in the same basis convention
    used by the external ``RST_cart.generate_plane_points`` helper.
    """

    width = np.asarray(plane.width, dtype=float)
    resolution = np.asarray(plane.resolution, dtype=int)
    if width.ndim == 0 and resolution.ndim == 0:
        if points_module is not None:
            points = points_module.generate_plane_points(
                np.asarray(plane.center, dtype=float),
                np.asarray(plane.normal, dtype=float),
                float(width),
                int(resolution),
            )
        else:
            points = _generate_rectangular_plane_points(
                plane.center,
                plane.normal,
                (float(width), float(width)),
                (int(resolution), int(resolution)),
            )
        return points, (int(resolution), int(resolution))

    if width.shape != (2,) or resolution.shape != (2,):
        raise ValueError("plane.width and plane.resolution must both be scalar or length-2 values")
    points = _generate_rectangular_plane_points(
        plane.center,
        plane.normal,
        (float(width[0]), float(width[1])),
        (int(resolution[0]), int(resolution[1])),
    )
    return points, (int(resolution[1]), int(resolution[0]))


def _generate_rectangular_plane_points(center, normal, width, resolution) -> np.ndarray:
    center = np.asarray(center, dtype=float)
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    a = np.array([0.0, 0.0, 1.0])
    if abs(float(np.dot(a, normal))) > 0.99:
        a = np.array([1.0, 0.0, 0.0])
    u = np.cross(normal, a)
    u = u / np.linalg.norm(u)
    v = np.cross(normal, u)
    v = v / np.linalg.norm(v)
    lin_u = np.linspace(-width[0] / 2.0, width[0] / 2.0, resolution[0])
    lin_v = np.linspace(-width[1] / 2.0, width[1] / 2.0, resolution[1])
    uu, vv = np.meshgrid(lin_u, lin_v)
    return center + np.outer(uu.ravel(), u) + np.outer(vv.ravel(), v)


@contextmanager
def import_dda_tools(dda_root: str | Path | None = None):
    """Import the external ``RST_cart`` DDA package.

    ``physical_module.J_cal`` opens AIA responses through a relative path, so
    this context manager temporarily switches to ``dda_root`` while invoking it.
    """

    dda_root = Path(dda_root or os.environ.get("SOLARUTILS_DDA_ROOT", "")).expanduser()
    if not str(dda_root):
        dda_root = None
    if dda_root is not None:
        dda_root = dda_root.resolve()
        if not (dda_root / "RST_cart").exists():
            raise FileNotFoundError(f"RST_cart was not found under DDA root: {dda_root}")
        sys.path.insert(0, str(dda_root))

    cwd = Path.cwd()
    try:
        from RST_cart import DDA as dda
        from RST_cart import generate_plane_points as points
        from RST_cart import physical_module

        if dda_root is not None:
            os.chdir(dda_root)
        yield {"dda": dda, "points": points, "physical_module": physical_module}
    except ImportError as exc:
        message = (
            "Could not import RST_cart. Install the DDA tools or pass "
            "dda_root=/path/to/DDA_method1, or set SOLARUTILS_DDA_ROOT."
        )
        raise ImportError(message) from exc
    finally:
        os.chdir(cwd)
        if dda_root is not None:
            try:
                sys.path.remove(str(dda_root))
            except ValueError:
                pass


def log_intensity_image(image, *, fill_value: float = 0.0) -> np.ndarray:
    """Return ``log10(image)`` for positive values and ``fill_value`` elsewhere."""

    image = np.asarray(image, dtype=float)
    out = np.full(image.shape, fill_value, dtype=float)
    mask = image > 0.0
    out[mask] = np.log10(image[mask])
    return np.nan_to_num(out, nan=fill_value, posinf=fill_value, neginf=fill_value)


def orient_image(image, *, flipud: bool = True, crop_pixels: int = 0) -> np.ndarray:
    """Apply the notebook display convention: optional vertical flip and side crop."""

    out = np.asarray(image)
    if flipud:
        out = np.flipud(out)
    if crop_pixels:
        crop_pixels = int(crop_pixels)
        if crop_pixels < 0:
            raise ValueError("crop_pixels must be non-negative")
        if 2 * crop_pixels >= out.shape[1]:
            raise ValueError("crop_pixels removes the whole image width")
        out = out[:, crop_pixels:-crop_pixels]
    return out


def plot_radiation_image(
    image,
    *,
    wavelength: int | str = 171,
    frame: int | None = None,
    output: str | Path | None = None,
    figsize=(12, 8),
    extent=(-1.5, 4.5, 0.0, 12.0),
    vmin: float | None = 0.0,
    vmax: float | None = None,
    crop_pixels: int = 128,
    image_is_log: bool = False,
    dpi: int = 200,
):
    """Plot one synthetic radiation image and optionally save it."""

    plot_data = np.asarray(image, dtype=float)
    if not image_is_log:
        plot_data = log_intensity_image(plot_data)
    plot_data = orient_image(plot_data, flipud=True, crop_pixels=crop_pixels)

    cmap = _aia_colormap(wavelength)
    fig, ax = plt.subplots(figsize=figsize, dpi=dpi)
    im = ax.imshow(
        plot_data,
        cmap=cmap,
        origin="lower",
        extent=extent,
        vmin=vmin,
        vmax=vmax if vmax is not None else float(np.nanmax(plot_data)),
    )
    label = f"AIA {wavelength} A" if isinstance(wavelength, int) else str(wavelength)
    if frame is not None:
        label = f"{label} - Frame {frame:04d}"
    ax.set_title(label)
    ax.set_xlabel("X [Mm]")
    ax.set_ylabel("Y [Mm]")
    ax.grid(True, alpha=0.2, linestyle="--")
    cbar = fig.colorbar(im, ax=ax, shrink=1.0, pad=0.02)
    cbar.set_label("log10(Intensity)", rotation=270, labelpad=15)
    fig.tight_layout()

    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(output, bbox_inches="tight")
    return fig, ax


def sample_pixel_line(image, start_pixel, end_pixel, n_samples: int | None = None, *, order: int = 1):
    """Sample a 2D array along a pixel line."""

    from solarutils.obs.time_distance import sample_line_data

    x0, y0 = start_pixel
    x1, y1 = end_pixel
    if n_samples is None:
        n_samples = max(2, int(np.hypot(x1 - x0, y1 - y0)))
    return sample_line_data(image, (x0, y0), (x1, y1), n_samples, order=order)


def time_distance_from_npy(
    paths,
    start_pixel,
    end_pixel,
    *,
    baseline=None,
    crop_pixels: int = 128,
    image_is_oriented: bool = False,
    output: str | Path | None = None,
):
    """Build a time-distance array from saved ``.npy`` radiation images."""

    paths = [Path(path) for path in paths]
    if not paths:
        raise ValueError("paths must contain at least one image")
    if baseline is not None:
        baseline = np.asarray(baseline, dtype=float)

    rows = []
    for path in paths:
        image = np.load(path)
        if not image_is_oriented:
            image = orient_image(image, flipud=True, crop_pixels=crop_pixels)
        if baseline is not None:
            image = image - baseline
        rows.append(sample_pixel_line(image, start_pixel, end_pixel))

    data = np.vstack(rows)
    if output is not None:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            output,
            paths=np.asarray([str(path) for path in paths]),
            start_pixel=np.asarray(start_pixel, dtype=float),
            end_pixel=np.asarray(end_pixel, dtype=float),
            data=data,
        )
    return data


def rectangular_flux_series(paths, *, y_slice: tuple[int, int], x_slice: tuple[int, int]) -> np.ndarray:
    """Return summed image intensity inside one pixel rectangle for each ``.npy`` file."""

    y0, y1 = y_slice
    x0, x1 = x_slice
    values = []
    for path in paths:
        image = np.load(path)
        values.append(float(np.nansum(orient_image(image, flipud=True)[y0:y1, x0:x1])))
    return np.asarray(values, dtype=float)


def _aia_colormap(wavelength: int | str):
    if not isinstance(wavelength, int):
        return "gray"
    cmap_name = f"sdoaia{wavelength}"
    try:
        import sunpy.visualization.colormaps  # noqa: F401
    except ImportError:
        pass
    try:
        plt.colormaps[cmap_name]
    except KeyError:
        return "gray"
    return cmap_name

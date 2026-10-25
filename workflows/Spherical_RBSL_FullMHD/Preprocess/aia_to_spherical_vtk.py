#!/usr/bin/env python
"""Sample an AIA map on a spherical lon/lat grid and write a VTS surface."""

from __future__ import annotations

import argparse
from pathlib import Path

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import SkyCoord
from pyevtk.hl import gridToVTK
from scipy.ndimage import map_coordinates
from sunpy.coordinates import frames
from sunpy.map import Map

try:
    from .rbsl_core.spherical import lon_lat_centers, lon_lat_to_cartesian
except ImportError:  # Direct execution: python Preprocess/aia_to_spherical_vtk.py
    from rbsl_core.spherical import lon_lat_centers, lon_lat_to_cartesian


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)

    aia_map = Map(args.aia)
    lon, lat = lon_lat_centers(
        tuple(args.lon_range),
        tuple(args.lat_range),
        args.resolution_degree,
        nghost=args.nghost,
    )
    data = sample_aia_on_stonyhurst(
        aia_map,
        lon,
        lat,
        exposure_normalize=args.exposure_normalize,
    )
    if args.log10:
        data = np.log10(np.clip(data, a_min=0.0, a_max=None) + 1.0)

    vtk_path = write_spherical_aia_vts(
        output,
        data,
        lon,
        lat,
        radius=args.radius,
        field_name=args.field_name,
    )
    if args.quicklook:
        plot_quicklook(data, lon, lat, output.with_suffix(".png"), log_scaled=args.log10)

    print(f"input: {args.aia}")
    print(f"grid: nlon={lon.size} nlat={lat.size} radius={args.radius}")
    print(f"lon range: {lon[0]:.6g} .. {lon[-1]:.6g}")
    print(f"lat range: {lat[0]:.6g} .. {lat[-1]:.6g}")
    print(f"wrote: {vtk_path}")
    return 0


def parse_args(argv: list[str] | None = None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--aia", required=True, help="Input AIA FITS path.")
    parser.add_argument("--output", required=True, help="Output VTS path, with or without .vts extension.")
    parser.add_argument("--lon-range", nargs=2, type=float, metavar=("LON0", "LON1"), required=True)
    parser.add_argument("--lat-range", nargs=2, type=float, metavar=("LAT0", "LAT1"), required=True)
    parser.add_argument("--resolution-degree", type=float, default=0.2)
    parser.add_argument("--nghost", type=int, default=0)
    parser.add_argument("--radius", type=float, default=1.001)
    parser.add_argument("--field-name", default="aia_intensity")
    parser.add_argument("--exposure-normalize", action="store_true", help="Divide AIA data by exposure_time.")
    parser.add_argument("--log10", action="store_true", help="Write log10(data + 1) instead of raw intensity.")
    parser.add_argument("--quicklook", action="store_true", help="Also write a PNG lon/lat quicklook.")
    return parser.parse_args(argv)


def sample_aia_on_stonyhurst(aia_map, lon_degree, lat_degree, *, exposure_normalize: bool = False):
    """Interpolate an AIA map onto Stonyhurst lon/lat centers."""

    lon2, lat2 = np.meshgrid(lon_degree, lat_degree)
    coords_hgs = SkyCoord(
        lon2.ravel() * u.degree,
        lat2.ravel() * u.degree,
        frame=frames.HeliographicStonyhurst,
        obstime=aia_map.date,
        rsun=aia_map.rsun_meters,
        observer="earth",
    )
    coords_hpc = coords_hgs.transform_to(aia_map.coordinate_frame)
    pixels = aia_map.world_to_pixel(coords_hpc)
    pixel_coords = np.array([pixels[1].value, pixels[0].value])
    source = np.asarray(aia_map.data, dtype=float)
    if exposure_normalize:
        source = source / aia_map.exposure_time.to_value(u.s)
    sampled = map_coordinates(source, pixel_coords, order=1, mode="constant", cval=np.nan)
    return sampled.reshape(lat2.shape)


def write_spherical_aia_vts(
    output,
    data,
    lon_degree,
    lat_degree,
    *,
    radius: float,
    field_name: str,
):
    """Write AIA lon/lat data as a Cartesian spherical-surface `.vts`."""

    output = Path(output)
    path_base = output.with_suffix("") if output.suffix.lower() == ".vts" else output
    data = np.asarray(data, dtype=np.float64)
    lon_degree = np.asarray(lon_degree, dtype=float)
    lat_degree = np.asarray(lat_degree, dtype=float)
    if data.shape != (lat_degree.size, lon_degree.size):
        raise ValueError("data shape must be (len(lat_degree), len(lon_degree))")

    lon2, lat2 = np.meshgrid(lon_degree, lat_degree)
    x, y, z = lon_lat_to_cartesian(lon2, lat2, radius=radius)
    vtk_data = data[None, :, :]
    return gridToVTK(
        str(path_base),
        x[None, :, :].astype(np.float64),
        y[None, :, :].astype(np.float64),
        z[None, :, :].astype(np.float64),
        pointData={field_name: vtk_data},
    )


def plot_quicklook(data, lon, lat, output, *, log_scaled: bool):
    fig, ax = plt.subplots(figsize=(7.0, 4.2), dpi=180)
    image = ax.imshow(
        data,
        origin="lower",
        extent=[lon[0], lon[-1], lat[0], lat[-1]],
        cmap="sdoaia193" if "sdoaia193" in plt.colormaps() else "gray",
        aspect="auto",
    )
    ax.set_xlabel("Stonyhurst longitude [deg]")
    ax.set_ylabel("Stonyhurst latitude [deg]")
    ax.set_title("AIA on spherical lon/lat grid" + (" (log10)" if log_scaled else ""))
    fig.colorbar(image, ax=ax, label="AIA intensity")
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())

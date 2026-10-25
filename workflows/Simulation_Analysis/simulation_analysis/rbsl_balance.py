"""RBSL thin-torus balance-current estimates on spherical AMR backgrounds."""

from __future__ import annotations

import gc
import mmap
import struct
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage import map_coordinates
from scipy.signal import savgol_filter


_VTU_RAW_TAG = b'<AppendedData encoding="raw">'


def read_rbsl_axis(path: str | Path) -> np.ndarray:
    """Read an RBSL axis stored as ASCII or raw float64 lon/lat/r triples."""
    source = Path(path).expanduser().resolve()
    try:
        axis = np.loadtxt(source)
        if axis.ndim == 2 and axis.shape[1] == 3:
            return axis
    except (UnicodeDecodeError, ValueError):
        pass
    return np.fromfile(source, dtype=np.float64).reshape(-1, 3)


def spherical_unit_vectors(theta: np.ndarray, phi: np.ndarray):
    """Return Cartesian spherical basis vectors for arrays of angles."""
    sin_theta, cos_theta = np.sin(theta), np.cos(theta)
    sin_phi, cos_phi = np.sin(phi), np.cos(phi)
    e_r = np.column_stack(
        (sin_theta * cos_phi, sin_theta * sin_phi, cos_theta)
    )
    e_theta = np.column_stack(
        (cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta)
    )
    e_phi = np.column_stack((-sin_phi, cos_phi, np.zeros_like(phi)))
    return e_r, e_theta, e_phi


def rbsl_current_from_flux(
    flux_mx: float,
    minor_radius_mm: float,
    unit_length_mm: float,
    unit_b_gauss: float,
    helicity_sign: int = -1,
) -> float:
    """Convert axial flux to the normalized current used by RBSL_flux_rope."""
    if helicity_sign not in {-1, 1}:
        raise ValueError("helicity_sign must be -1 or 1")
    flux_unit_mx = unit_b_gauss * (unit_length_mm * 1.0e8) ** 2
    flux_code = float(flux_mx) / flux_unit_mx
    radius_code = float(minor_radius_mm) / unit_length_mm
    return float(
        helicity_sign
        * 5.0
        * np.sqrt(2.0)
        * flux_code
        / (12.0 * np.pi * radius_code)
    )


def flux_from_rbsl_current(
    current_code: np.ndarray | float,
    minor_radius_mm: float,
    unit_length_mm: float,
    unit_b_gauss: float,
) -> np.ndarray:
    """Return the positive axial-flux magnitude represented by an RBSL current."""
    flux_unit_mx = unit_b_gauss * (unit_length_mm * 1.0e8) ** 2
    radius_code = float(minor_radius_mm) / unit_length_mm
    conversion = flux_unit_mx * 12.0 * np.pi * radius_code / (5.0 * np.sqrt(2.0))
    return np.abs(np.asarray(current_code)) * conversion


def thin_torus_balance_current(
    strap_field_code: np.ndarray,
    major_radius_mm: float,
    minor_radius_mm: float,
    unit_length_mm: float,
) -> tuple[np.ndarray, float]:
    """Evaluate the local thin-torus balance current along an axis."""
    denominator = np.log(8.0 * major_radius_mm / minor_radius_mm) - 1.0
    if denominator <= 0.0:
        raise ValueError("Thin-torus denominator must be positive")
    major_radius_code = major_radius_mm / unit_length_mm
    return major_radius_code * np.asarray(strap_field_code) / denominator, float(denominator)


def _resolve_magnetic_indices(dataset, config: dict) -> list[int]:
    lower = {str(name).lower(): index for index, name in enumerate(dataset.wnames)}
    indices = []
    for logical in ("b1", "b2", "b3"):
        candidates = config["fields"][logical]
        if isinstance(candidates, str):
            candidates = [candidates]
        for candidate in candidates:
            if candidate.lower() in lower:
                indices.append(lower[candidate.lower()])
                break
        else:
            raise KeyError(
                f"Cannot resolve {logical}; tried {candidates}, DAT fields={dataset.wnames}"
            )
    return indices


def sample_amr_trilinear(dataset, points: np.ndarray, field_indices: list[int]) -> np.ndarray:
    """Sample selected fields at arbitrary points using AMR ghost-cell interpolation."""
    if int(dataset.ghost_width) < 1:
        raise ValueError("Dataset must be opened with ghost_width >= 1")
    dataset.load_data(field_indices=field_indices)
    blocks = dataset.blocks(include_ghosts=True)
    nodes = np.asarray(dataset.mesh.rnode)
    lower = nodes[:, :3]
    upper = nodes[:, 3:6]
    spacing = nodes[:, 6:9]
    sampled = np.full((len(points), len(field_indices)), np.nan)

    for point_index, point in enumerate(points):
        matches = np.flatnonzero(np.all((point >= lower) & (point < upper), axis=1))
        if not matches.size:
            continue
        block_index = int(matches[0])
        cell_coordinate = (
            (point - lower[block_index]) / spacing[block_index]
            + dataset.ghost_width
            - 0.5
        )
        sampled[point_index] = [
            map_coordinates(
                blocks[block_index, field_position],
                cell_coordinate[:, None],
                order=1,
                mode="nearest",
            )[0]
            for field_position in range(len(field_indices))
        ]
    return sampled


def _clustered_unique(values: np.ndarray, tolerance: float = 2.0e-6) -> np.ndarray:
    ordered = np.sort(np.asarray(values, dtype=np.float64))
    starts = np.r_[0, np.flatnonzero(np.diff(ordered) > tolerance) + 1]
    ends = np.r_[starts[1:], len(ordered)]
    return np.asarray(
        [ordered[start:end].mean() for start, end in zip(starts, ends)]
    )


def _read_vtu_data_array(
    mapped_file: mmap.mmap,
    binary_base: int,
    element: ET.Element,
) -> np.ndarray:
    offset = int(element.attrib["offset"])
    position = binary_base + offset
    byte_count = struct.unpack_from("<I", mapped_file, position)[0]
    dtype = {
        "Float32": "<f4",
        "Float64": "<f8",
        "Int32": "<i4",
        "UInt8": "u1",
    }[element.attrib["type"]]
    return np.frombuffer(
        mapped_file,
        dtype=dtype,
        count=byte_count // np.dtype(dtype).itemsize,
        offset=position + 4,
    ).copy()


def _parse_raw_vtu_header(mapped_file: mmap.mmap) -> tuple[ET.Element, int]:
    tag_start = mapped_file.find(_VTU_RAW_TAG)
    if tag_start < 0:
        raise ValueError("VTU is not an AMRVAC raw-appended binary file")
    marker = mapped_file.find(b"_", tag_start + len(_VTU_RAW_TAG))
    if marker < 0:
        raise ValueError("VTU raw-appended data marker was not found")
    root = ET.fromstring(mapped_file[:tag_start] + b"</VTKFile>")
    return root, marker + 1


def _cartesian_to_spherical(points: np.ndarray) -> np.ndarray:
    radius = np.linalg.norm(points, axis=1)
    theta = np.arccos(np.clip(points[:, 2] / radius, -1.0, 1.0))
    phi = np.mod(np.arctan2(points[:, 1], points[:, 0]), 2.0 * np.pi)
    return np.column_stack((radius, theta, phi))


def _interpolate_vtu_block(
    native_points: np.ndarray,
    field_values: np.ndarray,
    targets: np.ndarray,
) -> np.ndarray:
    axes = [_clustered_unique(native_points[:, dim]) for dim in range(3)]
    if any(len(axis) < 2 for axis in axes):
        raise RuntimeError(f"Invalid VTU block axes: {[len(axis) for axis in axes]}")
    indices = np.column_stack(
        [
            np.argmin(
                np.abs(native_points[:, dim, None] - axis[None, :]), axis=1
            )
            for dim, axis in enumerate(axes)
        ]
    )
    grid = np.full(
        tuple(len(axis) for axis in axes) + (field_values.shape[1],), np.nan
    )
    grid[indices[:, 0], indices[:, 1], indices[:, 2], :] = field_values
    if np.any(~np.isfinite(grid)):
        raise RuntimeError(
            "Incomplete structured VTU block: "
            f"axes={[len(axis) for axis in axes]}, "
            f"missing={np.count_nonzero(~np.isfinite(grid[..., 0]))}"
        )
    return RegularGridInterpolator(
        tuple(axes),
        grid,
        method="linear",
        bounds_error=False,
        fill_value=np.nan,
    )(targets)


def sample_raw_vtu_cartesian(
    vtu_path: str | Path,
    spherical_points: np.ndarray,
) -> tuple[np.ndarray, dict]:
    """Sample Cartesian B from an AMRVAC raw-binary VTU on its true grid."""
    source = Path(vtu_path).expanduser().resolve()
    sampled = np.full((len(spherical_points), 3), np.nan)
    source_piece = np.full(len(spherical_points), -1, dtype=int)
    matched_pieces = 0

    with source.open("rb") as handle:
        mapped_file = mmap.mmap(handle.fileno(), 0, access=mmap.ACCESS_READ)
        try:
            root, binary_base = _parse_raw_vtu_header(mapped_file)
            pieces = root.findall("./UnstructuredGrid/Piece")
            for piece_index, piece in enumerate(pieces):
                points_element = piece.find("./Points/DataArray")
                point_data = {
                    element.attrib.get("Name"): element
                    for element in piece.findall("./PointData/DataArray")
                }
                if points_element is None or not all(
                    name in point_data for name in ("b1", "b2", "b3")
                ):
                    continue

                points = _read_vtu_data_array(
                    mapped_file, binary_base, points_element
                ).reshape(-1, 3)
                native = _cartesian_to_spherical(points)
                lower = native.min(axis=0) - 3.0e-6
                upper = native.max(axis=0) + 3.0e-6
                candidates = np.flatnonzero(
                    (source_piece < 0)
                    & np.all(spherical_points >= lower, axis=1)
                    & np.all(spherical_points <= upper, axis=1)
                )
                if not len(candidates):
                    continue

                fields = np.column_stack(
                    [
                        _read_vtu_data_array(
                            mapped_file, binary_base, point_data[name]
                        )
                        for name in ("b1", "b2", "b3")
                    ]
                )
                values = _interpolate_vtu_block(
                    native, fields, spherical_points[candidates]
                )
                finite = np.all(np.isfinite(values), axis=1)
                sampled[candidates[finite]] = values[finite]
                source_piece[candidates[finite]] = piece_index
                matched_pieces += 1
        finally:
            mapped_file.close()

    return sampled, {
        "vtu_piece_count": int(len(pieces)),
        "matched_vtu_pieces": int(matched_pieces),
        "sampled_points": int(np.count_nonzero(source_piece >= 0)),
        "unsampled_points": int(np.count_nonzero(source_piece < 0)),
    }


def estimate_rbsl_balance(
    field_path: str | Path,
    path_dat: str | Path,
    path_info: dict,
    config: dict,
    *,
    minor_radius_mm: float,
    current_flux_mx: float | None = None,
    helicity_sign: int = -1,
    central_fraction: tuple[float, float] = (0.1, 0.9),
) -> tuple[dict, dict[str, np.ndarray]]:
    """Sample an external AMR field and estimate its RBSL balance-current profile.

    Uniform-r DAT files are read through simesh. Geometrically stretched grids
    must be supplied as AMRVAC VTU so the converter's true coordinates are used.
    """

    lower_fraction, upper_fraction = central_fraction
    if not 0.0 <= lower_fraction < upper_fraction <= 1.0:
        raise ValueError("central_fraction must satisfy 0 <= lower < upper <= 1")

    solar_radius_mm = float(path_info["solar_radius_mm"])
    depth_mm = float(path_info["depth_mm"])
    major_radius_mm = float(path_info["major_radius_mm"])
    unit_length_mm = float(config["units"]["length_cm"]) / 1.0e8
    unit_b_gauss = float(config["units"]["magnetic_field_G"])

    raw_axis = read_rbsl_axis(path_dat)
    phi = np.radians(raw_axis[:, 0])
    latitude = np.radians(raw_axis[:, 1])
    theta = 0.5 * np.pi - latitude
    radius = raw_axis[:, 2]
    e_r, e_theta, e_phi = spherical_unit_vectors(theta, phi)
    sample_points = np.column_stack(
        (radius, theta, np.mod(phi, 2.0 * np.pi))
    )
    source = Path(field_path).expanduser().resolve()
    dataset = None
    source_metadata: dict = {}

    if source.suffix.lower() == ".vtu":
        b_cartesian_code, source_metadata = sample_raw_vtu_cartesian(
            source, sample_points
        )
        b_spherical_code = np.column_stack(
            (
                np.sum(b_cartesian_code * e_r, axis=1),
                np.sum(b_cartesian_code * e_theta, axis=1),
                np.sum(b_cartesian_code * e_phi, axis=1),
            )
        )
        sampling_description = (
            "AMRVAC true-coordinate raw VTU block trilinear interpolation"
        )
    elif source.suffix.lower() == ".dat":
        if config["grid"]["radial_type"] != "uniform":
            raise ValueError(
                "Direct DAT balance sampling is unsafe for a geometrically "
                "stretched radial grid because current simesh reconstructs block "
                "radii uniformly. Convert the DAT with AMRVAC -convert and pass "
                "the resulting raw-binary VTU with --vtu."
            )
        from simesh.amrvac import open_dataset

        dataset = open_dataset(str(source), ghost_width=2)
        if int(dataset.ndim) != 3 or "spherical" not in str(dataset.geometry).lower():
            raise ValueError(
                "RBSL balance analysis requires a 3D spherical AMRVAC DAT"
            )
        magnetic_indices = _resolve_magnetic_indices(dataset, config)
        b_spherical_code = sample_amr_trilinear(
            dataset, sample_points, magnetic_indices
        )
        b_cartesian_code = (
            b_spherical_code[:, 0, None] * e_r
            + b_spherical_code[:, 1, None] * e_theta
            + b_spherical_code[:, 2, None] * e_phi
        )
        source_metadata = {
            "dat_time_code": float(dataset.metadata.get("time", 0.0)),
            "dat_iteration": int(dataset.metadata.get("it", 0)),
            "amr_max_level": int(dataset.levmax),
            "amr_leaf_blocks": int(dataset.nleafs),
        }
        sampling_description = (
            "AMR leaf-block trilinear interpolation with two ghost cells"
        )
    else:
        raise ValueError("field_path must end in .dat or .vtu")

    position_mm = radius[:, None] * solar_radius_mm * e_r

    endpoint_midpoint = e_r[0] + e_r[-1]
    endpoint_midpoint /= np.linalg.norm(endpoint_midpoint)
    circle_center_mm = (solar_radius_mm - depth_mm) * endpoint_midpoint

    window = min(201, len(position_mm) // 2 * 2 - 1)
    smooth_position_mm = savgol_filter(
        position_mm, window_length=window, polyorder=3, axis=0
    )
    segments_mm = np.linalg.norm(np.diff(smooth_position_mm, axis=0), axis=1)
    arc_length_mm = np.concatenate(([0.0], np.cumsum(segments_mm)))
    tangent = np.gradient(smooth_position_mm, arc_length_mm, axis=0)
    tangent /= np.linalg.norm(tangent, axis=1)[:, None]
    inward = circle_center_mm - smooth_position_mm
    inward -= np.sum(inward * tangent, axis=1)[:, None] * tangent
    inward /= np.linalg.norm(inward, axis=1)[:, None]
    binormal = np.cross(tangent, inward)
    binormal /= np.linalg.norm(binormal, axis=1)[:, None]

    strap_field_code = -np.sum(b_cartesian_code * binormal, axis=1)
    balance_current_code, denominator = thin_torus_balance_current(
        strap_field_code, major_radius_mm, minor_radius_mm, unit_length_mm
    )
    local_balance_flux_mx = flux_from_rbsl_current(
        balance_current_code, minor_radius_mm, unit_length_mm, unit_b_gauss
    )

    fraction = arc_length_mm / arc_length_mm[-1]
    valid = (
        np.all(np.isfinite(b_spherical_code), axis=1)
        & np.isfinite(balance_current_code)
        & (fraction >= lower_fraction)
        & (fraction <= upper_fraction)
    )
    if not np.any(valid):
        raise RuntimeError("No valid samples in the requested central axis interval")
    weights = np.gradient(arc_length_mm)[valid]
    currents = balance_current_code[valid]
    strap_gauss = strap_field_code[valid] * unit_b_gauss
    global_current = float(np.average(currents, weights=weights))
    global_strap_gauss = float(np.average(strap_gauss, weights=weights))
    global_flux_mx = float(
        flux_from_rbsl_current(
            global_current, minor_radius_mm, unit_length_mm, unit_b_gauss
        )
    )
    apex = int(np.argmax(radius))

    current_code = None
    if current_flux_mx is not None:
        current_code = rbsl_current_from_flux(
            current_flux_mx,
            minor_radius_mm,
            unit_length_mm,
            unit_b_gauss,
            helicity_sign,
        )

    summary = {
        "field_source": str(source),
        "sampling": sampling_description,
        **source_metadata,
        "path_dat": str(Path(path_dat).expanduser().resolve()),
        "path_points": int(len(raw_axis)),
        "solar_radius_mm": solar_radius_mm,
        "depth_mm": depth_mm,
        "major_radius_mm": major_radius_mm,
        "minor_radius_mm": float(minor_radius_mm),
        "r_over_a": major_radius_mm / minor_radius_mm,
        "log_denominator": denominator,
        "central_fraction": [lower_fraction, upper_fraction],
        "valid_points": int(np.count_nonzero(valid)),
        "fraction_with_balance_current_matching_helicity": float(
            np.mean(np.sign(currents) == helicity_sign)
        ),
        "b_strap_gauss_global_mean": global_strap_gauss,
        "b_strap_gauss_apex": float(strap_field_code[apex] * unit_b_gauss),
        "b_strap_gauss_median": float(np.median(strap_gauss)),
        "i_balance_code_global_mean": global_current,
        "i_balance_code_apex": float(balance_current_code[apex]),
        "i_balance_code_median": float(np.median(currents)),
        "i_balance_code_p16": float(np.percentile(currents, 16.0)),
        "i_balance_code_p84": float(np.percentile(currents, 84.0)),
        "flux_balance_global_mx": global_flux_mx,
        "flux_balance_apex_mx": float(local_balance_flux_mx[apex]),
        "flux_balance_median_mx": float(
            flux_from_rbsl_current(
                np.median(currents), minor_radius_mm, unit_length_mm, unit_b_gauss
            )
        ),
    }
    if current_flux_mx is not None and current_code is not None:
        summary.update(
            {
                "current_flux_mx": float(current_flux_mx),
                "helicity_sign": int(helicity_sign),
                "i_current_code": current_code,
                "i_current_over_global_balance": current_code / global_current,
                "current_flux_over_balance_global": current_flux_mx / global_flux_mx,
            }
        )

    profile = {
        "index": np.arange(len(radius)),
        "fraction": fraction,
        "arc_length_mm": arc_length_mm,
        "phi_input_deg": raw_axis[:, 0],
        "latitude_deg": raw_axis[:, 1],
        "radius_rsun": radius,
        "br_gauss": b_spherical_code[:, 0] * unit_b_gauss,
        "btheta_gauss": b_spherical_code[:, 1] * unit_b_gauss,
        "bphi_gauss": b_spherical_code[:, 2] * unit_b_gauss,
        "b_strap_gauss": strap_field_code * unit_b_gauss,
        "i_balance_code": balance_current_code,
        "flux_balance_local_mx": local_balance_flux_mx,
    }
    if dataset is not None:
        del dataset
    gc.collect()
    return summary, profile

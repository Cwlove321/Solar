"""Helpers for spherical simulation grids."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.ndimage import laplace
from scipy.interpolate import RegularGridInterpolator


DEGREE = np.pi / 180.0


@dataclass(frozen=True)
class SphericalGrid:
    """A spherical simulation grid with explicit axis arrays.

    `theta` is colatitude in radians unless documented otherwise by a workflow.
    The expected data shape is `(r, theta, phi)`.
    """

    r: np.ndarray
    theta: np.ndarray
    phi: np.ndarray
    radius_unit: str = ""
    angle_unit: str = "rad"

    @property
    def shape(self) -> tuple[int, int, int]:
        """Return the expected data shape in `(r, theta, phi)` order."""

        return (len(self.r), len(self.theta), len(self.phi))

    @property
    def latitude(self) -> np.ndarray:
        """Return latitude values derived from colatitude."""

        theta = np.asarray(self.theta, dtype=float)
        return np.pi / 2.0 - theta


def spherical_to_cartesian(r, theta, phi):
    """Convert spherical coordinates to Cartesian coordinates.

    `theta` is colatitude and `phi` is longitude, both in radians.
    """

    r = np.asarray(r, dtype=float)
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def lon_lat_to_cartesian(lon, lat, radius=1.0, *, degrees: bool = True):
    """Convert longitude/latitude coordinates to Cartesian coordinates."""

    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    radius = np.asarray(radius, dtype=float)
    if degrees:
        lon = np.radians(lon)
        lat = np.radians(lat)

    cos_lat = np.cos(lat)
    x = radius * cos_lat * np.cos(lon)
    y = radius * cos_lat * np.sin(lon)
    z = radius * np.sin(lat)
    return x, y, z


def cartesian_to_lon_lat_r(x, y, z, *, degrees: bool = True):
    """Convert Cartesian coordinates to longitude/latitude/radius."""

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    z = np.asarray(z, dtype=float)
    radius = np.sqrt(x**2 + y**2 + z**2)
    safe_radius = np.where(radius == 0.0, np.nan, radius)
    lon = np.arctan2(y, x)
    lat = np.arcsin(np.clip(z / safe_radius, -1.0, 1.0))
    if degrees:
        lon = np.degrees(lon)
        lat = np.degrees(lat)
    return lon, lat, radius


def flux_rope_boundary_fieldline(path_lonlatr, rope_radius_mm, twist, *, solar_radius_mm: float = 696.0):
    """Return one helical boundary field line around an axis path.

    `path_lonlatr` has columns `(lon_deg, lat_deg, radius_rsun)`. The returned
    array has the same column convention. The first/last axis points are kept as
    the phase endpoints, while the boundary itself is offset by `rope_radius_mm`.
    """

    path_lonlatr = np.asarray(path_lonlatr, dtype=float)
    if path_lonlatr.ndim != 2 or path_lonlatr.shape[1] != 3:
        raise ValueError("path_lonlatr must have shape (N, 3)")
    if path_lonlatr.shape[0] < 2:
        raise ValueError("path_lonlatr must contain at least two points")
    rope_radius_mm = float(rope_radius_mm)
    twist = float(twist)
    solar_radius_mm = float(solar_radius_mm)
    if rope_radius_mm <= 0.0:
        raise ValueError("rope_radius_mm must be positive")
    if solar_radius_mm <= 0.0:
        raise ValueError("solar_radius_mm must be positive")

    x, y, z = lon_lat_to_cartesian(
        path_lonlatr[:, 0],
        path_lonlatr[:, 1],
        path_lonlatr[:, 2] * solar_radius_mm,
    )
    axis = np.column_stack([x, y, z])
    tangent = np.gradient(axis, axis=0)
    tangent = _normalize_vectors(tangent)

    radial = _normalize_vectors(axis)
    normal = radial - np.sum(radial * tangent, axis=1)[:, None] * tangent
    degenerate = np.linalg.norm(normal, axis=1) < 1e-10
    if np.any(degenerate):
        fallback = np.tile(np.array([0.0, 0.0, 1.0]), (axis.shape[0], 1))
        fallback_parallel = np.abs(np.sum(fallback * tangent, axis=1)) > 0.95
        fallback[fallback_parallel] = np.array([0.0, 1.0, 0.0])
        normal[degenerate] = fallback[degenerate] - (
            np.sum(fallback[degenerate] * tangent[degenerate], axis=1)[:, None] * tangent[degenerate]
        )
    normal = _normalize_vectors(normal)
    binormal = _normalize_vectors(np.cross(tangent, normal))

    segment_lengths = np.linalg.norm(np.diff(axis, axis=0), axis=1)
    distance = np.concatenate([[0.0], np.cumsum(segment_lengths)])
    phase_fraction = np.zeros_like(distance) if distance[-1] == 0.0 else distance / distance[-1]
    phase = 2.0 * np.pi * twist * phase_fraction

    boundary = axis + rope_radius_mm * (np.cos(phase)[:, None] * normal + np.sin(phase)[:, None] * binormal)
    lon, lat, radius_mm = cartesian_to_lon_lat_r(boundary[:, 0], boundary[:, 1], boundary[:, 2])
    return np.column_stack([lon, lat, radius_mm / solar_radius_mm])


def _normalize_vectors(vectors):
    vectors = np.asarray(vectors, dtype=float)
    norms = np.linalg.norm(vectors, axis=1)
    safe_norms = np.where(norms == 0.0, 1.0, norms)
    return vectors / safe_norms[:, None]


def angular_distance(lon1, lat1, lon2, lat2, *, degrees: bool = True):
    """Return great-circle angular distance between lon/lat points."""

    lon1 = np.asarray(lon1, dtype=float)
    lat1 = np.asarray(lat1, dtype=float)
    lon2 = np.asarray(lon2, dtype=float)
    lat2 = np.asarray(lat2, dtype=float)
    if degrees:
        lon1 = np.radians(lon1)
        lat1 = np.radians(lat1)
        lon2 = np.radians(lon2)
        lat2 = np.radians(lat2)

    cos_distance = np.sin(lat1) * np.sin(lat2) + np.cos(lat1) * np.cos(lat2) * np.cos(lon2 - lon1)
    distance = np.arccos(np.clip(cos_distance, -1.0, 1.0))
    return np.degrees(distance) if degrees else distance


def great_circle_points(lon1, lat1, lon2, lat2, n: int, *, degrees: bool = True):
    """Return `n` equally spaced great-circle points between two lon/lat points."""

    if n < 1:
        raise ValueError("n must be at least 1")
    if n == 1:
        return np.array([[lon1, lat1]], dtype=float)

    lon1_rad, lat1_rad, lon2_rad, lat2_rad = (lon1, lat1, lon2, lat2)
    if degrees:
        lon1_rad, lat1_rad, lon2_rad, lat2_rad = np.radians([lon1, lat1, lon2, lat2])

    p1 = np.array(lon_lat_to_cartesian(lon1_rad, lat1_rad, degrees=False))
    p2 = np.array(lon_lat_to_cartesian(lon2_rad, lat2_rad, degrees=False))
    sigma = float(np.arccos(np.clip(np.dot(p1, p2), -1.0, 1.0)))
    if sigma < 1e-12:
        point = np.array([lon1, lat1], dtype=float)
        return np.repeat(point[None, :], n, axis=0)

    t = np.linspace(0.0, 1.0, n)
    sin_sigma = np.sin(sigma)
    weights1 = np.sin((1.0 - t) * sigma) / sin_sigma
    weights2 = np.sin(t * sigma) / sin_sigma
    points = weights1[:, None] * p1 + weights2[:, None] * p2
    lon = np.arctan2(points[:, 1], points[:, 0])
    lat = np.arcsin(np.clip(points[:, 2] / np.linalg.norm(points, axis=1), -1.0, 1.0))
    if degrees:
        lon = np.degrees(lon)
        lat = np.degrees(lat)
    return np.column_stack([lon, lat])


def resample_great_circle_path(points, n: int, *, degrees: bool = True):
    """Resample a lon/lat polyline to `n` equally spaced great-circle points."""

    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")
    if points.shape[0] < 2:
        raise ValueError("points must contain at least two points")
    if n < 2:
        raise ValueError("n must be at least 2")

    if degrees:
        distances = np.radians(
            angular_distance(points[:-1, 0], points[:-1, 1], points[1:, 0], points[1:, 1])
        )
    else:
        distances = angular_distance(points[:-1, 0], points[:-1, 1], points[1:, 0], points[1:, 1], degrees=False)
    cumulative = np.concatenate([[0.0], np.cumsum(distances)])
    total = float(cumulative[-1])
    if total < 1e-12:
        return np.repeat(points[:1], n, axis=0), 0.0

    targets = np.linspace(0.0, total, n)
    out = np.empty((n, 2), dtype=float)
    for index, target in enumerate(targets):
        segment = np.searchsorted(cumulative, target, side="right") - 1
        segment = min(max(segment, 0), points.shape[0] - 2)
        segment_length = cumulative[segment + 1] - cumulative[segment]
        fraction = 0.0 if segment_length == 0 else (target - cumulative[segment]) / segment_length
        if fraction == 0.0:
            out[index] = points[segment]
        elif fraction == 1.0:
            out[index] = points[segment + 1]
        else:
            out[index] = _great_circle_fraction(points[segment], points[segment + 1], fraction, degrees=degrees)

    return out, np.degrees(total) if degrees else total


def spherical_cap_flux(lon, lat, br, center_lon, center_lat, radius, *, solar_radius=696.0):
    """Estimate flux through a circular cap on a lon/lat `Br` map.

    `radius` and `solar_radius` must use the same length unit. `Br` is assumed to
    be in Gauss, so the returned flux is in Mx when lengths are in Mm.
    """

    lon = np.asarray(lon, dtype=float)
    lat = np.asarray(lat, dtype=float)
    br = np.asarray(br, dtype=float)
    lon2, lat2 = np.meshgrid(lon, lat)
    if br.shape != lon2.shape:
        raise ValueError("br shape must be (len(lat), len(lon))")

    cap_radius = float(radius) / float(solar_radius)
    distance = angular_distance(center_lon, center_lat, lon2, lat2, degrees=False)
    mask = distance <= cap_radius
    if not np.any(mask):
        return {"flux_mx": 0.0, "mean_br": 0.0, "pixel_count": 0}

    mean_br = float(np.nanmean(br[mask]))
    area_cm2 = np.pi * (float(radius) * 1e8) ** 2
    return {
        "flux_mx": float(mean_br * area_cm2),
        "mean_br": mean_br,
        "pixel_count": int(np.count_nonzero(mask)),
    }


def flux_rope_axis_from_points(points, height_mm, *, npoints: int = 6000, solar_radius_mm: float = 696.0):
    """Build a 3D flux-rope axis from clicked lon/lat path points.

    Parameters
    ----------
    points
        `(N, 2)` array of Stonyhurst longitude/latitude points in degrees.
    height_mm
        Apex height above the photosphere in Mm.
    npoints
        Number of points written to `path.dat`; AMRVAC currently expects 6000.
    solar_radius_mm
        Solar radius in Mm.

    Returns
    -------
    dict
        Contains `path` with columns `(lon_deg, lat_deg, radius_rsun)`, plus
        depth/major-radius diagnostics in Mm.
    """

    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2:
        raise ValueError("points must have shape (N, 2)")
    if points.shape[0] < 2:
        raise ValueError("at least two path points are required")
    if npoints < 2:
        raise ValueError("npoints must be at least 2")

    height_mm = float(height_mm)
    solar_radius_mm = float(solar_radius_mm)
    if height_mm <= 0:
        raise ValueError("height_mm must be positive")
    if solar_radius_mm <= 0:
        raise ValueError("solar_radius_mm must be positive")

    path_lonlat, rope_length_degree = resample_great_circle_path(points, npoints)
    endpoint_angle = np.radians(angular_distance(*points[0], *points[-1]))
    half_angle = endpoint_angle / 2.0
    one_minus_cos = 1.0 - np.cos(half_angle)
    depth_mm = (
        2.0 * solar_radius_mm**2 * one_minus_cos - height_mm**2
    ) / (2.0 * (height_mm + solar_radius_mm * one_minus_cos))
    major_radius_mm = depth_mm + height_mm

    axis_angle = np.abs(np.linspace(-half_angle, half_angle, npoints))
    inside = major_radius_mm**2 - ((solar_radius_mm - depth_mm) * np.sin(axis_angle)) ** 2
    if np.any(inside < -1e-8):
        raise ValueError("height/footpoint geometry gives an invalid circular axis")
    radius_rsun = (
        (solar_radius_mm - depth_mm) * np.cos(axis_angle) + np.sqrt(np.clip(inside, 0.0, None))
    ) / solar_radius_mm

    path = np.column_stack([path_lonlat[:, 0], path_lonlat[:, 1], radius_rsun])
    return {
        "path": path,
        "depth_mm": float(depth_mm),
        "major_radius_mm": float(major_radius_mm),
        "endpoint_angle_rad": float(endpoint_angle),
        "projected_length_mm": float(np.radians(rope_length_degree) * solar_radius_mm),
    }


def write_path_dat(path, axis, *, add_lon_180: bool = True):
    """Write AMRVAC `path.dat` as float64 `(lon, lat, radius)` records."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    axis = np.asarray(axis, dtype=float)
    if axis.ndim != 2 or axis.shape[1] != 3:
        raise ValueError("axis must have shape (N, 3)")
    output = axis.copy()
    if add_lon_180:
        output[:, 0] += 180.0
    np.asarray(output, dtype=np.float64).tofile(path)
    return path


def read_path_dat(path, *, subtract_lon_180: bool = False):
    """Read AMRVAC `path.dat` into an `(N, 3)` float64 array."""

    data = np.fromfile(path, dtype=np.float64)
    if data.size % 3:
        raise ValueError("path.dat length is not divisible by 3")
    axis = data.reshape(-1, 3)
    if subtract_lon_180:
        axis = axis.copy()
        axis[:, 0] -= 180.0
    return axis


def _great_circle_fraction(start, end, fraction: float, *, degrees: bool):
    start = np.asarray(start, dtype=float)
    end = np.asarray(end, dtype=float)
    start_rad = np.radians(start) if degrees else start
    end_rad = np.radians(end) if degrees else end
    p1 = np.array(lon_lat_to_cartesian(start_rad[0], start_rad[1], degrees=False))
    p2 = np.array(lon_lat_to_cartesian(end_rad[0], end_rad[1], degrees=False))
    sigma = float(np.arccos(np.clip(np.dot(p1, p2), -1.0, 1.0)))
    if sigma < 1e-12:
        return start.copy()

    sin_sigma = np.sin(sigma)
    point = (
        np.sin((1.0 - fraction) * sigma) / sin_sigma * p1
        + np.sin(fraction * sigma) / sin_sigma * p2
    )
    lon = np.arctan2(point[1], point[0])
    lat = np.arcsin(np.clip(point[2] / np.linalg.norm(point), -1.0, 1.0))
    result = np.array([lon, lat])
    return np.degrees(result) if degrees else result


def hmi_disambiguate_azimuth(azimuth, disambig, method: int = 2):
    """Apply the HMI disambiguation bit mask to an azimuth array.

    `method` follows the HMI convention used in the notebooks: 0, 1, or 2 select
    the corresponding bit from the disambiguation segment. Pixels with a set bit
    are rotated by 180 degrees.
    """

    method = 2 if method < 0 or method > 2 else int(method)
    azimuth = np.asarray(azimuth, dtype=float)
    disambig = np.nan_to_num(disambig, nan=0, posinf=0, neginf=0).astype(np.int32)
    if azimuth.shape != disambig.shape:
        raise ValueError("azimuth and disambig must have the same shape")

    corrected = azimuth.copy()
    corrected[((disambig >> method) & 1).astype(bool)] += 180.0
    return corrected


def hmi_vector_components(field, inclination, azimuth):
    """Convert HMI vector field segments into image-plane `Bx, By, Bz`.

    Inputs are in the usual HMI units: field strength in Gauss and inclination /
    azimuth in degrees. The returned components match the older RBSL notebooks.
    """

    field = np.asarray(field, dtype=float)
    inclination = np.asarray(inclination, dtype=float)
    azimuth = np.asarray(azimuth, dtype=float)
    if field.shape != inclination.shape or field.shape != azimuth.shape:
        raise ValueError("field, inclination, and azimuth must have the same shape")

    inc = np.radians(inclination)
    azi = np.radians(azimuth + 270.0)
    bx = field * np.sin(inc) * np.cos(azi)
    by = field * np.sin(inc) * np.sin(azi)
    bz = field * np.cos(inc)
    return bx, by, bz


def image_magnetic_to_spherical_components(x, y, bx, by, bz, theta_deg, rsun):
    """Project image-plane magnetic vectors to spherical `Br, Btheta, Bphi`.

    `x` and `y` are helioprojective coordinates in arcsec, `theta_deg` is the
    observer latitude offset used by the RBSL notebooks, and `rsun` is the
    apparent solar radius in arcsec.
    """

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    bx = np.asarray(bx, dtype=float)
    by = np.asarray(by, dtype=float)
    bz = np.asarray(bz, dtype=float)

    theta = np.radians(theta_deg)
    rotation = np.array(
        [
            [0.0, -np.cos(theta), -np.sin(theta)],
            [1.0, 0.0, 0.0],
            [0.0, -np.sin(theta), np.cos(theta)],
        ]
    )

    z = np.sqrt(np.maximum(rsun**2 - x**2 - y**2, 0.0))
    position = np.column_stack([x, y, z])
    vector = np.column_stack([bx, by, bz])
    position_rot = position @ rotation.T
    vector_rot = vector @ rotation.T

    r = np.linalg.norm(position_rot, axis=1)
    theta_s = np.arccos(np.clip(position_rot[:, 2] / r, -1.0, 1.0))
    phi_s = np.arctan2(position_rot[:, 1], position_rot[:, 0])

    sin_theta = np.sin(theta_s)
    cos_theta = np.cos(theta_s)
    sin_phi = np.sin(phi_s)
    cos_phi = np.cos(phi_s)
    r_hat = np.column_stack([sin_theta * cos_phi, sin_theta * sin_phi, cos_theta])
    theta_hat = np.column_stack([cos_theta * cos_phi, cos_theta * sin_phi, -sin_theta])
    phi_hat = np.column_stack([-sin_phi, cos_phi, np.zeros_like(phi_s)])

    br = np.einsum("ij,ij->i", vector_rot, r_hat)
    bt = np.einsum("ij,ij->i", vector_rot, theta_hat)
    bp = np.einsum("ij,ij->i", vector_rot, phi_hat)
    return r, np.degrees(theta_s), np.degrees(phi_s), br, bt, bp


def lon_lat_centers(
    lon_range: tuple[float, float],
    lat_range: tuple[float, float],
    resolution_degree: float,
    *,
    nghost: int = 0,
):
    """Return equally spaced longitude/latitude centers in degrees.

    The upper edge is not included, matching `np.linspace(..., endpoint=False)`
    in the source notebooks.
    """

    resolution_degree = float(resolution_degree)
    ghost = int(nghost) * resolution_degree
    lon0, lon1 = lon_range
    lat0, lat1 = lat_range
    nlon = int(np.round((lon1 - lon0 + 2.0 * ghost) / resolution_degree))
    nlat = int(np.round((lat1 - lat0 + 2.0 * ghost) / resolution_degree))
    lon = np.linspace(lon0 - ghost, lon1 + ghost, nlon, endpoint=False)
    lat = np.linspace(lat0 - ghost, lat1 + ghost, nlat, endpoint=False)
    return lon, lat


def spherical_flux_balance(br, latitude_degree):
    """Calculate signed, positive, negative, and unsigned spherical flux weights.

    The returned values are in field-strength times relative area weight. For a
    true physical flux, multiply by the surface area represented by one grid
    cell. `balance` is `signed / unsigned`, matching the notebook diagnostic.
    """

    br = np.asarray(br, dtype=float)
    latitude = np.asarray(latitude_degree, dtype=float)
    if latitude.ndim == 1:
        weights = np.cos(np.radians(latitude))[:, None]
    else:
        weights = np.cos(np.radians(latitude))
    weighted = br * weights
    finite = np.isfinite(weighted)
    positive = float(np.sum(weighted[finite & (weighted > 0)]))
    negative = float(np.sum(weighted[finite & (weighted < 0)]))
    unsigned = float(np.sum(np.abs(weighted[finite])))
    signed = float(np.sum(weighted[finite]))
    balance = signed / unsigned if unsigned else np.nan
    return {
        "positive": positive,
        "negative": negative,
        "unsigned": unsigned,
        "signed": signed,
        "balance": float(balance),
    }


def preprocess_spherical_field(
    br,
    bt,
    bp,
    theta,
    phi,
    *,
    mu3: float = 0.01,
    mu4: float = 0.1,
    scl: float = 1.0,
    max_iter: int = 5000,
    dt: int = 200,
    tol: float = 1e-4,
    progress=None,
):
    """Preprocess a spherical vector magnetogram by minimizing force/torque terms.

    This is a direct, vectorized translation of the RBSL notebook routine. It
    returns `(Br, Btheta, Bphi)` with the original field scale restored.
    """

    br = np.asarray(br, dtype=float)
    bt = np.asarray(bt, dtype=float)
    bp = np.asarray(bp, dtype=float)
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    if br.shape != bt.shape or br.shape != bp.shape:
        raise ValueError("br, bt, and bp must have the same shape")
    if br.shape != (theta.size, phi.size):
        raise ValueError("field shape must be (len(theta), len(phi))")

    mu1 = 1.0 * scl
    mu2 = 1.0 * scl
    mu3 = mu3 * scl
    mu4 = mu4 * scl

    bave = np.sqrt(br**2 + bt**2 + bp**2).mean()
    if not np.isfinite(bave) or bave == 0:
        raise ValueError("mean field strength must be finite and non-zero")
    br0, bt0, bp0 = br / bave, bt / bave, bp / bave
    bro, bto, bpo = br0.copy(), bt0.copy(), bp0.copy()

    theta2, phi2 = _make_mesh(theta, phi)
    sinth = np.sin(theta2)
    eb = np.sum(sinth * (bt0**2 + bp0**2 + br0**2))

    l12_old = 0.0
    l3_old = 0.0
    l4_old = 0.0
    old_l = np.inf

    for iteration in range(max_iter + 1):
        terms = _preprocess_terms(br0, bt0, bp0, theta2, phi2, sinth, eb)
        sin_th = terms["sin_th"]
        cos_th = terms["cos_th"]
        cos_ph = terms["cos_ph"]
        sin_ph = terms["sin_ph"]
        term1a = terms["term1a"]
        term1b = terms["term1b"]
        term1c = terms["term1c"]
        term2a = terms["term2a"]
        term2b = terms["term2b"]
        term2c = terms["term2c"]
        term3a = (br0 - bro) / eb
        term3b = (bt0 - bto) / eb
        term3c = (bp0 - bpo) / eb
        term4a = _laplace2(br0) / eb
        term4b = _laplace2(bt0) / eb
        term4c = _laplace2(bp0) / eb

        l1 = term1a**2 + term1b**2 + term1c**2
        l2 = term2a**2 + term2b**2 + term2c**2
        l12 = l1 + l2
        l3 = np.sum(term3a**2 + term3b**2 + term3c**2)
        l4 = np.sum(term4a**2 + term4b**2 + term4c**2)
        loss = mu1 * l12 + mu3 * l3 + mu4 * l4

        if iteration == 0:
            dloss = 1.0
        else:
            dloss = (
                _safe_relative_change(l12, l12_old)
                + _safe_relative_change(l3, l3_old)
                + _safe_relative_change(l4, l4_old)
            )
        l12_old, l3_old, l4_old = l12, l3, l4

        if progress is not None and (iteration % dt == 0 or iteration == max_iter):
            progress(
                {
                    "iteration": iteration,
                    "dloss": float(dloss),
                    "l12": float(l12),
                    "l3": float(l3),
                    "l4": float(l4),
                    "loss": float(loss),
                    "eps_f": float(np.abs([term1a, term1b, term1c]).sum()),
                    "eps_t": float(np.abs([term2a, term2b, term2c]).sum()),
                }
            )

        if iteration > 0 and dloss < tol and loss < old_l:
            break
        old_l = loss

        br1 = br0 - 2.0 * mu3 * term3a - 2.0 * mu4 * _laplace2(term4a)
        bt1 = (
            bt0
            - 2.0 * mu1 * term1a * (bt0 * sin_th * sin_th * cos_ph - br0 * sin_th * cos_th * cos_ph)
            - 2.0 * mu1 * term1b * (bt0 * sin_th * sin_th * sin_ph - br0 * sin_th * cos_th * sin_ph)
            - 2.0 * mu1 * term1c * (bt0 * sin_th * cos_th + br0 * sin_th * sin_th)
            - 2.0 * mu2 * (term2a * br0 * sin_th * sin_ph - term2b * br0 * sin_th * cos_ph)
            - 2.0 * mu3 * term3b
            - 2.0 * mu4 * _laplace2(term4b)
        )
        bp1 = (
            bp0
            - 2.0 * mu1 * term1a * (bp0 * sin_th * sin_th * cos_ph + br0 * sin_th * sin_ph)
            - 2.0 * mu1 * term1b * (bp0 * sin_th * sin_th * sin_ph - br0 * sin_th * cos_ph)
            - 2.0 * mu1 * term1c * (bp0 * sin_th * cos_th)
            - 2.0
            * mu2
            * (
                term2a * br0 * cos_th * cos_ph * sin_th
                + term2b * br0 * cos_th * sin_ph * sin_th
                + term2c * br0 * sin_th * sin_th
            )
            - 2.0 * mu3 * term3c
            - 2.0 * mu4 * _laplace2(term4c)
        )
        br0, bt0, bp0 = br1, bt1, bp1

    return br0 * bave, bt0 * bave, bp0 * bave


def spherical_preprocess_metrics(br, bt, bp, theta, phi):
    """Return force/torque and smoothing diagnostics used by preprocessing.

    `Br`, `Btheta`, and `Bphi` are normalized internally in the same way as
    `preprocess_spherical_field`, so these metrics can be compared before and
    after preprocessing.
    """

    br = np.asarray(br, dtype=float)
    bt = np.asarray(bt, dtype=float)
    bp = np.asarray(bp, dtype=float)
    theta = np.asarray(theta, dtype=float)
    phi = np.asarray(phi, dtype=float)
    if br.shape != bt.shape or br.shape != bp.shape:
        raise ValueError("br, bt, and bp must have the same shape")
    if br.shape != (theta.size, phi.size):
        raise ValueError("field shape must be (len(theta), len(phi))")

    bave = np.sqrt(br**2 + bt**2 + bp**2).mean()
    if not np.isfinite(bave) or bave == 0:
        raise ValueError("mean field strength must be finite and non-zero")
    br0, bt0, bp0 = br / bave, bt / bave, bp / bave
    theta2, phi2 = _make_mesh(theta, phi)
    sinth = np.sin(theta2)
    eb = np.sum(sinth * (bt0**2 + bp0**2 + br0**2))
    terms = _preprocess_terms(br0, bt0, bp0, theta2, phi2, sinth, eb)
    term4a = _laplace2(br0) / eb
    term4b = _laplace2(bt0) / eb
    term4c = _laplace2(bp0) / eb

    l1 = terms["term1a"] ** 2 + terms["term1b"] ** 2 + terms["term1c"] ** 2
    l2 = terms["term2a"] ** 2 + terms["term2b"] ** 2 + terms["term2c"] ** 2
    l4 = np.sum(term4a**2 + term4b**2 + term4c**2)
    return {
        "l1_force": float(l1),
        "l2_torque": float(l2),
        "l12_force_torque": float(l1 + l2),
        "l4_smoothing": float(l4),
        "eps_force": float(np.abs([terms["term1a"], terms["term1b"], terms["term1c"]]).sum()),
        "eps_torque": float(np.abs([terms["term2a"], terms["term2b"], terms["term2c"]]).sum()),
    }


def resample_mrdailysynframe(data, resolution_degree: float):
    """Resample an HMI `mrdailysynframe_polfil` map onto regular lon/lat centers.

    The HMI map latitude coordinate is uniform in `sin(latitude)`. The returned
    `lon` and `lat` arrays are in degrees and the data shape is `(lat, lon)`.
    """

    data = np.asarray(data, dtype=float)
    ny, nx = data.shape
    coord_sinlat = np.arcsin(np.linspace(-1.0, 1.0, ny))
    coord_lon = (np.arange(nx) - nx / 2.0 + 0.5) * (360.0 / nx) * np.pi / 180.0
    interpolator = RegularGridInterpolator(
        (coord_sinlat, coord_lon),
        data,
        bounds_error=False,
        fill_value=np.nan,
    )

    lon = np.arange(-180.0, 180.0, resolution_degree) + resolution_degree / 2.0
    lat = np.arange(-90.0, 90.0, resolution_degree) + resolution_degree / 2.0
    lon_rad = np.radians(lon)
    lat_rad = np.radians(lat)
    lon_grid, lat_grid = np.meshgrid(lon_rad, lat_rad)
    points = np.column_stack([lat_grid.ravel(), lon_grid.ravel()])
    resampled = interpolator(points).reshape(lat_grid.shape)
    return resampled, lon, lat


def patch_bounds_from_lon_lat(
    lon_range: tuple[float, float],
    lat_range: tuple[float, float],
    resolution_degree: float,
):
    """Return integer `(x0, y0, nx, ny)` for a lon/lat patch on a global grid."""

    lon0, lon1 = lon_range
    lat0, lat1 = lat_range
    x0 = int(np.round((lon0 + 180.0) / resolution_degree))
    y0 = int(np.round((lat0 + 90.0) / resolution_degree))
    nx = int(np.round((lon1 - lon0) / resolution_degree))
    ny = int(np.round((lat1 - lat0) / resolution_degree))
    return x0, y0, nx, ny


def replace_synoptic_patch(global_br, patch_br, x0: int, y0: int):
    """Replace one rectangular patch in a global synoptic map."""

    global_br = np.asarray(global_br, dtype=float).copy()
    patch_br = np.asarray(patch_br, dtype=float)
    ny, nx = patch_br.shape
    if y0 < 0 or x0 < 0 or y0 + ny > global_br.shape[0] or x0 + nx > global_br.shape[1]:
        raise ValueError("patch is outside the global map")
    global_br[y0 : y0 + ny, x0 : x0 + nx] = patch_br
    return global_br


def match_template_location(big, template, *, method: str = "ccoeff_normed"):
    """Return the best `(x, y, score)` location of `template` within `big`.

    This mirrors the OpenCV `matchTemplate` step used in the RBSL notebooks.
    """

    try:
        import cv2
    except ImportError as exc:
        raise ImportError("match_template_location requires OpenCV (`cv2`)") from exc

    method_map = {
        "ccoeff_normed": cv2.TM_CCOEFF_NORMED,
        "sqdiff_normed": cv2.TM_SQDIFF_NORMED,
        "ccorr_normed": cv2.TM_CCORR_NORMED,
    }
    if method not in method_map:
        raise ValueError(f"method must be one of {tuple(method_map)}")

    big = np.nan_to_num(np.asarray(big, dtype=np.float32), nan=0.0)
    template = np.nan_to_num(np.asarray(template, dtype=np.float32), nan=0.0)
    if big.ndim != 2 or template.ndim != 2:
        raise ValueError("big and template must be 2D arrays")
    if template.shape[0] > big.shape[0] or template.shape[1] > big.shape[1]:
        raise ValueError("template must fit inside big")

    result = cv2.matchTemplate(big, template, method_map[method])
    if method_map[method] in (cv2.TM_SQDIFF, cv2.TM_SQDIFF_NORMED):
        min_val, _, min_loc, _ = cv2.minMaxLoc(result)
        return int(min_loc[0]), int(min_loc[1]), float(1.0 - min_val)

    _, max_val, _, max_loc = cv2.minMaxLoc(result)
    return int(max_loc[0]), int(max_loc[1]), float(max_val)


def pfss_mod_usr_parameters(global_shape, x0: int, y0: int, nx: int, ny: int):
    """Return index limits for `pfss_mod_usr.t`."""

    pfss_ny, _ = global_shape
    return {
        "lat0": int(pfss_ny - (y0 + ny - 1) - 1),
        "lat1": int(pfss_ny - y0 - 1),
        "lon0": int(x0),
        "lon1": int(x0 + nx - 1),
    }


def pfss_amrvac_parameters(
    global_shape,
    x0: int,
    y0: int,
    nx: int,
    ny: int,
    *,
    nr: int,
    length_r: float,
    refine_level: int = 1,
):
    """Return domain settings for `pfss_amrvac.par`."""

    pfss_ny, pfss_nx = global_shape
    divisor = 2 ** (int(refine_level) - 1)
    return {
        "domain_nx1": ny / divisor,
        "domain_nx2": ny / divisor,
        "domain_nx3": nx / divisor,
        "xprobmin1": 1.0 + 0.5 * length_r / nr,
        "xprobmax1": 1.0 + 0.5 * length_r / nr + length_r,
        "xprobmin2": (pfss_ny - y0 - ny) / pfss_ny * 0.5,
        "xprobmax2": (pfss_ny - y0) / pfss_ny * 0.5,
        "xprobmin3": x0 / pfss_nx,
        "xprobmax3": (x0 + nx) / pfss_nx,
    }


def write_pfss_boundary(path, pfss_br, lon_degree, lat_degree, *, scale: float = 1.0):
    """Write `pfss_boundary.dat` in the RBSL/PFSS binary format."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    pfss_br = np.asarray(pfss_br, dtype=float)
    lon_degree = np.asarray(lon_degree, dtype=float)
    lat_degree = np.asarray(lat_degree, dtype=float)
    if pfss_br.shape != (lat_degree.size, lon_degree.size):
        raise ValueError("pfss_br shape must be (len(lat_degree), len(lon_degree))")

    coord_t = np.pi / 2.0 - np.radians(lat_degree)
    coord_p = np.pi + np.radians(lon_degree)
    with path.open("wb") as handle:
        np.array(pfss_br.shape[1], dtype=np.int32).tofile(handle)
        np.array(pfss_br.shape[0], dtype=np.int32).tofile(handle)
        np.asarray(coord_t, dtype=np.float64).tofile(handle)
        np.asarray(coord_p, dtype=np.float64).tofile(handle)
        np.asarray(pfss_br * scale, dtype=np.float64).tofile(handle)
    return path


def write_cutout(path, br, bt, bp, *, scale: float = 1.0, rotate: bool = True):
    """Write local vector boundary `cutout.dat` in the RBSL binary format."""

    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    br = np.asarray(br, dtype=float)
    bt = np.asarray(bt, dtype=float)
    bp = np.asarray(bp, dtype=float)
    if br.shape != bt.shape or br.shape != bp.shape:
        raise ValueError("br, bt, and bp must have the same shape")
    if rotate:
        br, bt, bp = np.rot90(br, k=-1), np.rot90(bt, k=-1), np.rot90(bp, k=-1)

    with path.open("wb") as handle:
        np.array(br.shape[1], dtype=np.int32).tofile(handle)
        np.array(br.shape[0], dtype=np.int32).tofile(handle)
        np.asarray(br * scale, dtype=np.float64).tofile(handle)
        np.asarray(bt * scale, dtype=np.float64).tofile(handle)
        np.asarray(bp * scale, dtype=np.float64).tofile(handle)
    return path


def _laplace2(field):
    field_ext = np.hstack([field, field[:, :1]])
    return laplace(field_ext)[:, :-1]


def _make_mesh(theta, phi):
    return np.tile(theta.reshape(-1, 1), (1, phi.size)), np.tile(phi.reshape(1, -1), (theta.size, 1))


def _preprocess_terms(br0, bt0, bp0, theta2, phi2, sinth, eb):
    ebm = 0.5 * (bt0**2 + bp0**2 - br0**2)
    cos_th = np.cos(theta2)
    sin_th = np.sin(theta2)
    cos_ph = np.cos(phi2)
    sin_ph = np.sin(phi2)

    b1 = bt0 * cos_th * cos_ph - bp0 * sin_ph
    b2 = bt0 * cos_th * sin_ph + bp0 * cos_ph
    b3 = bp0 * cos_th * cos_ph + bt0 * sin_ph
    b4 = bp0 * cos_th * sin_ph - bt0 * cos_ph

    term1a = np.sum(sinth * (ebm * sin_th * cos_ph - br0 * b1)) / eb
    term1b = np.sum(sinth * (ebm * sin_th * sin_ph - br0 * b2)) / eb
    term1c = np.sum(sinth * (ebm * cos_th + br0 * bt0 * sin_th)) / eb
    term2a = np.sum(sinth * br0 * b3) / eb
    term2b = np.sum(sinth * br0 * b4) / eb
    term2c = np.sum(sinth * sin_th * br0 * bp0) / eb
    return {
        "sin_th": sin_th,
        "cos_th": cos_th,
        "sin_ph": sin_ph,
        "cos_ph": cos_ph,
        "term1a": term1a,
        "term1b": term1b,
        "term1c": term1c,
        "term2a": term2a,
        "term2b": term2b,
        "term2c": term2c,
    }


def _safe_relative_change(value, previous):
    denominator = abs(value)
    if denominator == 0:
        return 0.0 if previous == 0 else np.inf
    return abs(value - previous) / denominator

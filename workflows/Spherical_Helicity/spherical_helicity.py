#!/usr/bin/env python3
"""Finite-volume relative magnetic helicity in a spherical wedge.

The implementation follows Moraitis et al. (2018, Solar Physics 293, 92):

* construct the potential reference field from a spherical Laplace problem
  whose Neumann data match the magnetic field on all six boundaries;
* construct vector potentials in the spherical DeVore gauge, ``A_r = 0``;
* integrate ``H_V``, ``H_J`` and ``H_PJ`` with the spherical volume element.

MPI-AMRVAC DAT snapshots are sampled onto a uniform cell-centred spherical
grid with the bundled ``simesh`` package.  Passing the original ``amrvac.par``
is mandatory so that stretched radial coordinates are reconstructed before
sampling.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
from scipy import sparse
from scipy.integrate import cumulative_trapezoid
from scipy.sparse.linalg import LinearOperator, cg


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
LOCAL_SIMESH_SRC = REPOSITORY_ROOT / "Packages" / "simesh" / "src"


@dataclass(frozen=True)
class SphericalGrid:
    """Uniform cell-centred ``(r, theta, phi)`` grid and exact cell geometry."""

    r_edges: np.ndarray
    theta_edges: np.ndarray
    phi_edges: np.ndarray

    @classmethod
    def from_bounds(
        cls,
        lower: Iterable[float],
        upper: Iterable[float],
        resolution: Iterable[int],
    ) -> "SphericalGrid":
        lower = np.asarray(tuple(lower), dtype=float)
        upper = np.asarray(tuple(upper), dtype=float)
        resolution = np.asarray(tuple(resolution), dtype=int)
        if lower.shape != (3,) or upper.shape != (3,) or resolution.shape != (3,):
            raise ValueError("lower, upper and resolution must each contain three entries")
        if np.any(upper <= lower) or np.any(resolution < 3):
            raise ValueError("grid bounds must increase and each resolution must be at least 3")
        if lower[1] <= 0.0 or upper[1] >= np.pi:
            raise ValueError("the present wedge solver excludes the theta=0 and theta=pi poles")
        return cls(
            np.linspace(lower[0], upper[0], resolution[0] + 1),
            np.linspace(lower[1], upper[1], resolution[1] + 1),
            np.linspace(lower[2], upper[2], resolution[2] + 1),
        )

    @property
    def r(self) -> np.ndarray:
        return 0.5 * (self.r_edges[1:] + self.r_edges[:-1])

    @property
    def theta(self) -> np.ndarray:
        return 0.5 * (self.theta_edges[1:] + self.theta_edges[:-1])

    @property
    def phi(self) -> np.ndarray:
        return 0.5 * (self.phi_edges[1:] + self.phi_edges[:-1])

    @property
    def shape(self) -> tuple[int, int, int]:
        return (len(self.r), len(self.theta), len(self.phi))

    @property
    def cell_volumes(self) -> np.ndarray:
        radial = np.diff(self.r_edges**3) / 3.0
        polar = np.cos(self.theta_edges[:-1]) - np.cos(self.theta_edges[1:])
        azimuth = np.diff(self.phi_edges)
        return radial[:, None, None] * polar[None, :, None] * azimuth[None, None, :]

    @property
    def analytic_volume(self) -> float:
        return float(
            (self.r_edges[-1] ** 3 - self.r_edges[0] ** 3)
            / 3.0
            * (np.cos(self.theta_edges[0]) - np.cos(self.theta_edges[-1]))
            * (self.phi_edges[-1] - self.phi_edges[0])
        )

    def cell_points(self) -> np.ndarray:
        rr, tt, pp = np.meshgrid(self.r, self.theta, self.phi, indexing="ij")
        return np.column_stack((rr.ravel(), tt.ravel(), pp.ravel()))


@dataclass
class BoundaryFluxes:
    """Corrected outward-normal boundary fields and their face areas."""

    q_rmin: np.ndarray
    q_rmax: np.ndarray
    q_tmin: np.ndarray
    q_tmax: np.ndarray
    q_pmin: np.ndarray
    q_pmax: np.ndarray
    area_rmin: np.ndarray
    area_rmax: np.ndarray
    area_tmin: np.ndarray
    area_tmax: np.ndarray
    area_pmin: np.ndarray
    area_pmax: np.ndarray
    raw_net_flux: float
    raw_unsigned_flux: float
    correction: float

    @property
    def relative_imbalance(self) -> float:
        return float(abs(self.raw_net_flux) / max(self.raw_unsigned_flux, np.finfo(float).tiny))

    @property
    def corrected_net_flux(self) -> float:
        pairs = (
            (self.q_rmin, self.area_rmin),
            (self.q_rmax, self.area_rmax),
            (self.q_tmin, self.area_tmin),
            (self.q_tmax, self.area_tmax),
            (self.q_pmin, self.area_pmin),
            (self.q_pmax, self.area_pmax),
        )
        return float(sum(np.sum(q * area) for q, area in pairs))


def _boundary_face_areas(grid: SphericalGrid) -> dict[str, np.ndarray]:
    dphi = np.diff(grid.phi_edges)
    dtheta = np.diff(grid.theta_edges)
    polar_solid_angle = (
        np.cos(grid.theta_edges[:-1]) - np.cos(grid.theta_edges[1:])
    )[:, None] * dphi[None, :]
    area_rmin = grid.r_edges[0] ** 2 * polar_solid_angle
    area_rmax = grid.r_edges[-1] ** 2 * polar_solid_angle

    shell = 0.5 * np.diff(grid.r_edges**2)
    area_tmin = shell[:, None] * np.sin(grid.theta_edges[0]) * dphi[None, :]
    area_tmax = shell[:, None] * np.sin(grid.theta_edges[-1]) * dphi[None, :]
    meridional = shell[:, None] * dtheta[None, :]
    return {
        "rmin": area_rmin,
        "rmax": area_rmax,
        "tmin": area_tmin,
        "tmax": area_tmax,
        "pmin": meridional,
        "pmax": meridional.copy(),
    }


def _extrapolate_cell_centres_to_face(
    field: np.ndarray,
    coordinates: np.ndarray,
    edge: float,
    *,
    axis: int,
    upper: bool,
) -> np.ndarray:
    """Quadratically extrapolate three cell centres to one physical face."""
    selection = np.arange(len(coordinates) - 3, len(coordinates)) if upper else np.arange(3)
    nodes = coordinates[selection]
    weights = np.ones(3, dtype=float)
    for i in range(3):
        for j in range(3):
            if i != j:
                weights[i] *= (edge - nodes[j]) / (nodes[i] - nodes[j])
    return sum(
        weights[i] * np.take(field, int(selection[i]), axis=axis)
        for i in range(3)
    )


def prepare_boundary_fluxes(
    grid: SphericalGrid,
    br: np.ndarray,
    btheta: np.ndarray,
    bphi: np.ndarray,
    face_components: dict[str, np.ndarray] | None = None,
) -> BoundaryFluxes:
    """Return minimally flux-balanced outward boundary data.

    Numerical data rarely satisfy the Neumann compatibility condition exactly.
    Following the finite-volume interpretation, a single uniform outward-field
    offset is removed from all boundary faces.  The offset and the original
    imbalance are retained as mandatory diagnostics.
    """
    for name, field in (("br", br), ("btheta", btheta), ("bphi", bphi)):
        if np.shape(field) != grid.shape:
            raise ValueError(f"{name} has shape {np.shape(field)}, expected {grid.shape}")
        if not np.all(np.isfinite(field)):
            raise ValueError(f"{name} contains non-finite values")

    areas = _boundary_face_areas(grid)
    if face_components is None:
        br_min = _extrapolate_cell_centres_to_face(
            br, grid.r, grid.r_edges[0], axis=0, upper=False
        )
        br_max = _extrapolate_cell_centres_to_face(
            br, grid.r, grid.r_edges[-1], axis=0, upper=True
        )
        bt_min = _extrapolate_cell_centres_to_face(
            btheta, grid.theta, grid.theta_edges[0], axis=1, upper=False
        )
        bt_max = _extrapolate_cell_centres_to_face(
            btheta, grid.theta, grid.theta_edges[-1], axis=1, upper=True
        )
        bp_min = _extrapolate_cell_centres_to_face(
            bphi, grid.phi, grid.phi_edges[0], axis=2, upper=False
        )
        bp_max = _extrapolate_cell_centres_to_face(
            bphi, grid.phi, grid.phi_edges[-1], axis=2, upper=True
        )
    else:
        required = {"rmin", "rmax", "tmin", "tmax", "pmin", "pmax"}
        missing = required.difference(face_components)
        if missing:
            raise ValueError(f"face_components is missing {sorted(missing)}")
        br_min, br_max = face_components["rmin"], face_components["rmax"]
        bt_min, bt_max = face_components["tmin"], face_components["tmax"]
        bp_min, bp_max = face_components["pmin"], face_components["pmax"]
    raw = {
        "rmin": -np.asarray(br_min, dtype=float),
        "rmax": np.asarray(br_max, dtype=float),
        "tmin": -np.asarray(bt_min, dtype=float),
        "tmax": np.asarray(bt_max, dtype=float),
        "pmin": -np.asarray(bp_min, dtype=float),
        "pmax": np.asarray(bp_max, dtype=float),
    }
    net = float(sum(np.sum(raw[key] * areas[key]) for key in raw))
    unsigned = float(sum(np.sum(np.abs(raw[key]) * areas[key]) for key in raw))
    total_area = float(sum(np.sum(areas[key]) for key in raw))
    correction = net / total_area
    corrected = {key: raw[key] - correction for key in raw}
    return BoundaryFluxes(
        q_rmin=corrected["rmin"],
        q_rmax=corrected["rmax"],
        q_tmin=corrected["tmin"],
        q_tmax=corrected["tmax"],
        q_pmin=corrected["pmin"],
        q_pmax=corrected["pmax"],
        area_rmin=areas["rmin"],
        area_rmax=areas["rmax"],
        area_tmin=areas["tmin"],
        area_tmax=areas["tmax"],
        area_pmin=areas["pmin"],
        area_pmax=areas["pmax"],
        raw_net_flux=net,
        raw_unsigned_flux=unsigned,
        correction=correction,
    )


def _add_pair_conductances(
    rows: list[np.ndarray],
    cols: list[np.ndarray],
    values: list[np.ndarray],
    left: np.ndarray,
    right: np.ndarray,
    conductance: np.ndarray,
) -> None:
    pair_shape = left.shape
    conductance = np.broadcast_to(conductance, pair_shape).ravel()
    left = left.ravel()
    right = right.ravel()
    rows.extend((left, right, left, right))
    cols.extend((left, right, right, left))
    values.extend((conductance, conductance, -conductance, -conductance))


def build_neumann_matrix(grid: SphericalGrid) -> sparse.csr_matrix:
    """Build the symmetric finite-volume matrix for ``-div(grad(phi))``."""
    nr, nt, np_ = grid.shape
    index = np.arange(nr * nt * np_, dtype=np.int64).reshape(grid.shape)
    rows: list[np.ndarray] = []
    cols: list[np.ndarray] = []
    values: list[np.ndarray] = []

    polar = np.cos(grid.theta_edges[:-1]) - np.cos(grid.theta_edges[1:])
    dphi = np.diff(grid.phi_edges)
    radial_area = (
        grid.r_edges[1:-1, None, None] ** 2
        * polar[None, :, None]
        * dphi[None, None, :]
    )
    radial_c = radial_area / np.diff(grid.r)[:, None, None]
    _add_pair_conductances(
        rows, cols, values, index[:-1], index[1:], radial_c
    )

    shell = 0.5 * np.diff(grid.r_edges**2)
    theta_area = (
        shell[:, None, None]
        * np.sin(grid.theta_edges[1:-1])[None, :, None]
        * dphi[None, None, :]
    )
    theta_distance = grid.r[:, None, None] * np.diff(grid.theta)[None, :, None]
    theta_c = theta_area / theta_distance
    _add_pair_conductances(
        rows, cols, values, index[:, :-1], index[:, 1:], theta_c
    )

    phi_area = shell[:, None, None] * np.diff(grid.theta_edges)[None, :, None]
    phi_distance = (
        grid.r[:, None, None]
        * np.sin(grid.theta)[None, :, None]
        * np.diff(grid.phi)[None, None, :]
    )
    phi_c = phi_area / phi_distance
    _add_pair_conductances(
        rows, cols, values, index[:, :, :-1], index[:, :, 1:], phi_c
    )

    row = np.concatenate(rows)
    col = np.concatenate(cols)
    data = np.concatenate([np.asarray(item, dtype=float).ravel() for item in values])
    matrix = sparse.coo_matrix((data, (row, col)), shape=(index.size, index.size))
    return matrix.tocsr()


def _neumann_rhs(grid: SphericalGrid, boundary: BoundaryFluxes) -> np.ndarray:
    rhs = np.zeros(grid.shape, dtype=float)
    # K stores -integral(div(grad(phi))) for the unknown internal faces.
    # Therefore K phi equals the prescribed outward flux on physical faces.
    rhs[0] += boundary.q_rmin * boundary.area_rmin
    rhs[-1] += boundary.q_rmax * boundary.area_rmax
    rhs[:, 0, :] += boundary.q_tmin * boundary.area_tmin
    rhs[:, -1, :] += boundary.q_tmax * boundary.area_tmax
    rhs[:, :, 0] += boundary.q_pmin * boundary.area_pmin
    rhs[:, :, -1] += boundary.q_pmax * boundary.area_pmax
    # Roundoff only: the physical compatibility correction was already made on
    # the boundary fields themselves and is reported to the caller.
    rhs -= np.sum(rhs) / rhs.size
    return rhs


def _solve_gauge_fixed_system(
    matrix: sparse.csr_matrix,
    rhs: np.ndarray,
    *,
    rtol: float,
    maxiter: int,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Solve a compatible singular SPD system after fixing its first value."""
    rhs = np.asarray(rhs, dtype=float).ravel().copy()
    rhs -= np.sum(rhs) / rhs.size
    reduced = matrix[1:, 1:]
    reduced_rhs = rhs[1:]
    diagonal = reduced.diagonal()
    if np.any(diagonal <= 0.0):
        raise RuntimeError("finite-volume matrix is not positive definite after gauge fixing")
    preconditioner = LinearOperator(
        reduced.shape, matvec=lambda value: value / diagonal, dtype=float
    )
    iterations = 0

    def count_iteration(_value: np.ndarray) -> None:
        nonlocal iterations
        iterations += 1

    try:
        solution, info = cg(
            reduced,
            reduced_rhs,
            rtol=rtol,
            atol=0.0,
            maxiter=maxiter,
            M=preconditioner,
            callback=count_iteration,
        )
    except TypeError:  # SciPy < 1.12 compatibility
        solution, info = cg(
            reduced,
            reduced_rhs,
            tol=rtol,
            maxiter=maxiter,
            M=preconditioner,
            callback=count_iteration,
        )
    if info != 0:
        raise RuntimeError(f"conjugate-gradient solve failed: info={info}")
    full = np.zeros(rhs.size, dtype=float)
    full[1:] = solution
    residual = matrix @ full - rhs
    return full, {
        "iterations": iterations,
        "relative_linear_residual": float(
            np.linalg.norm(residual) / max(np.linalg.norm(rhs), np.finfo(float).tiny)
        ),
    }


def solve_scalar_potential(
    grid: SphericalGrid,
    boundary: BoundaryFluxes,
    *,
    rtol: float = 1.0e-10,
    maxiter: int = 4000,
) -> tuple[np.ndarray, dict[str, float | int]]:
    """Solve the all-Neumann spherical Laplace problem with a fixed zero gauge."""
    matrix = build_neumann_matrix(grid)
    rhs = _neumann_rhs(grid, boundary).ravel()
    scalar, metrics = _solve_gauge_fixed_system(
        matrix, rhs, rtol=rtol, maxiter=maxiter
    )
    scalar -= np.average(scalar, weights=grid.cell_volumes.ravel())
    metrics["corrected_boundary_flux"] = boundary.corrected_net_flux
    return scalar.reshape(grid.shape), metrics


def _all_face_areas(
    grid: SphericalGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    polar = np.cos(grid.theta_edges[:-1]) - np.cos(grid.theta_edges[1:])
    dphi = np.diff(grid.phi_edges)
    area_r = (
        grid.r_edges[:, None, None] ** 2
        * polar[None, :, None]
        * dphi[None, None, :]
    )
    shell = 0.5 * np.diff(grid.r_edges**2)
    area_theta = (
        shell[:, None, None]
        * np.sin(grid.theta_edges)[None, :, None]
        * dphi[None, None, :]
    )
    area_phi = (
        shell[:, None, None]
        * np.diff(grid.theta_edges)[None, :, None]
        * np.ones((1, 1, len(grid.phi_edges)))
    )
    return area_r, area_theta, area_phi


def _cell_fields_to_faces(
    grid: SphericalGrid,
    field: tuple[np.ndarray, np.ndarray, np.ndarray],
    boundary: BoundaryFluxes,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    br, bt, bp = field
    fr = np.empty((grid.shape[0] + 1, grid.shape[1], grid.shape[2]), dtype=float)
    ft = np.empty((grid.shape[0], grid.shape[1] + 1, grid.shape[2]), dtype=float)
    fp = np.empty((grid.shape[0], grid.shape[1], grid.shape[2] + 1), dtype=float)
    fr[0], fr[-1] = -boundary.q_rmin, boundary.q_rmax
    ft[:, 0, :], ft[:, -1, :] = -boundary.q_tmin, boundary.q_tmax
    fp[:, :, 0], fp[:, :, -1] = -boundary.q_pmin, boundary.q_pmax
    fr[1:-1] = 0.5 * (br[:-1] + br[1:])
    ft[:, 1:-1] = 0.5 * (bt[:, :-1] + bt[:, 1:])
    fp[:, :, 1:-1] = 0.5 * (bp[:, :, :-1] + bp[:, :, 1:])
    return fr, ft, fp


def _faces_to_cell_fields(
    faces: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    fr, ft, fp = faces
    return (
        0.5 * (fr[:-1] + fr[1:]),
        0.5 * (ft[:, :-1] + ft[:, 1:]),
        0.5 * (fp[:, :, :-1] + fp[:, :, 1:]),
    )


def finite_volume_flux_residual(
    grid: SphericalGrid,
    faces: tuple[np.ndarray, np.ndarray, np.ndarray],
) -> np.ndarray:
    """Net outward magnetic flux of every uniform spherical cell."""
    fr, ft, fp = faces
    ar, at, ap = _all_face_areas(grid)
    return (
        fr[1:] * ar[1:] - fr[:-1] * ar[:-1]
        + ft[:, 1:] * at[:, 1:] - ft[:, :-1] * at[:, :-1]
        + fp[:, :, 1:] * ap[:, :, 1:] - fp[:, :, :-1] * ap[:, :, :-1]
    )


def project_solenoidal(
    grid: SphericalGrid,
    magnetic: tuple[np.ndarray, np.ndarray, np.ndarray],
    boundary: BoundaryFluxes,
    *,
    rtol: float,
    maxiter: int,
) -> tuple[tuple[np.ndarray, np.ndarray, np.ndarray], dict[str, float | int]]:
    """Finite-volume Hodge projection preserving all corrected boundary fluxes."""
    faces = _cell_fields_to_faces(grid, magnetic, boundary)
    raw_residual = finite_volume_flux_residual(grid, faces)
    matrix = build_neumann_matrix(grid)
    correction, metrics = _solve_gauge_fixed_system(
        matrix, -raw_residual.ravel(), rtol=rtol, maxiter=maxiter
    )
    correction = correction.reshape(grid.shape)

    fr, ft, fp = (item.copy() for item in faces)
    fr[1:-1] -= (
        correction[1:] - correction[:-1]
    ) / np.diff(grid.r)[:, None, None]
    ft[:, 1:-1] -= (
        correction[:, 1:] - correction[:, :-1]
    ) / (grid.r[:, None, None] * np.diff(grid.theta)[None, :, None])
    fp[:, :, 1:-1] -= (
        correction[:, :, 1:] - correction[:, :, :-1]
    ) / (
        grid.r[:, None, None]
        * np.sin(grid.theta)[None, :, None]
        * np.diff(grid.phi)[None, None, :]
    )
    cleaned_faces = (fr, ft, fp)
    cleaned = _faces_to_cell_fields(cleaned_faces)
    cleaned_residual = finite_volume_flux_residual(grid, cleaned_faces)
    volumes = grid.cell_volumes
    raw_strength = _field_strength_l2(magnetic, volumes)
    correction_strength = _field_strength_l2(
        tuple(raw - clean for raw, clean in zip(magnetic, cleaned)), volumes
    )
    metrics.update(
        {
            "raw_flux_residual_l2": float(np.linalg.norm(raw_residual.ravel())),
            "cleaned_flux_residual_l2": float(np.linalg.norm(cleaned_residual.ravel())),
            "relative_field_correction_l2": float(
                correction_strength / max(raw_strength, np.finfo(float).tiny)
            ),
        }
    )
    return cleaned, metrics


def scalar_gradient(
    scalar: np.ndarray,
    grid: SphericalGrid,
    boundary: BoundaryFluxes | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Finite-volume spherical gradient reconstructed from face derivatives."""
    if boundary is None:
        raise ValueError("finite-volume scalar_gradient requires Neumann boundary data")
    nr, nt, np_ = grid.shape
    fr = np.empty((nr + 1, nt, np_), dtype=float)
    ft = np.empty((nr, nt + 1, np_), dtype=float)
    fp = np.empty((nr, nt, np_ + 1), dtype=float)
    fr[0], fr[-1] = -boundary.q_rmin, boundary.q_rmax
    ft[:, 0, :], ft[:, -1, :] = -boundary.q_tmin, boundary.q_tmax
    fp[:, :, 0], fp[:, :, -1] = -boundary.q_pmin, boundary.q_pmax
    fr[1:-1] = (scalar[1:] - scalar[:-1]) / np.diff(grid.r)[:, None, None]
    ft[:, 1:-1] = (scalar[:, 1:] - scalar[:, :-1]) / (
        grid.r[:, None, None] * np.diff(grid.theta)[None, :, None]
    )
    fp[:, :, 1:-1] = (scalar[:, :, 1:] - scalar[:, :, :-1]) / (
        grid.r[:, None, None]
        * np.sin(grid.theta)[None, :, None]
        * np.diff(grid.phi)[None, None, :]
    )
    return _faces_to_cell_fields((fr, ft, fp))


def divergence_spherical(
    br: np.ndarray,
    btheta: np.ndarray,
    bphi: np.ndarray,
    grid: SphericalGrid,
) -> np.ndarray:
    r = grid.r[:, None, None]
    sin_theta = np.sin(grid.theta)[None, :, None]
    return (
        np.gradient(r**2 * br, grid.r, axis=0, edge_order=2) / r**2
        + np.gradient(sin_theta * btheta, grid.theta, axis=1, edge_order=2)
        / (r * sin_theta)
        + np.gradient(bphi, grid.phi, axis=2, edge_order=2) / (r * sin_theta)
    )


def curl_spherical_de_vore(
    atheta: np.ndarray,
    aphi: np.ndarray,
    grid: SphericalGrid,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Curl of ``(0, A_theta, A_phi)`` on the cell-centred spherical grid."""
    r = grid.r[:, None, None]
    sin_theta = np.sin(grid.theta)[None, :, None]
    br = (
        np.gradient(sin_theta * aphi, grid.theta, axis=1, edge_order=2)
        - np.gradient(atheta, grid.phi, axis=2, edge_order=2)
    ) / (r * sin_theta)
    btheta = -np.gradient(r * aphi, grid.r, axis=0, edge_order=2) / r
    bphi = np.gradient(r * atheta, grid.r, axis=0, edge_order=2) / r
    return br, btheta, bphi


def _integral_from_reference(
    values: np.ndarray, coordinates: np.ndarray, reference: str
) -> np.ndarray:
    if reference == "bottom":
        return cumulative_trapezoid(values, coordinates, axis=0, initial=0.0)
    if reference == "top":
        reversed_result = cumulative_trapezoid(
            values[::-1], coordinates[::-1], axis=0, initial=0.0
        )
        return reversed_result[::-1]
    raise ValueError("reference must be 'bottom' or 'top'")


def de_vore_vector_potential(
    br: np.ndarray,
    btheta: np.ndarray,
    bphi: np.ndarray,
    grid: SphericalGrid,
    *,
    reference: str = "bottom",
    c: float = 0.5,
) -> tuple[np.ndarray, np.ndarray]:
    """Spherical DeVore-simple vector potential from Moraitis et al. (2018)."""
    if not 0.0 <= c <= 1.0:
        raise ValueError("the DVS gauge parameter c must lie in [0, 1]")
    index = 0 if reference == "bottom" else -1
    r0 = float(grid.r[index])
    surface_br = br[index]
    sin_theta_2d = np.sin(grid.theta)[:, None]

    theta_integral = cumulative_trapezoid(
        sin_theta_2d * surface_br,
        grid.theta,
        axis=0,
        initial=0.0,
    )
    phi_integral = cumulative_trapezoid(
        surface_br,
        grid.phi,
        axis=1,
        initial=0.0,
    )
    surface_aphi = c * r0 * theta_integral / sin_theta_2d
    surface_atheta = -(1.0 - c) * r0 * sin_theta_2d * phi_integral

    radial_bphi = _integral_from_reference(
        grid.r[:, None, None] * bphi, grid.r, reference
    )
    radial_btheta = _integral_from_reference(
        grid.r[:, None, None] * btheta, grid.r, reference
    )
    radius = grid.r[:, None, None]
    atheta = (r0 * surface_atheta[None, :, :] + radial_bphi) / radius
    aphi = (r0 * surface_aphi[None, :, :] - radial_btheta) / radius
    return atheta, aphi


def weighted_relative_l2(
    calculated: tuple[np.ndarray, np.ndarray, np.ndarray],
    expected: tuple[np.ndarray, np.ndarray, np.ndarray],
    weights: np.ndarray,
    *,
    interior: bool = True,
) -> float:
    selection = np.s_[1:-1, 1:-1, 1:-1] if interior else np.s_[:, :, :]
    numerator = sum((a[selection] - b[selection]) ** 2 for a, b in zip(calculated, expected))
    denominator = sum(b[selection] ** 2 for b in expected)
    return float(
        np.sqrt(
            np.sum(numerator * weights[selection])
            / max(np.sum(denominator * weights[selection]), np.finfo(float).tiny)
        )
    )


def calculate_helicities(
    grid: SphericalGrid,
    magnetic: tuple[np.ndarray, np.ndarray, np.ndarray],
    potential: tuple[np.ndarray, np.ndarray, np.ndarray],
    *,
    reference: str,
    dvs_c: float = 0.5,
) -> tuple[dict[str, float], dict[str, np.ndarray]]:
    br, bt, bp = magnetic
    pbr, pbt, pbp = potential
    at, ap = de_vore_vector_potential(br, bt, bp, grid, reference=reference, c=dvs_c)
    pat, pap = de_vore_vector_potential(
        pbr, pbt, pbp, grid, reference=reference, c=dvs_c
    )
    jr, jt, jp = br - pbr, bt - pbt, bp - pbp
    volumes = grid.cell_volumes
    hv_density = (at + pat) * jt + (ap + pap) * jp
    hj_density = (at - pat) * jt + (ap - pap) * jp
    hpj_density = 2.0 * (pat * jt + pap * jp)
    hv = float(np.sum(hv_density * volumes))
    hj = float(np.sum(hj_density * volumes))
    hpj = float(np.sum(hpj_density * volumes))

    curl_a = curl_spherical_de_vore(at, ap, grid)
    curl_ap = curl_spherical_de_vore(pat, pap, grid)
    metrics = {
        "H_V": hv,
        "H_J": hj,
        "H_PJ": hpj,
        "decomposition_closure": float(
            abs(hv - hj - hpj) / max(abs(hv), abs(hj) + abs(hpj), np.finfo(float).tiny)
        ),
        "curl_A_relative_l2_interior": weighted_relative_l2(
            curl_a, magnetic, volumes, interior=True
        ),
        "curl_Ap_relative_l2_interior": weighted_relative_l2(
            curl_ap, potential, volumes, interior=True
        ),
    }
    arrays = {
        "A_theta": at,
        "A_phi": ap,
        "Ap_theta": pat,
        "Ap_phi": pap,
        "hV_density": hv_density,
        "hJ_density": hj_density,
        "hPJ_density": hpj_density,
    }
    return metrics, arrays


def _resolve_magnetic_indices(dataset) -> list[int]:
    names = list(dataset.metadata["w_names"])
    missing = [name for name in ("b1", "b2", "b3") if name not in names]
    if missing:
        raise KeyError(f"DAT is missing spherical magnetic components: {missing}")
    return [names.index(name) for name in ("b1", "b2", "b3")]


def load_uniform_spherical_dat(
    dat_path: Path,
    grid_config: Path,
    resolution: tuple[int, int, int],
    *,
    interpolation: str = "linear",
    batch_size: int = 200_000,
) -> tuple[
    SphericalGrid,
    tuple[np.ndarray, np.ndarray, np.ndarray],
    dict[str, np.ndarray],
    dict,
]:
    """Sample a stretched spherical DAT on a uniform spherical grid."""
    if LOCAL_SIMESH_SRC.is_dir():
        sys.path.insert(0, str(LOCAL_SIMESH_SRC))
    try:
        from simesh.amrvac import open_dataset
    except (ImportError, ModuleNotFoundError) as exc:
        raise ModuleNotFoundError(
            "simesh is required; install the bundled Packages/simesh package"
        ) from exc

    dataset = open_dataset(str(dat_path), ghost_width=0, grid_config=str(grid_config))
    if int(dataset.ndim) != 3 or "spherical" not in str(dataset.geometry).lower():
        raise ValueError("input must be a three-dimensional spherical MPI-AMRVAC DAT")
    if not dataset.uses_stretched_coordinates:
        raise ValueError(
            "grid_config did not enable stretched coordinates; verify that it is the DAT's original amrvac.par"
        )
    grid = SphericalGrid.from_bounds(
        dataset.physical_domain[0], dataset.physical_domain[1], resolution
    )
    indices = _resolve_magnetic_indices(dataset)
    dataset.load_data(field_indices=indices)
    points = grid.cell_points()
    sampled = np.empty((len(points), 3), dtype=float)
    for start in range(0, len(points), batch_size):
        stop = min(start + batch_size, len(points))
        sampled[start:stop] = dataset.sample_points(
            points[start:stop], field_indices=indices, interpolation=interpolation
        )
    if not np.all(np.isfinite(sampled)):
        count = int(np.count_nonzero(~np.isfinite(sampled)))
        raise RuntimeError(f"simesh returned {count} non-finite magnetic samples")
    fields = tuple(
        np.ascontiguousarray(sampled[:, component].reshape(grid.shape))
        for component in range(3)
    )

    def sample_face(
        first: np.ndarray,
        second: np.ndarray,
        third: np.ndarray,
        component: int,
    ) -> np.ndarray:
        aa, bb, cc = np.meshgrid(first, second, third, indexing="ij")
        face_points = np.column_stack((aa.ravel(), bb.ravel(), cc.ravel()))
        values = dataset.sample_points(
            face_points,
            field_indices=[indices[component]],
            interpolation=interpolation,
        )[:, 0]
        if not np.all(np.isfinite(values)):
            raise RuntimeError("simesh returned non-finite samples on a physical boundary")
        return values.reshape(aa.shape).squeeze(axis=component)

    face_components = {
        "rmin": sample_face(
            np.array([grid.r_edges[0]]), grid.theta, grid.phi, 0
        ),
        "rmax": sample_face(
            np.array([grid.r_edges[-1]]), grid.theta, grid.phi, 0
        ),
        "tmin": sample_face(
            grid.r, np.array([grid.theta_edges[0]]), grid.phi, 1
        ),
        "tmax": sample_face(
            grid.r, np.array([grid.theta_edges[-1]]), grid.phi, 1
        ),
        "pmin": sample_face(
            grid.r, grid.theta, np.array([grid.phi_edges[0]]), 2
        ),
        "pmax": sample_face(
            grid.r, grid.theta, np.array([grid.phi_edges[-1]]), 2
        ),
    }
    metadata = {
        "geometry": str(dataset.geometry),
        "time_code": float(dataset.metadata.get("time", 0.0)),
        "iteration": int(dataset.metadata.get("it", 0)),
        "nleafs": int(dataset.nleafs),
        "levmax": int(dataset.levmax),
        "domain_nx": [int(value) for value in dataset.domain_nx],
        "block_nx": [int(value) for value in dataset.block_nx],
        "coordinate_source": dataset.metadata.get("coordinate_source"),
        "uses_stretched_coordinates": bool(dataset.uses_stretched_coordinates),
        "interpolation": interpolation,
    }
    return grid, fields, face_components, metadata


def _field_strength_l2(
    field: tuple[np.ndarray, np.ndarray, np.ndarray], volumes: np.ndarray
) -> float:
    return float(np.sqrt(np.sum(sum(component**2 for component in field) * volumes)))


def run_calculation(args: argparse.Namespace) -> tuple[dict, dict[str, np.ndarray]]:
    grid, magnetic, face_components, dat_metadata = load_uniform_spherical_dat(
        args.dat,
        args.grid_config,
        args.resolution,
        interpolation=args.interpolation,
        batch_size=args.batch_size,
    )
    raw_magnetic = magnetic
    raw_br, raw_btheta, raw_bphi = raw_magnetic
    boundary = prepare_boundary_fluxes(
        grid, raw_br, raw_btheta, raw_bphi, face_components=face_components
    )
    if boundary.relative_imbalance > args.max_relative_flux_imbalance:
        raise RuntimeError(
            "boundary flux imbalance is too large for a gauge-invariant result: "
            f"{boundary.relative_imbalance:.6e} > {args.max_relative_flux_imbalance:.6e}"
        )
    magnetic, cleaning_metrics = project_solenoidal(
        grid,
        raw_magnetic,
        boundary,
        rtol=args.solver_rtol,
        maxiter=args.solver_maxiter,
    )
    br, btheta, bphi = magnetic
    scalar, solve_metrics = solve_scalar_potential(
        grid, boundary, rtol=args.solver_rtol, maxiter=args.solver_maxiter
    )
    potential = scalar_gradient(scalar, grid, boundary)
    volumes = grid.cell_volumes

    divergence = divergence_spherical(br, btheta, bphi, grid)
    characteristic_length = grid.analytic_volume ** (1.0 / 3.0)
    divergence_fraction = float(
        characteristic_length
        * np.sqrt(np.sum(divergence**2 * volumes))
        / max(_field_strength_l2(magnetic, volumes), np.finfo(float).tiny)
    )

    references: dict[str, dict[str, float]] = {}
    saved_arrays: dict[str, np.ndarray] = {}
    for reference in args.references:
        values, arrays = calculate_helicities(
            grid, magnetic, potential, reference=reference, dvs_c=args.dvs_c
        )
        references[reference] = values
        for name, array in arrays.items():
            saved_arrays[f"{reference}_{name}"] = array

    gauge_spread: dict[str, float] = {}
    if len(references) > 1:
        for quantity in ("H_V", "H_J", "H_PJ"):
            values = np.array([result[quantity] for result in references.values()])
            gauge_spread[quantity] = float(
                (np.max(values) - np.min(values))
                / max(np.max(np.abs(values)), np.finfo(float).tiny)
            )

    inner_unsigned_flux = float(
        np.sum(np.abs(boundary.q_rmin) * boundary.area_rmin)
    )
    for values in references.values():
        flux_squared = max(inner_unsigned_flux**2, np.finfo(float).tiny)
        values["H_V_over_inner_unsigned_flux_squared"] = values["H_V"] / flux_squared
        values["H_J_over_inner_unsigned_flux_squared"] = values["H_J"] / flux_squared
        values["H_PJ_over_inner_unsigned_flux_squared"] = values["H_PJ"] / flux_squared

    preferred_reference = min(
        references,
        key=lambda name: references[name]["curl_A_relative_l2_interior"]
        + references[name]["curl_Ap_relative_l2_interior"],
    )
    preferred_curl_error = max(
        references[preferred_reference]["curl_A_relative_l2_interior"],
        references[preferred_reference]["curl_Ap_relative_l2_interior"],
    )
    warnings: list[str] = []
    if cleaning_metrics["relative_field_correction_l2"] > 0.05:
        warnings.append(
            "uniform resampling required a solenoidal field correction larger than 5%; perform a resolution study"
        )
    if preferred_curl_error > 0.10:
        warnings.append(
            "the best DeVore curl reconstruction error exceeds 10%; scalar helicities are diagnostic only"
        )
    quality_status = "validated" if not warnings else "diagnostic_only"
    total_energy = float(
        np.sum((br**2 + btheta**2 + bphi**2) * volumes) / (8.0 * np.pi)
    )
    potential_energy = float(
        np.sum(sum(component**2 for component in potential) * volumes) / (8.0 * np.pi)
    )
    result = {
        "method": {
            "reference": "Moraitis et al. 2018, Solar Physics 293:92",
            "geometry": "finite spherical wedge",
            "potential_field": "cell-centred finite-volume all-Neumann Laplace solve",
            "vector_potential": f"DeVore simple gauge, A_r=0, c={args.dvs_c}",
            "helicity_definition": "Finn-Antonsen relative helicity",
        },
        "input": {
            "dat": str(args.dat.resolve()),
            "grid_config": str(args.grid_config.resolve()),
            **dat_metadata,
        },
        "uniform_grid": {
            "resolution": list(grid.shape),
            "r_bounds": [float(grid.r_edges[0]), float(grid.r_edges[-1])],
            "theta_bounds": [float(grid.theta_edges[0]), float(grid.theta_edges[-1])],
            "phi_bounds": [float(grid.phi_edges[0]), float(grid.phi_edges[-1])],
            "volume_code": grid.analytic_volume,
        },
        "boundary_flux": {
            "raw_net": boundary.raw_net_flux,
            "raw_unsigned": boundary.raw_unsigned_flux,
            "relative_imbalance": boundary.relative_imbalance,
            "uniform_outward_field_correction": boundary.correction,
            "corrected_net": boundary.corrected_net_flux,
        },
        "potential_solver": solve_metrics,
        "solenoidal_projection": cleaning_metrics,
        "quality": {
            "status": quality_status,
            "warnings": warnings,
            "preferred_reference": preferred_reference,
            "cell_centered_differential_divergence_l2": divergence_fraction,
            "gauge_spread": gauge_spread,
        },
        "integrals_code_units": {
            "total_magnetic_energy_over_8pi": total_energy,
            "potential_magnetic_energy_over_8pi": potential_energy,
            "free_magnetic_energy_over_8pi": total_energy - potential_energy,
            "inner_unsigned_flux": inner_unsigned_flux,
            "references": references,
        },
    }
    if args.unit_length_cm is not None and args.unit_b_gauss is not None:
        helicity_scale = args.unit_b_gauss**2 * args.unit_length_cm**4
        flux_scale = args.unit_b_gauss * args.unit_length_cm**2
        result["physical_units"] = {
            "unit_length_cm": args.unit_length_cm,
            "unit_b_gauss": args.unit_b_gauss,
            "inner_unsigned_flux_mx": inner_unsigned_flux * flux_scale,
            "references": {
                reference: {
                    f"{quantity}_mx2": value * helicity_scale
                    for quantity, value in values.items()
                    if quantity in {"H_V", "H_J", "H_PJ"}
                }
                for reference, values in references.items()
            },
        }

    saved_arrays.update(
        {
            "r": grid.r,
            "theta": grid.theta,
            "phi": grid.phi,
            "cell_volumes": volumes,
            "B_r": br,
            "B_theta": btheta,
            "B_phi": bphi,
            "B_raw_r": raw_br,
            "B_raw_theta": raw_btheta,
            "B_raw_phi": raw_bphi,
            "potential_scalar": scalar,
            "Bp_r": potential[0],
            "Bp_theta": potential[1],
            "Bp_phi": potential[2],
            "divB": divergence,
        }
    )
    return result, saved_arrays


def parse_resolution(text: str) -> tuple[int, int, int]:
    try:
        values = tuple(int(item) for item in text.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError("resolution must be nr,ntheta,nphi") from exc
    if len(values) != 3 or min(values) < 3:
        raise argparse.ArgumentTypeError("resolution must contain three integers >= 3")
    return values


def parse_references(text: str) -> tuple[str, ...]:
    values = tuple(item.strip().lower() for item in text.split(",") if item.strip())
    if not values or any(item not in {"bottom", "top"} for item in values):
        raise argparse.ArgumentTypeError("references must be bottom, top, or bottom,top")
    return tuple(dict.fromkeys(values))


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Compute relative magnetic helicity in a stretched spherical AMRVAC DAT."
    )
    parser.add_argument("dat", type=Path)
    parser.add_argument("--grid-config", type=Path, required=True, help="original amrvac.par")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--resolution", type=parse_resolution, default=(64, 32, 64))
    parser.add_argument("--interpolation", choices=("linear", "nearest"), default="linear")
    parser.add_argument("--batch-size", type=int, default=200_000)
    parser.add_argument("--references", type=parse_references, default=("bottom", "top"))
    parser.add_argument("--dvs-c", type=float, default=0.5)
    parser.add_argument("--solver-rtol", type=float, default=1.0e-10)
    parser.add_argument("--solver-maxiter", type=int, default=4000)
    parser.add_argument("--max-relative-flux-imbalance", type=float, default=5.0e-3)
    parser.add_argument("--unit-length-cm", type=float)
    parser.add_argument("--unit-b-gauss", type=float)
    parser.add_argument(
        "--no-save-fields",
        action="store_true",
        help="write only the JSON scalar result, not the compressed field archive",
    )
    return parser


def main() -> None:
    args = build_parser().parse_args()
    if args.batch_size < 1 or args.solver_rtol <= 0.0 or args.solver_maxiter < 1:
        raise ValueError("batch size and solver controls must be positive")
    if (args.unit_length_cm is None) != (args.unit_b_gauss is None):
        raise ValueError("physical conversion requires both --unit-length-cm and --unit-b-gauss")
    args.output.mkdir(parents=True, exist_ok=True)
    result, arrays = run_calculation(args)
    result_path = args.output / "spherical_helicity.json"
    result_path.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n")
    if not args.no_save_fields:
        np.savez_compressed(args.output / "spherical_helicity_fields.npz", **arrays)
    print(json.dumps(result, indent=2, ensure_ascii=False))
    print(f"result: {result_path}")


if __name__ == "__main__":
    main()

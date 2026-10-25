from __future__ import annotations

import importlib.util
from pathlib import Path
import sys

import numpy as np


MODULE_PATH = (
    Path(__file__).resolve().parents[1]
    / "workflows"
    / "Spherical_Helicity"
    / "spherical_helicity.py"
)
SPEC = importlib.util.spec_from_file_location("spherical_helicity", MODULE_PATH)
assert SPEC is not None and SPEC.loader is not None
sh = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = sh
SPEC.loader.exec_module(sh)


def uniform_cartesian_field(grid, vector=(0.4, -0.2, 0.8)):
    rr, theta, phi = np.meshgrid(grid.r, grid.theta, grid.phi, indexing="ij")
    bx, by, bz = vector
    br = (
        bx * np.sin(theta) * np.cos(phi)
        + by * np.sin(theta) * np.sin(phi)
        + bz * np.cos(theta)
    )
    btheta = (
        bx * np.cos(theta) * np.cos(phi)
        + by * np.cos(theta) * np.sin(phi)
        - bz * np.sin(theta)
    )
    bphi = -bx * np.sin(phi) + by * np.cos(phi)
    scalar = rr * br
    return scalar, (br, btheta, bphi)


def test_exact_cell_volume_matches_wedge_formula():
    grid = sh.SphericalGrid.from_bounds(
        (1.0, 0.7, 2.1), (1.5, 1.3, 2.8), (9, 8, 7)
    )
    assert np.isclose(np.sum(grid.cell_volumes), grid.analytic_volume, rtol=2.0e-15)


def test_neumann_potential_recovers_harmonic_linear_scalar():
    grid = sh.SphericalGrid.from_bounds(
        (1.0, 0.75, 2.2), (1.4, 1.3, 2.9), (14, 13, 12)
    )
    scalar_true, magnetic = uniform_cartesian_field(grid)
    boundary = sh.prepare_boundary_fluxes(grid, *magnetic)
    scalar, metrics = sh.solve_scalar_potential(grid, boundary, rtol=1.0e-11)
    scalar -= np.average(scalar, weights=grid.cell_volumes)
    scalar_true -= np.average(scalar_true, weights=grid.cell_volumes)
    relative = np.sqrt(
        np.sum((scalar - scalar_true) ** 2 * grid.cell_volumes)
        / np.sum(scalar_true**2 * grid.cell_volumes)
    )
    assert boundary.relative_imbalance < 2.0e-3
    assert metrics["relative_linear_residual"] < 1.0e-9
    assert relative < 0.04


def test_de_vore_curl_reconstructs_uniform_field_interior():
    grid = sh.SphericalGrid.from_bounds(
        (1.0, 0.8, 2.3), (1.45, 1.25, 2.85), (24, 22, 24)
    )
    _, magnetic = uniform_cartesian_field(grid)
    for reference in ("bottom", "top"):
        atheta, aphi = sh.de_vore_vector_potential(
            *magnetic, grid, reference=reference
        )
        reconstructed = sh.curl_spherical_de_vore(atheta, aphi, grid)
        error = sh.weighted_relative_l2(
            reconstructed, magnetic, grid.cell_volumes, interior=True
        )
        assert error < 0.035


def test_relative_and_current_helicity_vanish_when_field_is_reference():
    grid = sh.SphericalGrid.from_bounds(
        (1.0, 0.8, 2.3), (1.4, 1.2, 2.8), (12, 11, 10)
    )
    _, magnetic = uniform_cartesian_field(grid)
    values, _ = sh.calculate_helicities(
        grid, magnetic, magnetic, reference="bottom"
    )
    assert values["H_V"] == 0.0
    assert values["H_J"] == 0.0
    assert values["H_PJ"] == 0.0
    assert values["decomposition_closure"] == 0.0


def test_finite_volume_projection_removes_flux_divergence():
    grid = sh.SphericalGrid.from_bounds(
        (1.0, 0.8, 2.3), (1.4, 1.2, 2.8), (10, 9, 8)
    )
    generator = np.random.default_rng(42)
    magnetic = tuple(generator.normal(size=grid.shape) for _ in range(3))
    boundary = sh.prepare_boundary_fluxes(grid, *magnetic)
    _, metrics = sh.project_solenoidal(
        grid, magnetic, boundary, rtol=1.0e-10, maxiter=2000
    )
    assert metrics["relative_linear_residual"] < 1.0e-9
    assert metrics["cleaned_flux_residual_l2"] < 1.0e-8 * metrics["raw_flux_residual_l2"]

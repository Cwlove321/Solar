#!/usr/bin/env python3
"""Focused numerical tests for spherical_slice_analysis.py."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from types import SimpleNamespace

import numpy as np


SCRIPT = Path(__file__).with_name("spherical_slice_analysis.py")
SPEC = importlib.util.spec_from_file_location("spherical_slice_analysis", SCRIPT)
assert SPEC is not None and SPEC.loader is not None
TOOL = importlib.util.module_from_spec(SPEC)
sys.modules[SPEC.name] = TOOL
SPEC.loader.exec_module(TOOL)


def test_fortran_parameter_parser(tmp_path: Path) -> None:
    par = tmp_path / "amrvac.par"
    par.write_text(
        """
&mhd_list
 mhd_eta = 2.5d-6
 mhd_semirelativistic = .true.
/
&rc_list
 coolcurve = 'JCcorona'
/
""",
        encoding="utf-8",
    )
    parsed = TOOL.parse_amrvac_par(par)
    assert parsed["mhd_eta"] == 2.5e-6
    assert parsed["mhd_semirelativistic"] is True
    assert parsed["coolcurve"] == "JCcorona"
    assert parsed["mhd_list.mhd_eta"] == 2.5e-6


def test_spherical_cartesian_vector_round_trip_basis() -> None:
    theta = np.array([[0.7, 1.2]])
    phi = np.array([[2.1, 3.0]])
    radial = np.ones_like(theta)
    zero = np.zeros_like(theta)
    cart = TOOL.spherical_vectors_to_cartesian(radial, zero, zero, theta, phi)
    expected = TOOL.spherical_to_cartesian_points(radial, theta, phi)
    assert np.allclose(cart, expected)
    assert np.allclose(np.linalg.norm(cart, axis=-1), 1.0)


def test_cartesian_plane_derivatives() -> None:
    u = np.linspace(-0.2, 0.2, 9)
    v = np.linspace(-0.3, 0.3, 11)
    w = np.linspace(-0.1, 0.1, 5)
    ww, vv, uu = np.meshgrid(w, v, u, indexing="ij")
    # B=(-v,u,0) in an identity plane basis: curl(B)_w=2, div(B)=0.
    vector = np.stack((-vv, uu, np.zeros_like(uu)), axis=-1)
    dummy = TOOL.SampleGrid(
        mode="plane",
        points_spherical=np.empty(uu.shape + (3,)),
        u=u,
        v=v,
        w=w,
        r=np.ones_like(uu),
        theta=np.ones_like(uu),
        phi=np.ones_like(uu),
        axis_map=None,
        basis_u=np.array([1.0, 0.0, 0.0]),
        basis_v=np.array([0.0, 1.0, 0.0]),
        basis_w=np.array([0.0, 0.0, 1.0]),
        metadata={},
    )
    curl, divergence = TOOL.cartesian_plane_vector_derivatives(
        vector, dummy, rsun_mm=1.0
    )
    assert np.allclose(curl[..., 0], 0.0, atol=1.0e-12)
    assert np.allclose(curl[..., 1], 0.0, atol=1.0e-12)
    assert np.allclose(curl[..., 2], 2.0, atol=1.0e-12)
    assert np.allclose(divergence, 0.0, atol=1.0e-12)


def _spherical_grid(mode: str) -> TOOL.SampleGrid:
    r_values = np.linspace(1.02, 1.20, 7)
    theta_values = np.linspace(0.8, 1.2, 9)
    phi_values = np.linspace(2.0, 2.2, 5)
    if mode == "phi":
        phi, radius, theta = np.meshgrid(
            phi_values, r_values, theta_values, indexing="ij"
        )
        u, v, w = theta_values, r_values, phi_values
        axis_map = {"phi": 0, "r": 1, "theta": 2}
    elif mode == "theta":
        theta, radius, phi = np.meshgrid(
            theta_values[:5], r_values, np.linspace(2.0, 2.2, 9), indexing="ij"
        )
        u, v, w = np.linspace(2.0, 2.2, 9), r_values, theta_values[:5]
        axis_map = {"theta": 0, "r": 1, "phi": 2}
    else:
        radius, theta, phi = np.meshgrid(
            np.linspace(1.02, 1.20, 5),
            np.linspace(0.8, 1.2, 7),
            np.linspace(2.0, 2.2, 9),
            indexing="ij",
        )
        u, v, w = np.linspace(2.0, 2.2, 9), np.linspace(0.8, 1.2, 7), np.linspace(1.02, 1.20, 5)
        axis_map = {"r": 0, "theta": 1, "phi": 2}
    return TOOL.SampleGrid(
        mode=mode,
        points_spherical=np.stack((radius, theta, phi), axis=-1),
        u=u,
        v=v,
        w=w,
        r=radius,
        theta=theta,
        phi=phi,
        axis_map=axis_map,
        basis_u=None,
        basis_v=None,
        basis_w=None,
        metadata={},
    )


def test_spherical_derivatives_for_all_primary_modes() -> None:
    for mode in ("phi", "theta", "r"):
        grid = _spherical_grid(mode)
        # A radial B=C/r^2 has zero divergence and zero curl.
        radial = 1.0 / grid.r**2
        vector = np.stack(
            (radial, np.zeros_like(radial), np.zeros_like(radial)), axis=-1
        )
        curl, divergence = TOOL.spherical_vector_derivatives(vector, grid)
        assert np.allclose(curl, 0.0, atol=2.0e-12)
        assert np.allclose(divergence, 0.0, atol=2.0e-12)


def test_local_fit_recovers_linear_spherical_component_gradients(
    tmp_path: Path,
) -> None:
    import pytest
    simesh_amrvac = pytest.importorskip('simesh.amrvac', reason='optional simesh analysis extension is not installed')
    open_dataset = simesh_amrvac.open_dataset
    write_datfile_from_uniform = simesh_amrvac.write_datfile_from_uniform
    shape = (12, 12, 12)
    xmin = np.array([1.0, 0.6, 2.0])
    xmax = np.array([2.0, 1.2, 2.6])
    axes = [
        xmin[dim]
        + (np.arange(shape[dim]) + 0.5)
        * (xmax[dim] - xmin[dim])
        / shape[dim]
        for dim in range(3)
    ]
    radius, theta, phi = np.meshgrid(*axes, indexing="ij")
    br = 2.0 * radius - theta + 0.5 * phi
    btheta = -0.25 * radius + 1.5 * theta + 0.3 * phi
    bphi = 0.75 * radius + 0.2 * theta - 1.2 * phi
    one = np.ones(shape)
    zero = np.zeros(shape)
    data = np.stack(
        (one, one, zero, zero, zero, br, btheta, bphi), axis=-1
    )
    dat = tmp_path / "linear_spherical.dat"
    par = tmp_path / "amrvac.par"
    write_datfile_from_uniform(
        dat,
        data,
        ["rho", "e", "m1", "m2", "m3", "b1", "b2", "b3"],
        xmin=xmin,
        xmax=xmax,
        block_nx=np.array([4, 4, 4]),
        overwrite=True,
        geometry="spherical",
    )
    par.write_text(
        """
&meshlist
 stretch_dim(1)='none'
/
&mhd_list
 mhd_gamma=1.6666666666666667d0
 mhd_semirelativistic=.false.
/
""",
        encoding="utf-8",
    )
    dataset = open_dataset(str(dat), ghost_width=0, grid_config=str(par))
    r_values = np.array([1.42, 1.50, 1.58])
    theta_values = np.array([0.82, 0.90, 0.98])
    rr, tt = np.meshgrid(r_values, theta_values, indexing="ij")
    pp = np.full_like(rr, 2.30)
    points = np.stack((rr, tt, pp), axis=-1)[None, ...]
    grid = TOOL.SampleGrid(
        mode="phi",
        points_spherical=points,
        u=theta_values,
        v=r_values,
        w=np.array([2.30]),
        r=points[..., 0],
        theta=points[..., 1],
        phi=points[..., 2],
        axis_map={"phi": 0, "r": 1, "theta": 2},
        basis_u=None,
        basis_v=None,
        basis_w=None,
        metadata={},
    )
    raw, _ = TOOL.sample_raw_fields(dataset, grid, "linear", 1000)
    parameters = TOOL.parse_amrvac_par(par)
    arrays, details = TOOL.local_fit_derivatives(
        dataset,
        grid,
        raw,
        parameters,
        dict(TOOL.DEFAULT_UNITS),
        np.array([0.5, 0.75, 1.0]),
        residual_max=0.05,
        sensitivity_max=0.05,
    )
    center_b = np.stack(
        (
            2.0 * rr - tt + 0.5 * pp,
            -0.25 * rr + 1.5 * tt + 0.3 * pp,
            0.75 * rr + 0.2 * tt - 1.2 * pp,
        ),
        axis=-1,
    )
    exact_gradient = np.empty(rr.shape + (3, 3))
    exact_gradient[..., 0, :] = [2.0, -0.25, 0.75]
    exact_gradient[..., 1, :] = [-1.0, 1.5, 0.2]
    exact_gradient[..., 2, :] = [0.5, 0.3, -1.2]
    exact_current, exact_divergence = TOOL.spherical_curl_div_from_partials(
        center_b, exact_gradient, rr, tt
    )
    actual_current = np.stack(
        (arrays["j_r_code"], arrays["j_theta_code"], arrays["j_phi_code"]),
        axis=-1,
    )
    assert np.allclose(actual_current, exact_current, atol=2.0e-10)
    assert np.allclose(arrays["divb_code"], exact_divergence, atol=2.0e-10)
    assert details["cross_level_fraction"] == 0.0
    assert np.array_equal(
        arrays["derivative_reliable_mask"],
        arrays["derivative_strict_reliable_mask"],
    )


def test_cache_round_trip(tmp_path: Path) -> None:
    output = tmp_path / "slice.npz"
    arrays = {"u": np.arange(3.0), "temperature_K": np.ones((2, 3))}
    metadata = {"surface": {"mode": "phi"}, "value": "中文"}
    TOOL.write_cache(output, arrays, metadata)
    loaded, loaded_metadata = TOOL.load_cache(output)
    assert np.array_equal(loaded["u"], arrays["u"])
    assert np.array_equal(loaded["temperature_K"], arrays["temperature_K"])
    assert loaded_metadata == metadata


def test_native_amr_grid_plot_supports_phi_and_theta(tmp_path: Path) -> None:
    dat = tmp_path / "sample_0004.dat"
    par = tmp_path / "amrvac.par"
    dat.write_bytes(b"synthetic DAT identity")
    par.write_text("&meshlist\n/\n", encoding="utf-8")
    dataset = SimpleNamespace(
        ndim=3,
        tree_info=(np.array([1, 2], dtype=int),),
        block_node_coordinates=(
            np.array([[1.00, 1.01, 1.02], [1.00, 1.01, 1.02]]),
            np.array([[1.00, 1.10, 1.20], [1.20, 1.30, 1.40]]),
            np.array([[2.00, 2.10, 2.20], [2.00, 2.10, 2.20]]),
        ),
        block_nx=np.array([2, 2, 2], dtype=int),
        metadata={
            "xmin": np.array([1.00, 1.00, 2.00]),
            "xmax": np.array([1.02, 1.40, 2.20]),
            "time": 0.02,
            "it": 17,
            "geometry": "spherical",
            "coordinate_source": "synthetic true nodes",
        },
    )

    phi_output = tmp_path / "phi_grid.png"
    phi_report = TOOL.plot_native_amr_grid(
        dataset,
        dat,
        par,
        phi_output,
        mode="phi",
        value_deg=np.degrees(2.15),
        min_height_mm=0.0,
        max_height_mm=20.0,
        u_min=None,
        u_max=None,
        dpi=80,
        unit_length_cm=TOOL.DEFAULT_UNITS["unit_length_cm"],
        unit_time_s=TOOL.DEFAULT_UNITS["unit_time_s"],
        marker_u=20.0,
        marker_v=5.0,
        marker_label="test point",
    )
    assert phi_output.is_file()
    assert phi_output.with_suffix(".json").is_file()
    assert phi_report["slice_leaf_blocks"] == 2
    assert phi_report["slice_leaf_counts_by_level"] == {"1": 1, "2": 1}
    assert phi_report["slice"]["marker"]["label"] == "test point"
    assert json.loads(phi_output.with_suffix(".json").read_text())["tool"][
        "version"
    ] == "1.2.0"

    theta_output = tmp_path / "theta_grid.png"
    theta_report = TOOL.plot_native_amr_grid(
        dataset,
        dat,
        par,
        theta_output,
        mode="theta",
        value_deg=np.degrees(1.15),
        min_height_mm=0.0,
        max_height_mm=20.0,
        u_min=None,
        u_max=None,
        dpi=80,
        unit_length_cm=TOOL.DEFAULT_UNITS["unit_length_cm"],
        unit_time_s=TOOL.DEFAULT_UNITS["unit_time_s"],
    )
    assert theta_output.is_file()
    assert theta_report["slice_leaf_blocks"] == 1
    assert theta_report["slice_leaf_counts_by_level"] == {"1": 1, "2": 0}

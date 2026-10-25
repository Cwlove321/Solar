from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from simulation_analysis.config import DEFAULT_CONFIG, load_config
from simulation_analysis.geometry import axis_edges
from simulation_analysis.logs import fortran_float, read_rows
from simulation_analysis.physics import derive_state
from simulation_analysis.rbsl_balance import (
    flux_from_rbsl_current,
    rbsl_current_from_flux,
    thin_torus_balance_current,
)


class FakeDataset:
    physical_domain = np.array([[1.0, 1.0, 2.0], [1.4, 1.4, 2.8]])
    domain_nx = np.array([4, 4, 4])


def test_geometric_edges_close_domain():
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    edges = axis_edges(FakeDataset(), 0, 3, config)
    assert edges.size == 17
    assert np.isclose(edges[0], 1.0)
    assert np.isclose(edges[-1], 1.4)
    assert np.all(np.diff(edges) > 0.0)


def test_internal_energy_classical_state():
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    config["physics"]["momentum_model"] = "classical"
    shape = (2, 2)
    raw = {
        "rho": np.ones(shape),
        "m1": np.full(shape, 2.0),
        "m2": np.zeros(shape),
        "m3": np.zeros(shape),
        "energy": np.full(shape, 3.0),
        "b1": np.ones(shape),
        "b2": np.zeros(shape),
        "b3": np.zeros(shape),
    }
    state = derive_state(raw, 5.0 / 3.0, config)
    assert np.allclose(state["v1_code"], 2.0)
    assert np.allclose(state["pressure_code"], 2.0)
    assert np.allclose(state["temperature_K"], 2.0e6)
    assert np.allclose(state["beta"], 4.0)


def test_semirelativistic_parallel_velocity_is_unchanged():
    config = json.loads(json.dumps(DEFAULT_CONFIG))
    shape = (1,)
    raw = {
        "rho": np.ones(shape),
        "m1": np.full(shape, 2.0),
        "m2": np.zeros(shape),
        "m3": np.zeros(shape),
        "energy": np.ones(shape),
        "b1": np.ones(shape),
        "b2": np.zeros(shape),
        "b3": np.zeros(shape),
    }
    state = derive_state(raw, 5.0 / 3.0, config)
    assert np.allclose(state["v1_code"], 2.0)
    assert np.all(state["alfven_reduced_km_s"] <= config["physics"]["reduced_c_km_s"])


def test_config_and_fortran_log(tmp_path: Path):
    run = tmp_path / "run"
    (run / "data").mkdir(parents=True)
    config_path = tmp_path / "analysis.json"
    config_path.write_text(json.dumps({"run_dir": "run"}))
    config = load_config(config_path)
    assert Path(config["run_dir"]) == run
    assert fortran_float("1.23D-04") == 1.23e-4
    values = [
        "10", "1.0E-03", "2.0E-06", "3", "4", "5", "6", "7", "8", "9", "10",
        "11", "12", "13", "14", "15", "16", "17", "18"
    ]
    (run / "data" / "run_.log").write_text(" ".join(values) + " | ignored\n")
    rows = read_rows(config)
    assert len(rows) == 1
    assert rows[0]["it"] == 10
    assert rows[0]["dt"] == 2.0e-6
    assert rows[0]["momentum_extreme"] == 6.0


def test_rbsl_flux_current_round_trip():
    flux = 2.1e20
    current = rbsl_current_from_flux(flux, 15.0, 695.5, 1.9977023903979652, -1)
    recovered = flux_from_rbsl_current(current, 15.0, 695.5, 1.9977023903979652)
    assert current < 0.0
    assert np.isclose(recovered, flux)


def test_thin_torus_balance_flux_matches_rbsl14_scale():
    unit_b = 1.9977023903979652
    current, denominator = thin_torus_balance_current(
        np.array([-5.782502902955166 / unit_b]), 152.85967269004476, 15.0, 695.5
    )
    flux = flux_from_rbsl_current(current, 15.0, 695.5, unit_b)
    assert denominator > 0.0
    assert np.isclose(flux[0], 2.0785037296510874e20)

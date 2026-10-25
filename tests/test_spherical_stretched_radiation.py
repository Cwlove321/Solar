import importlib.util
from pathlib import Path

import numpy as np
import pytest


SCRIPT = (
    Path(__file__).resolve().parents[1]
    / "workflows"
    / "Synthesize_Radiation"
    / "synthesize_aia_spherical_stretched.py"
)
SPEC = importlib.util.spec_from_file_location("spherical_stretched_radiation", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class UnitVolume:
    left_edge = np.array([0.0, 0.0, 0.0])
    right_edge = np.array([1.0, 1.0, 1.0])
    dims = np.array([2, 2, 2])
    dds = np.array([0.5, 0.5, 0.5])
    idds = 1.0 / dds
    kappa_view = np.zeros((2, 2, 2))
    j_view = np.ones((2, 2, 2))


class Ray:
    def __init__(self):
        self.sum_val = 0.0
        self.tau = 0.0


@pytest.mark.parametrize(
    ("start", "direction", "expected_length"),
    [
        ([2.0, 0.5, 0.5], [-1.0, 0.0, 0.0], 1.0),
        ([2.0, 2.0, 2.0], -np.ones(3) / np.sqrt(3.0), np.sqrt(3.0)),
    ],
)
def test_centimetre_dda_integrates_full_box_length(start, direction, expected_length):
    ray = Ray()

    hit = MODULE.walk_volume_cm(
        UnitVolume(), np.asarray(start), np.asarray(direction), ray
    )

    assert hit
    assert ray.sum_val == pytest.approx(expected_length)
    assert ray.tau == 0.0


def test_centimetre_dda_rejects_ray_that_misses_box():
    ray = Ray()

    hit = MODULE.walk_volume_cm(
        UnitVolume(), np.array([2.0, 2.0, 0.5]), np.array([-1.0, 0.0, 0.0]), ray
    )

    assert not hit
    assert ray.sum_val == 0.0

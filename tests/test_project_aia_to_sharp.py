from contextlib import contextmanager
import importlib.util
from pathlib import Path
from types import SimpleNamespace

import numpy as np


SCRIPT = Path(__file__).parents[1] / "workflows" / "observations" / "project_aia_to_sharp.py"
SPEC = importlib.util.spec_from_file_location("project_aia_to_sharp", SCRIPT)
MODULE = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(MODULE)


class FakeAiaMap:
    def __init__(self):
        self.calls = []

    def reproject_to(self, wcs, **kwargs):
        self.calls.append((wcs, kwargs))
        return "projected"


def test_reproject_to_sharp_applies_differential_rotation(monkeypatch):
    entered_models = []

    @contextmanager
    def fake_surface_propagation(*, rotation_model):
        entered_models.append(rotation_model)
        yield

    monkeypatch.setattr(MODULE, "propagate_with_solar_surface", fake_surface_propagation)
    aia = FakeAiaMap()
    target = SimpleNamespace(wcs="sharp-cea-wcs", data=np.empty((12, 34)))

    result = MODULE.reproject_to_sharp(aia, target, rotation_model="snodgrass")

    assert result == "projected"
    assert entered_models == ["snodgrass"]
    assert aia.calls == [
        ("sharp-cea-wcs", {"shape_out": (12, 34), "algorithm": "interpolation"})
    ]


def test_reproject_to_sharp_can_skip_differential_rotation(monkeypatch):
    def unexpected_surface_propagation(*args, **kwargs):
        raise AssertionError("surface propagation should be disabled")

    monkeypatch.setattr(MODULE, "propagate_with_solar_surface", unexpected_surface_propagation)
    aia = FakeAiaMap()
    target = SimpleNamespace(wcs="sharp-cea-wcs", data=np.empty((5, 7)))

    result = MODULE.reproject_to_sharp(aia, target, differential_rotation=False)

    assert result == "projected"
    assert aia.calls == [
        ("sharp-cea-wcs", {"shape_out": (5, 7), "algorithm": "interpolation"})
    ]

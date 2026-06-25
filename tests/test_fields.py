import numpy as np

from solarutils.fields.flux import magnetic_flux_from_data
from solarutils.fields.helicity import volume_integral
from solarutils.fields.vector import curl


def test_curl_of_simple_cartesian_field_is_constant():
    z, y, x = np.mgrid[0:4, 0:5, 0:6]
    b1 = -0.5 * y
    b2 = 0.5 * x
    b3 = np.zeros_like(x, dtype=float)

    cx, cy, cz = curl(b1, b2, b3, axis_order="zyx")

    assert np.allclose(cx, 0)
    assert np.allclose(cy, 0)
    assert np.allclose(cz, 1)


def test_volume_integral_of_constant_field_matches_box_volume():
    array = np.ones((4, 5, 6))

    assert volume_integral(array, dx=2.0) == 3 * 4 * 5 * 8


def test_magnetic_flux_from_data_supports_common_modes():
    data = np.array([[-2.0, -1.0, np.nan], [0.0, 3.0, 4.0]])

    assert magnetic_flux_from_data(data, 10.0, mode="unsigned") == 100.0
    assert magnetic_flux_from_data(data, 10.0, mode="positive") == 70.0
    assert magnetic_flux_from_data(data, 10.0, mode="negative") == -30.0


def test_magnetic_flux_from_data_rejects_unknown_mode():
    try:
        magnetic_flux_from_data(np.ones((2, 2)), 1.0, mode="net")
    except ValueError as exc:
        assert "mode must be" in str(exc)
    else:
        raise AssertionError("unknown magnetic flux mode should fail")

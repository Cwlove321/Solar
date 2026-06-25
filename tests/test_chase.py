import numpy as np

import solarutils.obs.chase as chase


def test_fit_absorption_gaussian_center_recovers_line_center():
    wavelength = np.linspace(6561.5, 6564.0, 80)
    true_center = 6562.83
    profile = 1000.0 - 300.0 * np.exp(-0.5 * ((wavelength - true_center) / 0.28) ** 2)

    center = chase.fit_absorption_gaussian_center(wavelength, profile)

    assert np.isclose(center, true_center, atol=1e-3)


def test_centroid_line_center_uses_absorption_weights():
    wavelength = np.array([1.0, 2.0, 3.0])
    profile = np.array([10.0, 5.0, 10.0])

    center = chase.centroid_line_center(wavelength, profile)

    assert center == 2.0


def test_centroid_line_center_map_is_vectorized():
    wavelength = np.array([1.0, 2.0, 3.0])
    data = np.array(
        [
            [[10.0, 10.0], [10.0, 8.0]],
            [[5.0, 8.0], [10.0, 4.0]],
            [[10.0, 6.0], [8.0, 8.0]],
        ]
    )

    centers = chase.centroid_line_center_map(wavelength, data)

    assert centers.shape == (2, 2)
    assert np.isclose(centers[0, 0], 2.0)
    assert np.isclose(centers[0, 1], 8.0 / 3.0)
    assert np.isclose(centers[1, 1], 2.0)


def test_doppler_velocity_selects_nearest_wavelength_index():
    line_center = 6562.8
    wavelength = line_center + np.arange(-3, 4) * 0.05
    velocity = (wavelength[5] / line_center - 1.0) * chase.SPEED_OF_LIGHT_KMS

    index = chase.chase_channel_index_for_velocity(wavelength, line_center, velocity)

    assert index == 5


def test_chase_wavelength_axis_uses_fits_linear_axis():
    header = {"NAXIS3": 4, "CRVAL3": 10.0, "CDELT3": 0.25}

    wavelength = chase.chase_wavelength_axis(header)

    assert np.allclose(wavelength, [10.0, 10.25, 10.5, 10.75])


def test_differential_rotation_rate_matches_equatorial_coefficient():
    omega = chase.differential_rotation_rate(np.array([0.0]))
    expected = chase.SNODGRASS_ROTATION_COEFFS_DEG_PER_DAY[0] * np.pi / 180.0 / 86400.0

    assert np.isclose(omega[0], expected)

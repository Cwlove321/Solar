from pathlib import Path

import numpy as np
import pytest

from solarutils.sim.radiation import (
    K_B_CGS,
    MASS_H_CGS,
    RadiationPlane,
    clean_cube_for_dda,
    density_pressure_temperature_from_amrvac,
    generate_detector_points,
    log_intensity_image,
    orient_image,
    rectangular_flux_series,
    time_distance_from_npy,
)


def test_radiation_plane_notebook_defaults():
    plane = RadiationPlane.from_notebook_defaults(normal=(1, 0, 0), resolution=128)

    assert plane.resolution == 128
    assert np.allclose(plane.normal, [1, 0, 0])
    assert np.allclose(plane.left_edge, np.array([-3, -3, 0.1]) * 1e9)
    assert plane.width == pytest.approx(12e9)


def test_rectangular_detector_points_use_requested_shape():
    plane = RadiationPlane(
        left_edge=np.array([-1, -2, 0]),
        right_edge=np.array([1, 2, 3]),
        center=np.array([0, 0, 1.5]),
        normal=np.array([1, 0, 0]),
        width=(4.0, 3.0),
        resolution=(5, 4),
    )

    points, shape = generate_detector_points(plane)

    assert points.shape == (20, 3)
    assert shape == (4, 5)
    assert np.allclose(points[:, 0], 0.0)


def test_density_pressure_temperature_conversion_matches_notebook():
    rho = np.array([[[1.4 * MASS_H_CGS]]])
    e = np.array([[[1.5]]])

    number_density, pressure, temperature = density_pressure_temperature_from_amrvac(rho, e)

    assert number_density[0, 0, 0] == pytest.approx(2.341e-15)
    assert pressure[0, 0, 0] == pytest.approx(0.317538)
    assert temperature[0, 0, 0] == pytest.approx(0.317538 / (2.341e-15 * 2.3 * K_B_CGS))


def test_clean_cube_transposes_and_removes_invalid_values():
    cube = np.array([[[1.0, np.nan]], [[np.inf, -np.inf]]])

    cleaned = clean_cube_for_dda(cube)

    assert cleaned.shape == (2, 1, 2)
    assert np.all(np.isfinite(cleaned))
    assert cleaned[0, 0, 0] == 1.0


def test_log_intensity_image_only_logs_positive_values():
    image = np.array([[1.0, 10.0], [0.0, -1.0]])

    logged = log_intensity_image(image, fill_value=-99.0)

    assert np.allclose(logged, [[0.0, 1.0], [-99.0, -99.0]])


def test_orient_image_flip_and_crop():
    image = np.arange(12).reshape(3, 4)

    oriented = orient_image(image, crop_pixels=1)

    assert np.array_equal(oriented, np.array([[9, 10], [5, 6], [1, 2]]))


def test_time_distance_and_flux_from_npy(tmp_path: Path):
    first = np.arange(16, dtype=float).reshape(4, 4)
    second = first + 10
    paths = []
    for index, image in enumerate([first, second]):
        path = tmp_path / f"image_{index:04d}.npy"
        np.save(path, image)
        paths.append(path)

    td = time_distance_from_npy(paths, (0, 0), (3, 0), crop_pixels=0)
    flux = rectangular_flux_series(paths, y_slice=(0, 2), x_slice=(0, 2))

    assert td.shape == (2, 3)
    assert np.allclose(td[0], [12.0, 13.5, 15.0])
    assert np.allclose(flux, [42.0, 82.0])

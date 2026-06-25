import numpy as np

from solarutils.obs.time_distance import line_distance_km, sample_line_data, time_distance_series


def test_sample_line_data_includes_both_endpoints():
    data = np.arange(25, dtype=float).reshape(5, 5)

    sampled = sample_line_data(data, (0, 0), (4, 4), 5, order=0)

    assert np.array_equal(sampled, np.array([0, 6, 12, 18, 24], dtype=float))


def test_sample_line_data_rejects_too_few_samples():
    try:
        sample_line_data(np.ones((2, 2)), (0, 0), (1, 1), 1)
    except ValueError as exc:
        assert "at least 2" in str(exc)
    else:
        raise AssertionError("n_samples < 2 should fail")


def test_time_distance_series_rejects_empty_paths():
    try:
        time_distance_series([], (0, 0), (1, 1), 2)
    except ValueError as exc:
        assert "at least one" in str(exc)
    else:
        raise AssertionError("empty paths should fail")


def test_line_distance_km_uses_observer_distance():
    smap = type("MapLike", (), {"dsun": 1.0})()

    class Distance:
        def to_value(self, unit):
            assert unit == "km"
            return 1.0

    smap.dsun = Distance()

    assert np.isclose(line_distance_km(smap, (0, 0), (206264.80624709636, 0)), 1.0)

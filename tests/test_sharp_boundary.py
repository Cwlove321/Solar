from datetime import datetime

import pytest

from solarutils.obs.sharp_boundary import _observation_time


@pytest.mark.parametrize(
    ("metadata", "expected"),
    [
        ({"t_obs": "2024.06.01_18:22:45_TAI"}, datetime(2024, 6, 1, 18, 22, 45)),
        ({"t_rec": "2024.06.01_18:24:00_TAI"}, datetime(2024, 6, 1, 18, 24, 0)),
        ({"date-obs": "2024-06-01T18:22:45.000"}, datetime(2024, 6, 1, 18, 22, 45)),
    ],
)
def test_observation_time_fallbacks(metadata, expected):
    assert _observation_time(metadata) == expected


def test_observation_time_requires_supported_keyword():
    with pytest.raises(KeyError, match="T_OBS"):
        _observation_time({})

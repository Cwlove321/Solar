from pathlib import Path
from types import SimpleNamespace

import astropy.units as u
from astropy.time import Time
import pytest
from sunpy.coordinates import get_body_heliographic_stonyhurst

from workflows.Cartesian_NLFFF.Preprocess.make_cea_patch import infer_p_angle_degree
from workflows.Cartesian_NLFFF.Preprocess.prepare_hmi_b720_timeseries import (
    discover_records,
    reference_coordinate,
    rotate_center,
    select_reference_timestamp,
)


SEGMENTS = ("field", "inclination", "azimuth", "disambig")


def touch_record(root: Path, timestamp: str, *, series_letter: str = "B") -> None:
    for segment in SEGMENTS:
        (root / f"hmi.{series_letter}_720s.{timestamp}.{segment}.fits").touch()


def test_discover_records_accepts_series_letter_case_and_orders_records(tmp_path: Path) -> None:
    later = "20110907_001200_TAI"
    earlier = "20110906_221200_TAI"
    touch_record(tmp_path, later, series_letter="b")
    touch_record(tmp_path, earlier, series_letter="B")
    (tmp_path / "._hmi.B_720s.20110906_221200_TAI.field.fits").touch()

    records = discover_records(tmp_path)

    assert list(records) == [earlier, later]
    assert set(records[earlier]) == set(SEGMENTS)
    assert select_reference_timestamp(records, None) == earlier
    assert select_reference_timestamp(records, later.lower()) == later


def test_discover_records_rejects_an_incomplete_timestamp(tmp_path: Path) -> None:
    timestamp = "20110906_221200_TAI"
    touch_record(tmp_path, timestamp)
    (tmp_path / f"hmi.B_720s.{timestamp}.disambig.fits").unlink()

    with pytest.raises(ValueError, match="incomplete timestamps"):
        discover_records(tmp_path)


def test_rotate_center_tracks_one_day_of_solar_rotation() -> None:
    start = Time("2011-09-06T22:12:00")
    reference_map = SimpleNamespace(date=start, reference_date=start, rsun_meters=695700 * u.km)
    target_map = SimpleNamespace(
        date=start + 1 * u.day,
        reference_date=start + 1 * u.day,
        observer_coordinate=get_body_heliographic_stonyhurst("earth", start + 1 * u.day),
    )
    center = reference_coordinate(reference_map, center_lon=0.0, center_lat=20.0)

    longitude, latitude = rotate_center(center, target_map, model="howard")

    assert 12.0 < longitude < 14.5
    assert latitude == pytest.approx(20.0, abs=1e-10)


def test_hmi_p_angle_uses_clockwise_equivalent_of_crota2() -> None:
    smap = SimpleNamespace(meta={"crota2": 180.013443})

    assert infer_p_angle_degree(smap) == pytest.approx(-180.013443)

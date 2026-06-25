import astropy.units as u
from astropy.table import Table

from solarutils.obs.download import (
    HMI_SERIES_ALIASES,
    _sample_query_response,
    _sample_sharp_query_result,
    _with_sample,
    carrington_rotation_for_time,
    make_sharp_segment_recordset,
    make_sharp_time_recordset,
    resolve_hmi_series,
)


def test_resolve_hmi_series_aliases():
    assert resolve_hmi_series("hmi_m_45s") == "hmi.M_45s"
    assert resolve_hmi_series("hmi.45") == "hmi.M_45s"
    assert resolve_hmi_series("hmi_m_720s") == "hmi.M_720s"
    assert resolve_hmi_series("sharp") == "hmi.sharp_cea_720s"
    assert resolve_hmi_series("hmi_vector") == "hmi.B_720s"
    assert resolve_hmi_series("hmi_synoptic_polfil") == "hmi.mrdailysynframe_polfil_720s"
    assert resolve_hmi_series("hmi_synoptic_cr_polfil") == "hmi.synoptic_mr_polfil_720s"


def test_unknown_hmi_series_passes_through():
    assert resolve_hmi_series("custom.series") == "custom.series"


def test_hmi_alias_table_keeps_legacy_hmi75_spelling():
    assert HMI_SERIES_ALIASES["hmi.75"] == "hmi.M_720s"


def test_hmi_alias_table_includes_vector_and_synoptic_products():
    assert HMI_SERIES_ALIASES["hmi_b_720s"] == "hmi.B_720s"
    assert HMI_SERIES_ALIASES["hmi_mr_polfil_720s"] == "hmi.mrdailysynframe_polfil_720s"
    assert HMI_SERIES_ALIASES["hmi_synoptic_mr_polfil_720s"] == "hmi.synoptic_mr_polfil_720s"


def test_with_sample_leaves_query_unchanged_without_sample():
    assert _with_sample("query", object(), u, None) == "query"


def test_sample_query_response_keeps_near_cadence_rows():
    table = Table(
        {
            "Start time": [
                "2024-10-08T23:00:55.275",
                "2024-10-08T23:10:55.276",
                "2024-10-08T23:20:55.277",
                "2024-10-08T23:30:55.293",
                "2024-10-08T23:40:55.279",
                "2024-10-08T23:50:55.281",
            ],
        }
    )

    sampled = _sample_query_response([table], 10)

    assert len(sampled) == 6


def test_carrington_rotation_for_time_returns_integer_rotation():
    assert carrington_rotation_for_time("2024-10-08T12:00:00") == 2289


def test_make_sharp_time_recordset_converts_iso_range_to_jsoc_tai():
    recordset = make_sharp_time_recordset(
        "hmi.sharp_cea_720s",
        "2015-08-27T05:23:00",
        "2015-08-27T05:25:00",
        harpnum=5885,
    )

    assert recordset == "hmi.sharp_cea_720s[5885][2015.08.27_05:23:00_TAI-2015.08.27_05:25:00_TAI]"


def test_make_sharp_time_recordset_all_harps_uses_empty_primekey():
    recordset = make_sharp_time_recordset(
        "hmi.sharp_cea_720s",
        "2015.08.27_05:24:00_TAI",
        "2015.08.27_05:36:00_TAI",
    )

    assert recordset == "hmi.sharp_cea_720s[][2015.08.27_05:24:00_TAI-2015.08.27_05:36:00_TAI]"


def test_make_sharp_segment_recordset_selects_segments():
    recordset = make_sharp_segment_recordset(
        "hmi.sharp_cea_720s",
        5885,
        "2015.08.27_05:24:00_TAI",
        segments=("Br", "Bt", "Bp"),
    )

    assert recordset == "hmi.sharp_cea_720s[5885][2015.08.27_05:24:00_TAI]{Br,Bt,Bp}"


def test_sample_sharp_query_result_keeps_requested_cadence():
    table = Table(
        {
            "T_REC": [
                "2015.08.27_05:24:00_TAI",
                "2015.08.27_05:36:00_TAI",
                "2015.08.27_05:48:00_TAI",
            ],
            "HARPNUM": [5885, 5885, 5885],
        }
    ).to_pandas()

    sampled = _sample_sharp_query_result(table, 20)

    assert sampled["T_REC"].tolist() == [
        "2015.08.27_05:24:00_TAI",
        "2015.08.27_05:48:00_TAI",
    ]

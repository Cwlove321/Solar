"""Search and download helpers for common solar-observation data sources."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
from typing import Iterable, Sequence


HMI_SERIES_ALIASES = {
    "hmi_45s": "hmi.M_45s",
    "hmi_m_45s": "hmi.M_45s",
    "hmi.45": "hmi.M_45s",
    "hmi_720s": "hmi.M_720s",
    "hmi_m_720s": "hmi.M_720s",
    "hmi.720": "hmi.M_720s",
    "hmi.75": "hmi.M_720s",
    "hmi_b_720s": "hmi.B_720s",
    "hmi.b720": "hmi.B_720s",
    "hmi_b": "hmi.B_720s",
    "hmi_vector_720s": "hmi.B_720s",
    "hmi_vector": "hmi.B_720s",
    "hmi_mrdailysynframe_polfil_720s": "hmi.mrdailysynframe_polfil_720s",
    "hmi_mr_polfil_720s": "hmi.mrdailysynframe_polfil_720s",
    "hmi_synoptic_polfil": "hmi.mrdailysynframe_polfil_720s",
    "hmi_synoptic_mr_polfil_720s": "hmi.synoptic_mr_polfil_720s",
    "hmi_synoptic_cr_polfil": "hmi.synoptic_mr_polfil_720s",
    "sharp": "hmi.sharp_cea_720s",
    "sharp_cea": "hmi.sharp_cea_720s",
}

DEFAULT_AIA_SEGMENTS = ("image",)
DEFAULT_HMI_SEGMENTS = ("magnetogram",)
DEFAULT_HMI_VECTOR_SEGMENTS = ("field", "inclination", "azimuth", "disambig")
DEFAULT_HMI_MRDAILYSYNFRAME_SEGMENTS = ("Mr_polfil",)
DEFAULT_SHARP_SEGMENTS = ("Br", "Bp", "Bt", "bitmap", "conf_disambig")


@dataclass(frozen=True)
class SearchPlan:
    """A small record describing an observation-data search."""

    provider: str
    start: str
    end: str
    description: str
    output_dir: Path
    query_result: object
    fetcher: str = "fido"
    email: str | None = None
    export_recordsets: tuple[str, ...] = ()


def search_sdo_aia(
    start: str,
    end: str,
    email: str,
    wavelengths: Sequence[int] = (171, 193, 211, 304),
    segments: Sequence[str] = DEFAULT_AIA_SEGMENTS,
    sample_minutes: float | None = None,
    output_dir: str | Path = "data/raw/sdo/aia",
) -> SearchPlan:
    """Search JSOC for SDO/AIA level-1 EUV images."""

    Fido, a, u = _sunpy_imports()
    query = a.Time(start, end) & a.jsoc.Series("aia.lev1_euv_12s") & a.jsoc.Notify(email)
    query &= _or_attrs(a.jsoc.Segment, segments)
    query &= _or_attrs(lambda value: a.jsoc.PrimeKey("WAVELNTH", int(value)), wavelengths)
    query = _with_sample(query, a, u, sample_minutes)
    result = Fido.search(query)
    return SearchPlan("JSOC", start, end, f"SDO/AIA wavelengths={tuple(wavelengths)} Angstrom", Path(output_dir), result)


def search_hmi(
    start: str,
    end: str,
    email: str,
    series: str = "hmi_m_45s",
    segments: Sequence[str] = DEFAULT_HMI_SEGMENTS,
    sample_minutes: float | None = None,
    output_dir: str | Path = "data/raw/sdo/hmi",
) -> SearchPlan:
    """Search JSOC for HMI data.

    Common aliases include `hmi_m_45s`, `hmi_m_720s`, and `hmi_b_720s`.
    """

    Fido, a, _ = _sunpy_imports()
    resolved_series = resolve_hmi_series(series)
    query = a.Time(start, end) & a.jsoc.Series(resolved_series) & a.jsoc.Notify(email)
    if segments:
        query &= _or_attrs(a.jsoc.Segment, segments)
    query = _with_sample(query, a, _, sample_minutes)
    result = Fido.search(query)
    return SearchPlan("JSOC", start, end, f"HMI series={resolved_series} segments={tuple(segments)}", Path(output_dir), result)


def search_hmi_vector(
    start: str,
    end: str,
    email: str,
    segments: Sequence[str] = DEFAULT_HMI_VECTOR_SEGMENTS,
    sample_minutes: float | None = None,
    output_dir: str | Path = "data/raw/sdo/hmi_b_720s",
) -> SearchPlan:
    """Search JSOC for full-disk HMI vector magnetic-field data."""

    return search_hmi(
        start,
        end,
        email,
        series="hmi_b_720s",
        segments=segments,
        sample_minutes=sample_minutes,
        output_dir=output_dir,
    )


def search_hmi_mrdailysynframe_polfil(
    start: str,
    end: str,
    email: str,
    segments: Sequence[str] = DEFAULT_HMI_MRDAILYSYNFRAME_SEGMENTS,
    sample_minutes: float | None = None,
    output_dir: str | Path = "data/raw/sdo/hmi_mrdailysynframe_polfil_720s",
) -> SearchPlan:
    """Search JSOC for HMI daily synoptic-frame radial-field maps."""

    return search_hmi(
        start,
        end,
        email,
        series="hmi_mrdailysynframe_polfil_720s",
        segments=segments,
        sample_minutes=sample_minutes,
        output_dir=output_dir,
    )


def search_hmi_synoptic_mr_polfil(
    car_rotation: int,
    email: str,
    segments: Sequence[str] = DEFAULT_HMI_MRDAILYSYNFRAME_SEGMENTS,
    output_dir: str | Path = "data/raw/sdo/hmi_synoptic_mr_polfil_720s",
) -> SearchPlan:
    """Search JSOC for Carrington-rotation HMI synoptic radial-field maps."""

    Fido, a, _ = _sunpy_imports()
    resolved_series = resolve_hmi_series("hmi_synoptic_mr_polfil_720s")
    query = a.jsoc.Series(resolved_series) & a.jsoc.Notify(email)
    query &= a.jsoc.PrimeKey("CAR_ROT", int(car_rotation))
    if segments:
        query &= _or_attrs(a.jsoc.Segment, segments)
    result = Fido.search(query)
    return SearchPlan(
        "JSOC",
        str(car_rotation),
        str(car_rotation),
        f"HMI series={resolved_series} CAR_ROT={int(car_rotation)} segments={tuple(segments)}",
        Path(output_dir),
        result,
    )


def carrington_rotation_for_time(time: str) -> int:
    """Return the integer Carrington rotation containing `time`."""

    from sunpy.coordinates.sun import carrington_rotation_number

    return int(math.floor(float(carrington_rotation_number(time))))


def search_sharp(
    start: str,
    end: str,
    email: str,
    harpnum: int | None = None,
    series: str = "hmi.sharp_cea_720s",
    segments: Sequence[str] = DEFAULT_SHARP_SEGMENTS,
    sample_minutes: float | None = None,
    output_dir: str | Path = "data/raw/sdo/sharp",
) -> SearchPlan:
    """Search JSOC for SHARP vector-magnetic-field data through DRMS."""

    import drms

    resolved_series = resolve_hmi_series(series)
    recordset = make_sharp_time_recordset(resolved_series, start, end, harpnum=harpnum)
    client = drms.Client(email=email)
    result = client.query(recordset, key="T_REC,HARPNUM,NOAA_AR,CRVAL1,CRVAL2,CRLN_OBS,CRLT_OBS")
    result = _sample_sharp_query_result(result, sample_minutes)
    export_recordsets = tuple(
        make_sharp_segment_recordset(resolved_series, int(row.HARPNUM), str(row.T_REC), segments=segments)
        for row in result.itertuples(index=False)
    )
    harp_text = f" harpnum={harpnum}" if harpnum is not None else ""
    return SearchPlan(
        "JSOC-DRMS",
        start,
        end,
        f"SHARP series={resolved_series}{harp_text} records={len(export_recordsets)}",
        Path(output_dir),
        result,
        fetcher="drms_export",
        email=email,
        export_recordsets=export_recordsets,
    )


def search_stereo_a(
    start: str,
    end: str,
    instrument: str = "SECCHI",
    detector: str = "EUVI",
    wavelength: int | None = None,
    sample_minutes: float | None = None,
    output_dir: str | Path = "data/raw/stereo_a",
) -> SearchPlan:
    """Search VSO through SunPy/Fido for STEREO-A observations."""

    Fido, a, u = _sunpy_imports()
    query = a.Time(start, end) & a.Source("STEREO_A") & a.Instrument(instrument)
    if detector:
        query &= a.Detector(detector)
    if wavelength is not None:
        query &= a.Wavelength(int(wavelength) * u.AA)
    query = _with_sample(query, a, u, sample_minutes)
    result = Fido.search(query)
    return SearchPlan("VSO", start, end, f"STEREO-A {instrument}/{detector} wavelength={wavelength}", Path(output_dir), result)


def search_solar_orbiter(
    start: str,
    end: str,
    instrument: str = "EUI",
    product: str | None = None,
    level: str | None = "L1",
    sample_minutes: float | None = None,
    output_dir: str | Path = "data/raw/solar_orbiter",
) -> SearchPlan:
    """Search the Solar Orbiter Archive through the optional `sunpy-soar` plugin."""

    try:
        import sunpy_soar  # noqa: F401
    except ImportError as exc:
        raise ImportError(
            "Solar Orbiter searches require the optional dependency `sunpy-soar`. "
            "Install with `pip install sunpy-soar` or `pip install .[solar-orbiter]`."
        ) from exc

    Fido, a, _ = _sunpy_imports()
    query = a.Time(start, end) & a.Instrument(instrument)
    if level is not None:
        query &= a.Level(level)
    if product:
        query &= a.soar.Product(product)
    result = Fido.search(query)
    result = _sample_query_response(result, sample_minutes)
    product_text = f" product={product}" if product else ""
    return SearchPlan("SOAR", start, end, f"Solar Orbiter {instrument}{product_text}", Path(output_dir), result)


def fetch_results(plan: SearchPlan, overwrite: bool = False):
    """Download files from a `SearchPlan` into its output directory."""

    if plan.fetcher == "drms_export":
        return fetch_drms_export_results(plan, overwrite=overwrite)

    from sunpy.net import Fido

    plan.output_dir.mkdir(parents=True, exist_ok=True)
    path_template = str(plan.output_dir / "{file}")
    return Fido.fetch(plan.query_result, path=path_template, overwrite=overwrite)


def fetch_drms_export_results(
    plan: SearchPlan,
    overwrite: bool = False,
    *,
    method: str = "url",
    protocol: str = "fits",
    timeout: float = 900,
    sleep: float = 10,
):
    """Submit DRMS export requests and download the resulting files."""

    if not plan.email:
        raise ValueError("DRMS export downloads require a JSOC email address")
    if not plan.export_recordsets:
        return []

    import drms
    import pandas as pd

    plan.output_dir.mkdir(parents=True, exist_ok=True)
    client = drms.Client(email=plan.email)
    downloaded = []
    for recordset in plan.export_recordsets:
        request = client.export(recordset, method=method, protocol=protocol, email=plan.email)
        request.wait(timeout=timeout, sleep=sleep)
        if request.has_failed():
            raise RuntimeError(f"DRMS export failed for {recordset}: status={request.status}")
        existing_rows = []
        missing_indexes = []
        for index, row in request.data.iterrows():
            output_path = plan.output_dir / str(row["filename"])
            if output_path.exists() and not overwrite:
                url = request.urls.loc[index, "url"] if index in request.urls.index else ""
                existing_rows.append({"record": row["record"], "url": url, "download": str(output_path)})
                continue
            if output_path.exists() and overwrite:
                output_path.unlink()
            missing_indexes.append(index)

        if missing_indexes:
            files = request.download(str(plan.output_dir), index=missing_indexes, timeout=timeout)
            downloaded.append(files)
        if existing_rows:
            downloaded.append(pd.DataFrame(existing_rows))

    if not downloaded:
        return []
    return pd.concat(downloaded, ignore_index=True)


def resolve_hmi_series(series: str) -> str:
    """Resolve a short HMI alias to a JSOC series name."""

    return HMI_SERIES_ALIASES.get(series.lower(), series)


def make_sharp_time_recordset(series: str, start: str, end: str, *, harpnum: int | None = None) -> str:
    """Build a SHARP DRMS record-set query for a time interval."""

    harp_key = "" if harpnum is None else str(int(harpnum))
    return f"{series}[{harp_key}][{_to_jsoc_tai_time(start)}-{_to_jsoc_tai_time(end)}]"


def make_sharp_segment_recordset(
    series: str,
    harpnum: int,
    t_rec: str,
    *,
    segments: Sequence[str] = DEFAULT_SHARP_SEGMENTS,
) -> str:
    """Build a single-record SHARP export record-set with selected segments."""

    if not segments:
        raise ValueError("at least one SHARP segment is required")
    segment_text = ",".join(str(segment) for segment in segments)
    return f"{series}[{int(harpnum)}][{_to_jsoc_tai_time(t_rec)}]{{{segment_text}}}"


def _to_jsoc_tai_time(value: str) -> str:
    """Return a DRMS-friendly T_REC string."""

    text = str(value).strip()
    if text.endswith(("_TAI", "_UTC")):
        return text
    if "T" in text:
        date, time = text.split("T", 1)
    elif " " in text:
        date, time = text.split(" ", 1)
    else:
        return f"{text}_TAI" if "_" in text else text.replace("-", ".")
    return f"{date.replace('-', '.')}_{time}_TAI"


def _or_attrs(factory, values: Iterable):
    values = tuple(values)
    if not values:
        raise ValueError("at least one value is required")
    query = factory(values[0])
    for value in values[1:]:
        query |= factory(value)
    return query


def _with_sample(query, attrs, units, sample_minutes: float | None):
    if sample_minutes is None:
        return query
    return query & attrs.Sample(float(sample_minutes) * units.min)


def _sample_query_response(result, sample_minutes: float | None):
    """Greedily thin a Fido response when a provider does not support Sample."""

    if sample_minutes is None:
        return result

    import astropy.units as u
    from astropy.time import Time, TimeDelta

    cadence = TimeDelta(float(sample_minutes) * u.min)
    tolerance = TimeDelta(1 * u.s)
    sampled_tables = []
    for table in result:
        if len(table) == 0 or "Start time" not in table.colnames:
            sampled_tables.append(table)
            continue

        keep = []
        next_time = None
        for index, start_time in enumerate(table["Start time"]):
            current = Time(start_time)
            if next_time is None or current >= next_time - tolerance:
                keep.append(index)
                next_time = current + cadence
        sampled_tables.append(table[keep])

    if len(sampled_tables) == 1:
        return sampled_tables[0]
    return sampled_tables


def _sample_sharp_query_result(result, sample_minutes: float | None):
    """Thin a DRMS SHARP query result to an approximate cadence."""

    if sample_minutes is None or len(result) == 0:
        return result

    import astropy.units as u
    from astropy.time import Time, TimeDelta

    cadence = TimeDelta(float(sample_minutes) * u.min)
    tolerance = TimeDelta(1 * u.s)
    keep = []
    next_time = None
    for index, row in result.iterrows():
        t_rec = str(row["T_REC"]).replace("_TAI", "").replace("_", "T").replace(".", "-")
        current = Time(t_rec, format="isot", scale="tai")
        if next_time is None or current >= next_time - tolerance:
            keep.append(index)
            next_time = current + cadence
    return result.loc[keep].reset_index(drop=True)


def _sunpy_imports():
    from sunpy.net import Fido, attrs as a
    import astropy.units as u

    return Fido, a, u

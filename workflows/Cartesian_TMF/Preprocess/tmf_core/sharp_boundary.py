"""Read cropped SHARP vector-magnetogram components for local simulations."""

from __future__ import annotations

import warnings
from datetime import datetime
from pathlib import Path

import sunpy.map


def read_sharp_submap(
    bp_path: str | Path,
    bt_path: str | Path,
    br_path: str | Path,
    x1: int,
    x2: int,
    y1: int,
    y2: int,
):
    """Return local ``Bx, By, Bz, UTC`` for a 0-based exclusive-end crop.

    The Cartesian convention used throughout the TMF/NLFFF workflow is
    ``Bx=Bp``, ``By=-Bt``, and ``Bz=Br``.
    """

    warnings.filterwarnings(
        "ignore",
        category=UserWarning,
        module="sunpy.map.mapbase",
        message="Could not parse unit string.*",
    )
    br_map = sunpy.map.Map(br_path)
    bt_map = sunpy.map.Map(bt_path)
    bp_map = sunpy.map.Map(bp_path)
    bx = bp_map.data[y1:y2, x1:x2]
    by = -bt_map.data[y1:y2, x1:x2]
    bz = br_map.data[y1:y2, x1:x2]
    observed_at = _observation_time(br_map.meta)
    return bx, by, bz, observed_at


def _observation_time(metadata) -> datetime:
    for key in ("t_obs", "t_rec"):
        value = metadata.get(key)
        if value:
            return datetime.strptime(str(value)[:19], "%Y.%m.%d_%H:%M:%S")
    value = metadata.get("date-obs")
    if value:
        return datetime.fromisoformat(str(value)[:19])
    raise KeyError("SHARP FITS header has none of T_OBS, T_REC, or DATE-OBS")


# Public alias used by the Cartesian TMF driven-boundary workflow.
sub_data = read_sharp_submap

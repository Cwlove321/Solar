"""Reusable AMRVAC simulation-analysis helpers."""

from __future__ import annotations

import os
import sys
from pathlib import Path


def _add_local_simesh() -> None:
    """Make the repository's bundled simesh importable without installation."""
    candidates = []
    if os.environ.get("SIMESH_SRC"):
        candidates.append(Path(os.environ["SIMESH_SRC"]).expanduser())
    candidates.append(Path(__file__).resolve().parents[3] / "Packages/simesh/src")
    for candidate in candidates:
        if candidate.is_dir() and str(candidate) not in sys.path:
            sys.path.insert(0, str(candidate))
            break


_add_local_simesh()

from .config import load_config
from .snapshot import extract_axis_slice, inspect_snapshot

__all__ = ["extract_axis_slice", "inspect_snapshot", "load_config"]

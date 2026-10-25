"""Physical coordinate reconstruction for MPI-AMRVAC block grids.

MPI-AMRVAC DAT snapshots deliberately contain topology and field data, not the
``stretch_dim`` mesh-list settings.  This module recreates the exact coordinate
rules used by MPI-AMRVAC when those settings are supplied by the caller.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import re
from typing import Mapping

import numpy as np


def _fortran_float(value: str) -> float:
    return float(value.strip().replace("D", "e").replace("d", "e"))


def _indexed_values(text: str, name: str, ndim: int, cast, default):
    values = [default] * ndim
    # MPI-AMRVAC .par uses e.g. stretch_dim(1) = 'uni'.
    pattern = rf"\b{re.escape(name)}\s*\(\s*(\d+)\s*\)\s*=\s*([^!,/\n]+)"
    for index, value in re.findall(pattern, text, flags=re.IGNORECASE):
        i = int(index) - 1
        if 0 <= i < ndim:
            values[i] = cast(value.strip().strip("'\""))
    return values


def grid_config_from_par(path: str | Path, ndim: int = 3) -> dict:
    """Read the mesh stretching entries relevant to a DAT coordinate rebuild."""
    text = Path(path).read_text()
    stretch_type = _indexed_values(text, "stretch_dim", ndim, lambda v: v.lower(), "none")
    q = _indexed_values(text, "qstretch_baselevel", ndim, _fortran_float, None)
    nstretched = _indexed_values(text, "nstretchedblocks_baselevel", ndim, int, 0)
    return {
        "stretch_type": stretch_type,
        "qstretch_baselevel": q,
        "nstretchedblocks_baselevel": nstretched,
        "source": str(path),
    }


@dataclass(frozen=True)
class StretchConfig:
    """Validated coordinate settings supplied beside a DAT snapshot."""

    stretch_type: tuple[str, ...]
    qstretch_baselevel: tuple[float | None, ...]
    nstretchedblocks_baselevel: tuple[int, ...]

    @property
    def is_stretched(self) -> bool:
        return any(kind != "none" for kind in self.stretch_type)

    @classmethod
    def from_user(cls, value, ndim: int) -> "StretchConfig":
        if value is None:
            return cls(("none",) * ndim, (None,) * ndim, (0,) * ndim)
        if isinstance(value, (str, Path)):
            value = grid_config_from_par(value, ndim)
        if not isinstance(value, Mapping):
            raise TypeError("grid_config must be a mapping or an MPI-AMRVAC .par path")

        types = value.get("stretch_type", value.get("stretch_dim", ["none"] * ndim))
        qvalues = value.get("qstretch_baselevel", [None] * ndim)
        nvalues = value.get("nstretchedblocks_baselevel", [0] * ndim)
        if len(types) != ndim or len(qvalues) != ndim or len(nvalues) != ndim:
            raise ValueError(f"grid_config entries must each contain ndim={ndim} values")
        types = tuple(str(v).lower() for v in types)
        if any(v not in {"none", "uni", "symm"} for v in types):
            raise ValueError("stretch_type entries must be 'none', 'uni', or 'symm'")
        qvalues = tuple(None if v is None else float(v) for v in qvalues)
        nvalues = tuple(int(v) for v in nvalues)
        for i, kind in enumerate(types):
            if kind != "none" and (qvalues[i] is None or qvalues[i] <= 1.0):
                raise ValueError(f"stretched dimension {i} needs finite qstretch_baselevel > 1")
            if kind == "symm" and (nvalues[i] <= 0 or nvalues[i] % 2):
                raise ValueError(f"symmetric stretched dimension {i} needs even nstretchedblocks_baselevel")
        return cls(types, qvalues, nvalues)


class AMRVACCoordinates:
    """True per-block nodes and centers for uniform or MPI-AMRVAC stretched grids."""

    def __init__(self, xmin, xmax, domain_nx, block_nx, block_levels, block_ids, config: StretchConfig):
        self.xmin = np.asarray(xmin, dtype=float)
        self.xmax = np.asarray(xmax, dtype=float)
        self.domain_nx = np.asarray(domain_nx, dtype=int)
        self.block_nx = np.asarray(block_nx, dtype=int)
        self.block_levels = np.asarray(block_levels, dtype=int)
        self.block_ids = np.asarray(block_ids, dtype=int) - 1  # DAT tree indices are Fortran 1-based.
        self.ndim = len(self.domain_nx)
        self.config = config
        if self.block_ids.shape != (len(self.block_levels), self.ndim):
            raise ValueError("DAT tree block indices do not match coordinate dimensionality")

    @property
    def is_stretched(self) -> bool:
        return self.config.is_stretched

    def _level_nodes(self, level: int, dim: int) -> np.ndarray:
        """All cell nodes on one AMR level, matching MPI-AMRVAC's recursion."""
        n = int(self.domain_nx[dim]) * 2 ** (int(level) - 1)
        lo, hi = self.xmin[dim], self.xmax[dim]
        kind = self.config.stretch_type[dim]
        if kind == "none":
            return np.linspace(lo, hi, n + 1)

        q0 = self.config.qstretch_baselevel[dim]
        q = q0 ** (1.0 / 2 ** (int(level) - 1))
        # AMRVAC: dxfirst(1)=(xmax-xmin)(1-q)/(1-q**domain_nx), and each
        # refinement has dxfirst_l=dxfirst_(l-1)/(1+sqrt(q_(l-1))).
        if kind == "uni":
            first = (hi - lo) * (1.0 - q0) / (1.0 - q0 ** int(self.domain_nx[dim]))
            for _ in range(1, int(level)):
                first /= 1.0 + q0 ** (1.0 / 2 ** _)
            indices = np.arange(n + 1, dtype=float)
            return lo + first * np.expm1(indices * np.log(q)) / np.expm1(np.log(q))

        nstretch = self.config.nstretchedblocks_baselevel[dim] * 2 ** (int(level) - 1)
        edge = nstretch * int(self.block_nx[dim]) // 2
        if 2 * edge > n:
            raise ValueError("nstretchedblocks_baselevel exceeds the level-one domain")
        first1 = self._symmetric_first_cell(dim)
        first = first1
        for _ in range(1, int(level)):
            first /= 1.0 + q0 ** (1.0 / 2 ** _)
        middle_dx = first1 / 2 ** (int(level) - 1)
        widths = np.empty(n, dtype=float)
        widths[:edge] = first * q ** np.arange(edge)
        widths[edge:n-edge] = middle_dx
        widths[n-edge:] = widths[:edge][::-1]
        nodes = lo + np.concatenate(([0.0], np.cumsum(widths)))
        nodes[-1] = hi  # Remove harmless floating point accumulation drift.
        return nodes

    def _symmetric_first_cell(self, dim: int) -> float:
        n = int(self.domain_nx[dim])
        b = int(self.block_nx[dim])
        stretched_blocks = self.config.nstretchedblocks_baselevel[dim]
        m = stretched_blocks * b // 2
        q = self.config.qstretch_baselevel[dim]
        length = self.xmax[dim] - self.xmin[dim]
        xstretch = length / (2.0 + (n - stretched_blocks * b) * (1.0 - q) / (1.0 - q ** m))
        return xstretch * (1.0 - q) / (1.0 - q ** m)

    def block_nodes(self) -> tuple[np.ndarray, ...]:
        result = []
        for dim in range(self.ndim):
            out = np.empty((len(self.block_levels), int(self.block_nx[dim]) + 1), dtype=float)
            cache = {}
            for leaf, level in enumerate(self.block_levels):
                nodes = cache.setdefault((int(level), dim), self._level_nodes(int(level), dim))
                start = int(self.block_ids[leaf, dim]) * int(self.block_nx[dim])
                out[leaf] = nodes[start:start + int(self.block_nx[dim]) + 1]
            result.append(out)
        return tuple(result)

    def block_centers(self) -> tuple[np.ndarray, ...]:
        return tuple(0.5 * (nodes[:, 1:] + nodes[:, :-1]) for nodes in self.block_nodes())

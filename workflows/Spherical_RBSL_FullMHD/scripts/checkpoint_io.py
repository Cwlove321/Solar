"""Strict streaming validation of MPI-AMRVAC v5 restart snapshots.

Layout follows AMRVAC src/io/mod_input_output.t and simesh legacy datio.py:
native-endian fixed-width header; forest, levels, Morton indices, int64 block
offsets; per-block ghost extents, variable-major CC doubles and CT face doubles.
Only the 3-D MHD/internal-energy variant used by this workflow is accepted.
"""
from __future__ import annotations

import hashlib
import math
import struct
from pathlib import Path

import numpy as np


class CheckpointError(ValueError):
    pass


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(4 * 1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read(stream, fmt):
    size = struct.calcsize("=" + fmt)
    raw = stream.read(size)
    if len(raw) != size:
        raise CheckpointError(f"truncated header at byte {stream.tell()}")
    return struct.unpack("=" + fmt, raw)


def _name(stream):
    try:
        return _read(stream, "16s")[0].decode("ascii").strip(" \x00")
    except UnicodeDecodeError as exc:
        raise CheckpointError("non-ASCII DAT field name") from exc


def read_header(path: Path) -> dict:
    path = Path(path)
    with path.open("rb") as stream:
        version, = _read(stream, "i")
        if version != 5:
            raise CheckpointError(f"only AMRVAC v5 DAT is supported, got {version}")
        values = _read(stream, "9id")
        header = dict(zip(("offset_tree", "offset_blocks", "nw", "ndir", "ndim", "levmax", "nleafs", "nparents", "iteration", "time"), values))
        if header["ndim"] != 3 or header["ndir"] != 3:
            raise CheckpointError("this workflow requires 3-D/3-component DAT")
        if not 1 <= header["nw"] <= 1024 or not 1 <= header["levmax"] <= 30:
            raise CheckpointError("invalid variable count or maximum level")
        if not 0 < header["nleafs"] < 100_000_000 or header["nparents"] < 0:
            raise CheckpointError("invalid forest size")
        for key, fmt in (("xmin", "3d"), ("xmax", "3d"), ("domain_nx", "3i"), ("block_nx", "3i"), ("periodic", "3i")):
            header[key] = list(_read(stream, fmt))
        header["geometry"] = _name(stream)
        staggered, = _read(stream, "i")
        if staggered not in (0, 1):
            raise CheckpointError("invalid staggered logical")
        header["staggered"] = bool(staggered)
        header["w_names"] = [_name(stream).lower() for _ in range(header["nw"])]
        if len(set(header["w_names"])) != header["nw"]:
            raise CheckpointError("duplicate state variable names")
        header["physics_type"] = _name(stream)
        npars, = _read(stream, "i")
        if not 0 <= npars <= 1024:
            raise CheckpointError("invalid physics parameter count")
        header["parameters"] = list(_read(stream, f"{npars}d"))
        header["parameter_names"] = [_name(stream) for _ in range(npars)]
        header["snapshotnext"], header["slicenext"], header["collapsenext"] = _read(stream, "3i")
        header["header_end"] = stream.tell()
    scalars = [header["time"], *header["xmin"], *header["xmax"], *header["parameters"]]
    if not all(math.isfinite(v) for v in scalars) or header["time"] < 0 or header["iteration"] < 0:
        raise CheckpointError("non-finite or negative time/iteration in DAT header")
    if not all(a < b for a, b in zip(header["xmin"], header["xmax"])):
        raise CheckpointError("invalid DAT domain bounds")
    if any(n <= 0 for n in header["block_nx"] + header["domain_nx"]):
        raise CheckpointError("invalid DAT mesh dimensions")
    if any(d % b for d, b in zip(header["domain_nx"], header["block_nx"])):
        raise CheckpointError("DAT base mesh is not divisible by block size")
    if header["geometry"].lower() != "spherical_3d":
        raise CheckpointError("this workflow requires spherical_3D DAT geometry")
    if not header["staggered"]:
        raise CheckpointError("this workflow requires CT staggered face data")
    required = {"rho", "m1", "m2", "m3", "e", "b1", "b2", "b3"}
    if header["physics_type"].lower() != "mhd" or not required <= set(header["w_names"]):
        raise CheckpointError("complete MHD state rho,m1,m2,m3,e,b1,b2,b3 is required")
    if not header["header_end"] <= header["offset_tree"] < header["offset_blocks"] <= path.stat().st_size:
        raise CheckpointError("DAT header/tree/block offsets are out of bounds")
    if header["offset_tree"] + 28 * header["nleafs"] + 4 * header["nparents"] != header["offset_blocks"]:
        raise CheckpointError("DAT tree size does not match forest/leaf counts")
    header["version"] = version
    return header


def _array(stream, dtype, count):
    data = np.fromfile(stream, dtype=dtype, count=count)
    if len(data) != count:
        raise CheckpointError("truncated DAT tree")
    return data


def _check_doubles(stream, count, label, positive=False):
    smallest = math.inf
    largest = -math.inf
    consumed = 0
    while consumed < count:
        data = _array(stream, "=f8", min(count - consumed, 131072))
        mask = ~np.isfinite(data)
        if positive:
            mask |= data <= 0
        if np.any(mask):
            index = consumed + int(np.flatnonzero(mask)[0])
            value = float(data[int(np.flatnonzero(mask)[0])])
            raise CheckpointError(f"{label}: invalid value {value} at stored-cell index {index}")
        smallest = min(smallest, float(data.min()))
        largest = max(largest, float(data.max()))
        consumed += len(data)
    return smallest, largest


def validate_checkpoint(path: Path, *, target_time: float | None = None, expected_mesh: dict | None = None) -> dict:
    path = Path(path)
    h = read_header(path)
    file_size = path.stat().st_size
    if target_time is not None and h["time"] + max(1e-12, abs(target_time) * 1e-9) < target_time:
        raise CheckpointError(f"time {h['time']:.17g} has not reached target {target_time:.17g}")
    if expected_mesh:
        for key in ("domain_nx", "block_nx"):
            if key in expected_mesh and list(expected_mesh[key]) != h[key]:
                raise CheckpointError(f"DAT {key} differs from stage namelist")
        for key in ("xmin", "xmax"):
            if key in expected_mesh and not np.allclose(h[key], expected_mesh[key], rtol=1e-12, atol=1e-12):
                raise CheckpointError(f"DAT {key} differs from stage namelist")
    n = h["nleafs"]
    extents = np.array(h["domain_nx"]) // np.array(h["block_nx"])
    summary = {name: [math.inf, -math.inf] for name in h["w_names"]}
    ct_range = [math.inf, -math.inf]
    with path.open("rb") as stream:
        stream.seek(h["offset_tree"])
        forest = _array(stream, "=i4", n + h["nparents"])
        if np.any((forest != 0) & (forest != 1)) or np.count_nonzero(forest) != n:
            raise CheckpointError("invalid DAT forest flags")
        pending_levels = [1] * int(np.prod(extents))
        expected_levels = []
        for leaf in forest:
            if not pending_levels:
                raise CheckpointError("DAT forest has extra nodes")
            level = pending_levels.pop()
            if leaf:
                expected_levels.append(level)
            else:
                if level >= h["levmax"]:
                    raise CheckpointError("DAT forest refines beyond maximum level")
                pending_levels.extend([level + 1] * 8)
        if pending_levels:
            raise CheckpointError("DAT forest is incomplete")
        levels = _array(stream, "=i4", n)
        indices = _array(stream, "=i4", 3 * n).reshape(n, 3)
        offsets = _array(stream, "=i8", n)
        if np.any(levels < 1) or np.any(levels > h["levmax"]):
            raise CheckpointError("invalid DAT leaf level")
        if not np.array_equal(levels, expected_levels):
            raise CheckpointError("DAT leaf levels disagree with its forest topology")
        if np.any(indices < 1) or np.any(indices > extents * (2 ** (levels - 1))[:, None]):
            raise CheckpointError("invalid DAT leaf index")
        if len({(int(level), *map(int, index)) for level, index in zip(levels, indices)}) != n:
            raise CheckpointError("duplicate DAT leaf")
        if stream.tell() > h["offset_blocks"] or len(set(map(int, offsets))) != n:
            raise CheckpointError("overlapping DAT tree or duplicate block offsets")
        order = np.argsort(offsets)
        cursor = h["offset_blocks"]
        for leaf in order:
            offset = int(offsets[leaf])
            if offset != cursor:
                raise CheckpointError(f"DAT block gap/overlap at {offset}, expected {cursor}")
            stream.seek(offset)
            ghosts = np.array(_read(stream, "6i")).reshape(2, 3)
            if np.any(ghosts < 0) or np.any(ghosts > 64):
                raise CheckpointError("invalid DAT ghost-cell extents")
            shape = np.array(h["block_nx"]) + ghosts.sum(axis=0)
            cells = int(np.prod(shape))
            end = offset + 24 + 8 * (h["nw"] * cells + (3 * int(np.prod(shape + 1)) if h["staggered"] else 0))
            if end > file_size:
                raise CheckpointError("truncated DAT block")
            for name in h["w_names"]:
                lo, hi = _check_doubles(stream, cells, f"leaf {leaf + 1}, {name}", name in {"rho", "e"})
                summary[name][0] = min(summary[name][0], lo)
                summary[name][1] = max(summary[name][1], hi)
            if h["staggered"]:
                lo, hi = _check_doubles(stream, 3 * int(np.prod(shape + 1)), f"leaf {leaf + 1}, CT faces")
                ct_range = [min(ct_range[0], lo), max(ct_range[1], hi)]
            cursor = end
        if cursor != file_size:
            raise CheckpointError(f"DAT has {file_size - cursor} unrecognized trailing bytes")
    return {"path": str(path.resolve()), "sha256": sha256_file(path), "bytes": path.stat().st_size,
            "time": h["time"], "iteration": h["iteration"], "nleafs": n,
            "levels": {str(level): int(np.count_nonzero(levels == level)) for level in np.unique(levels)},
            "ranges": summary, "ct_range": ct_range if h["staggered"] else None,
            "validated": True, "format": "AMRVAC-v5-3D-MHD-internal-e"}


def select_checkpoint(paths, target_time: float, expected_mesh: dict | None = None) -> dict:
    """Select the latest actual code time, never a hard-coded snapshot number."""
    candidates = []
    for path in paths:
        header = read_header(Path(path))
        candidates.append((header["time"], header["iteration"], str(path), Path(path)))
    if not candidates:
        raise CheckpointError("stage produced no DAT snapshots")
    # Validate every produced DAT; an earlier nonfinite snapshot must not be
    # hidden by a later positive-state repair. Check newest first to fail fast
    # on incomplete final writes, never falling back to an older checkpoint.
    ordered = sorted(candidates, reverse=True)
    result = validate_checkpoint(ordered[0][3], target_time=target_time, expected_mesh=expected_mesh)
    reports = [result]
    reports += [validate_checkpoint(item[3], expected_mesh=expected_mesh) for item in ordered[1:]]
    result["validated_snapshots"] = [{key: report[key] for key in ("path", "sha256", "bytes", "time", "iteration")} for report in reports]
    return result

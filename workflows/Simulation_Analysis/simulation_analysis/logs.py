"""Parsing and compact reporting for AMRVAC per-step logs."""

from __future__ import annotations

import csv
import json
import math
import re
import subprocess
import time
from pathlib import Path


def fortran_float(value: str) -> float:
    value = value.strip().replace("D", "E").replace("d", "e")
    if "e" not in value.lower():
        value = re.sub(r"(?<=\d)([+-]\d{2,4})$", r"E\1", value)
    return float(value)


def read_rows(config: dict) -> list[dict[str, float | int]]:
    run = Path(config["run_dir"])
    paths = sorted(run.glob(config["log"]["pattern"]))
    if not paths:
        return []
    columns = list(config["log"]["columns"])
    rows = []
    for line in paths[-1].read_text(errors="replace").splitlines():
        values = line.split("|", 1)[0].split()
        if len(values) != len(columns) or not values[0].isdigit():
            continue
        try:
            row: dict[str, float | int] = {columns[0]: int(values[0])}
            for name, value in zip(columns[1:], values[1:]):
                number = fortran_float(value)
                row[name] = int(number) if name.startswith("n") and number.is_integer() else number
            momentum_names = [name for name in ("m1_extreme", "m2_extreme", "m3_extreme") if name in row]
            if momentum_names:
                row["momentum_extreme"] = max(abs(float(row[name])) for name in momentum_names)
            rows.append(row)
        except (ValueError, OverflowError):
            continue
    return rows


def slope(rows: list[dict[str, float | int]], field: str) -> float | None:
    finite = [
        (float(row["time"]), float(row[field]))
        for row in rows
        if field in row and math.isfinite(float(row[field]))
    ]
    if len(finite) < 2 or finite[-1][0] == finite[0][0]:
        return None
    return (finite[-1][1] - finite[0][1]) / (finite[-1][0] - finite[0][0])


def job_status(job_id: str | None) -> str | None:
    if not job_id:
        return None
    result = subprocess.run(
        ["bjobs", "-a", "-noheader", "-o", "stat", str(job_id)],
        capture_output=True,
        text=True,
        check=False,
    )
    return result.stdout.strip() or "UNKNOWN"


def build_report(config: dict, job_id: str | None = None) -> dict:
    run = Path(config["run_dir"])
    rows = read_rows(config)
    dat_files = sorted(run.glob(config["data_glob"]))
    errors = []
    for path in sorted((run / "log").glob("*.err")):
        if path.stat().st_size:
            lines = path.read_text(errors="replace").splitlines()
            errors.append({"file": path.name, "bytes": path.stat().st_size, "tail": lines[-20:]})
    finite_dt = [float(row["dt"]) for row in rows if "dt" in row and math.isfinite(float(row["dt"]))]
    nonfinite = None
    for row in rows:
        if any(isinstance(value, float) and not math.isfinite(value) for value in row.values()):
            nonfinite = int(row["it"])
            break
    tail = rows[-100:]
    return {
        "updated_unix": time.time(),
        "run": str(run),
        "job_id": job_id,
        "job_status": job_status(job_id),
        "rows": len(rows),
        "last": rows[-1] if rows else None,
        "dat_files": [{"name": path.name, "bytes": path.stat().st_size} for path in dat_files],
        "min_positive_dt": min((value for value in finite_dt if value > 0.0), default=None),
        "first_nonfinite_it": nonfinite,
        "stderr": errors,
        "last_100_slopes_per_code_time": {
            field: slope(tail, field)
            for field in ("dt", "rho_max", "momentum_extreme", "e_max")
            if any(field in row for row in tail)
        },
    }


def write_report(config: dict, report: dict) -> tuple[Path, Path]:
    output = Path(config["output_dir"])
    output.mkdir(parents=True, exist_ok=True)
    json_path = output / "run_monitor.json"
    temp = json_path.with_suffix(".json.tmp")
    temp.write_text(json.dumps(report, indent=2, allow_nan=True) + "\n", encoding="ascii")
    temp.replace(json_path)

    csv_path = output / "run_monitor.csv"
    rows = read_rows(config)
    if rows:
        temp_csv = csv_path.with_suffix(".csv.tmp")
        with temp_csv.open("w", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        temp_csv.replace(csv_path)
    return json_path, csv_path

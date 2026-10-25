#!/usr/bin/env python3
"""Populate inputs and synchronize event parameters across all stages."""

from __future__ import annotations

import argparse
import json
import math
import re
import shutil
from pathlib import Path
from typing import Any

from atmosphere_config import FORTRAN_KEYS, normalize_atmosphere


ROOT = Path(__file__).resolve().parents[1]
STAGES = ("01_pfss", "02_relax", "03_insert", "04_fullmhd")
GLOBAL_STAGE = "01_pfss_global_ct"


def fortran_value(value: Any) -> str:
    if isinstance(value, bool):
        return ".true." if value else ".false."
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError("namelist numbers must be finite")
        return f"{value:.16g}d0" if "e" not in f"{value:.16g}".lower() else f"{value:.16e}".replace("e", "d")
    if isinstance(value, str) and re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", value):
        return f"'{value}'"
    raise TypeError(f"unsupported namelist value: {value!r}")


def replace_key_text(text: str, path: Path, key: str, value: Any) -> str:
    pattern = re.compile(rf"(?m)^(\s*{re.escape(key)}\s*=\s*).*$")
    replacement = rf"\g<1>{fortran_value(value)}"
    updated, count = pattern.subn(replacement, text)
    if count != 1:
        raise ValueError(f"{path}: expected exactly one {key}= line, found {count}")
    return updated


def configure(config_path: Path, atmosphere_only: bool = False) -> None:
    config = json.loads(config_path.read_text(encoding="utf-8"))
    if not isinstance(config, dict):
        raise ValueError("case configuration must be a JSON object")
    atmosphere = normalize_atmosphere(config["atmosphere"]) if "atmosphere" in config else normalize_atmosphere()
    par_files = [ROOT / stage / "amrvac.par" for stage in STAGES]
    par_files.extend(sorted((ROOT / "04_fullmhd" / "profiles").glob("*.par")))
    atmosphere_files = list(par_files)
    if (ROOT / GLOBAL_STAGE).is_dir():
        atmosphere_files.append(ROOT / GLOBAL_STAGE / "amrvac.par")
    updates: dict[Path, str] = {}

    def queue_replace(path: Path, key: str, value: Any) -> None:
        if path not in updates:
            updates[path] = path.read_text(encoding="utf-8")
        updates[path] = replace_key_text(updates[path], path, key, value)

    # Prepare all atmosphere replacements before the first filesystem write.
    # The optional global PFSS stage receives atmosphere only, not local grids.
    for path in atmosphere_files:
        for key, fortran_key in FORTRAN_KEYS.items():
            queue_replace(path, fortran_key, atmosphere[key])
    if atmosphere_only:
        for path, updated in updates.items():
            path.write_text(updated, encoding="utf-8")
        print(f"configured atmosphere only in {len(atmosphere_files)} namelists")
        print(json.dumps(atmosphere, indent=2, ensure_ascii=False))
        return

    boundary_dir = Path(config["boundary_dir"]).expanduser().resolve()
    path_dat = Path(config["path_dat"]).expanduser().resolve()
    brll = boundary_dir / "Brll0000.fits"
    for required in (brll, path_dat):
        if not required.is_file():
            raise FileNotFoundError(required)

    npoints = int(config["rbsl"]["npoints"])
    expected_path_bytes = npoints * 3 * 8
    if path_dat.stat().st_size != expected_path_bytes:
        raise ValueError(
            f"{path_dat}: binary path.dat must contain {npoints} interleaved "
            f"float64 triples ({expected_path_bytes} bytes)"
        )

    grid = config["grid"]
    grid_keys = (
        "domain_nx1",
        "domain_nx2",
        "domain_nx3",
        "xprobmin1",
        "xprobmax1",
        "xprobmin2",
        "xprobmax2",
        "xprobmin3",
        "xprobmax3",
        "refine_max_level",
    )
    for path in par_files:
        for key in grid_keys:
            queue_replace(path, key, grid[key])
        queue_replace(path, "qstretch_baselevel(1)", grid["qstretch_baselevel1"])

    rbsl = config["rbsl"]
    rbsl_values = {
        "Bflux": rbsl["bflux_mx"],
        "Rmfr": rbsl["radius_mm"],
        "I_Helix": rbsl["positive_helicity"],
        "npoints": rbsl["npoints"],
        "flip_axis": rbsl["flip_axis"],
    }
    for path in par_files:
        for key, value in rbsl_values.items():
            queue_replace(path, key, value)

    filament = config["filament"]
    method = str(filament.get("method", "instant")).lower()
    if method not in {"instant", "liakh_source", "none"}:
        raise ValueError("filament.method must be instant, liakh_source, or none")
    instant_loading = method == "instant"
    liakh = filament.get("liakh_mass_loading", {})
    liakh_enabled = method == "liakh_source"

    liakh_values = {
        "liakh_mass_loading": liakh_enabled,
        "liakh_density_contrast": liakh.get("density_contrast", 30.0),
        "liakh_t_start_s": liakh.get("start_time_s", 0.0),
        "liakh_t_load_s": liakh.get("load_time_s", 100.0),
        "liakh_sigma_horizontal_Mm": liakh.get("sigma_horizontal_mm", 3.6),
        "liakh_sigma_radial_Mm": liakh.get("sigma_radial_mm", 4.8),
        "liakh_min_height_Mm": liakh.get("minimum_height_mm", 5.0),
    }
    insert_par = ROOT / "03_insert" / "amrvac.par"
    queue_replace(insert_par, "load_filament_density", instant_loading)
    queue_replace(insert_par, "filament_density_factor", filament["density_factor"])
    queue_replace(insert_par, "filament_core_fraction", filament["core_fraction"])
    queue_replace(insert_par, "filament_min_height_Mm", filament["minimum_height_mm"])
    for key, value in liakh_values.items():
        queue_replace(insert_par, key, value)

    insert_stage = config.get("insert_stage", {})
    queue_replace(insert_par, "time_max", insert_stage.get("time_max_code", 0.0))
    queue_replace(insert_par, "dtsave_dat", insert_stage.get("save_interval_code", 0.01))
    queue_replace(insert_par, "autoconvert", insert_stage.get("autoconvert", False))

    resolved = dict(config)
    resolved["boundary_dir"] = str(boundary_dir)
    resolved["path_dat"] = str(path_dat)
    resolved["atmosphere"] = atmosphere
    resolved_text = json.dumps(resolved, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    # Preflight succeeded: only now copy inputs and commit prepared parameters.
    (ROOT / "01_pfss" / "Bphot").mkdir(parents=True, exist_ok=True)
    shutil.copy2(brll, ROOT / "01_pfss" / "Bphot" / "Brll0000.fits")
    shutil.copy2(path_dat, ROOT / "01_pfss" / "path.dat")
    shutil.copy2(path_dat, ROOT / "03_insert" / "path.dat")
    for path, updated in updates.items():
        path.write_text(updated, encoding="utf-8")
    (ROOT / "case.resolved.json").write_text(
        resolved_text, encoding="utf-8"
    )
    print("configured all stages")
    print(f"Br input: {brll}")
    print(f"RBSL path: {path_dat}")
    print(f"Filament loading: {method}")
    print(f"03_insert time_max: {insert_stage.get('time_max_code', 0.0)}")
    print(f"resolved config: {ROOT / 'case.resolved.json'}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("config", type=Path, help="Case JSON based on case.example.json")
    parser.add_argument("--atmosphere-only", action="store_true", help="Synchronize only atmosphere settings; boundary/path inputs are not required")
    args = parser.parse_args()
    configure(args.config.expanduser().resolve(), atmosphere_only=args.atmosphere_only)


if __name__ == "__main__":
    main()

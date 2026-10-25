#!/usr/bin/env python3
"""Static checks for stage consistency and common production mistakes."""

from __future__ import annotations

import hashlib
import re
import sys
from pathlib import Path

from atmosphere_config import atmosphere_from_namelist


ROOT = Path(__file__).resolve().parents[1]
STAGES = ("01_pfss", "02_relax", "03_insert", "04_fullmhd")
GLOBAL_STAGE = "01_pfss_global_ct"


def source_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_assignments(path: Path) -> dict[str, str]:
    """Parse the supported, one-assignment-per-line template syntax strictly.

    A static checker must reject syntax it cannot interpret rather than miss a
    later Fortran namelist override. Arrays on one assignment line are allowed;
    groups, their closing '/', and assignments must occupy separate lines.
    """
    values: dict[str, str] = {}
    current = None
    seen_groups: set[str] = set()

    def visible_code(raw: str, number: int) -> tuple[str, str]:
        quote = None
        masked: list[str] = []
        index = 0
        while index < len(raw):
            char = raw[index]
            if quote:
                masked.append(" ")
                if char == quote:
                    if index + 1 < len(raw) and raw[index + 1] == quote:
                        masked.append(" ")
                        index += 1
                    else:
                        quote = None
            elif char in "\"'":
                quote = char
                masked.append(" ")
            elif char == "!":
                break
            else:
                masked.append(char)
            index += 1
        if quote:
            raise ValueError(f"{path}:{number}: unterminated string or multiline value is unsupported")
        return raw[:index].strip(), "".join(masked).strip()

    for number, raw in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        line, outside_quotes = visible_code(raw, number)
        if not line:
            continue
        if outside_quotes.startswith("$"):
            raise ValueError(f"{path}:{number}: $ namelist syntax is unsupported; use a standalone &group")
        if outside_quotes.startswith("&"):
            if current:
                raise ValueError(f"{path}:{number}: unclosed namelist {current} before another group")
            if not re.fullmatch(r"&[A-Za-z][A-Za-z0-9_]*", line):
                raise ValueError(f"{path}:{number}: group header must be standalone; inline assignments are unsupported")
            current = line[1:].strip().lower()
            if current in seen_groups:
                raise ValueError(f"{path}:{number}: duplicate namelist group {current}")
            seen_groups.add(current)
            continue
        if line == "/":
            if current is None:
                raise ValueError(f"{path}:{number}: closing / outside a namelist")
            current = None
            continue
        if current is None:
            raise ValueError(f"{path}:{number}: unexpected text outside a namelist")
        if outside_quotes.count("=") != 1 or any(char in outside_quotes for char in "&$/"):
            raise ValueError(f"{path}:{number}: expected one complete assignment; multiline/multiple assignments are unsupported")
        key, value = line.split("=", 1)
        key = key.strip().lower()
        if not re.fullmatch(r"[a-z][a-z0-9_]*(?:\(\s*\d+(?:\s*,\s*\d+)*\s*\))?", key):
            raise ValueError(f"{path}:{number}: unsupported assignment key {key!r}")
        value = value.strip().rstrip(",").strip()
        if not value:
            raise ValueError(f"{path}:{number}: empty or continued assignment is unsupported")
        normalized = f"{current}.{key}"
        if normalized in values:
            raise ValueError(f"{path}:{number}: duplicate key {normalized}")
        values[normalized] = value
    if current is not None:
        raise ValueError(f"{path}: unclosed namelist {current} at end of file")
    return values


def truth(value: str) -> bool:
    return value.lower() == ".true."


def require(condition: bool, message: str, errors: list[str]) -> None:
    if not condition:
        errors.append(message)


def normalized_subroutine(source: str, name: str) -> str:
    """Ignore comments/case/layout while checking the shared implementations."""
    uncommented = "\n".join(line.split("!", 1)[0] for line in source.splitlines())
    pattern = rf"(?ims)^\s*subroutine\s+{re.escape(name)}\b.*?^\s*end\s+subroutine\b[^\n]*"
    matches = re.findall(pattern, uncommented)
    if len(matches) != 1:
        raise ValueError(f"expected exactly one subroutine {name}, found {len(matches)}")
    return re.sub(r"\s+|&", "", matches[0]).lower()


def validate_atmosphere_sources(errors: list[str]) -> None:
    """Keep atmosphere physics identical without constraining unrelated code."""
    canonical: dict[str, str] = {}
    for stage in (*STAGES, GLOBAL_STAGE):
        source = (ROOT / stage / "mod_usr.t").read_text(encoding="utf-8")
        for name in ("inithdstatic", "read_atmosphere_params", "validate_atmosphere_parameters"):
            try:
                body = normalized_subroutine(source, name)
            except ValueError as exc:
                errors.append(f"{stage}: {exc}")
                continue
            if name not in canonical:
                canonical[name] = body
            require(body == canonical[name], f"{stage}: {name} differs from shared atmosphere implementation", errors)
            if name == "inithdstatic":
                require(
                    re.search(r"/\(pa\([^()]*\)-pa\([^()]*\)\)", body) is None,
                    f"{stage}: pressure-difference interpolation denominator has regressed",
                    errors,
                )


def validate_refinement_sources(errors: list[str]) -> None:
    """All five stages use the same bottom-boundary-only user refinement."""
    reference = normalized_subroutine(
        (ROOT / STAGES[0] / "mod_usr.t").read_text(encoding="utf-8"),
        "special_refine_grid",
    )
    for stage in (*STAGES[1:], GLOBAL_STAGE):
        try:
            actual = normalized_subroutine(
                (ROOT / stage / "mod_usr.t").read_text(encoding="utf-8"),
                "special_refine_grid",
            )
        except ValueError as exc:
            errors.append(f"{stage}: {exc}")
            continue
        require(
            actual == reference,
            f"{stage}: special_refine_grid differs from the shared bottom-boundary-only policy",
            errors,
        )


def main() -> int:
    errors: list[str] = []
    warnings: list[str] = []
    hashes = {stage: source_hash(ROOT / stage / "mod_usr.t") for stage in STAGES}
    background_hashes = {stage: hashes[stage] for stage in ("01_pfss", "02_relax", "04_fullmhd")}
    require(
        len(set(background_hashes.values())) == 1,
        f"01_pfss/02_relax/04_fullmhd mod_usr.t copies differ: {background_hashes}",
        errors,
    )
    shared_source = (ROOT / STAGES[0] / "mod_usr.t").read_text(encoding="utf-8")
    for required_text in (
        "access='stream',form='unformatted',convert='little_endian'",
        "expected_path_bytes=24_8*int(np,kind=8)",
    ):
        require(
            required_text in shared_source,
            f"shared mod_usr.t: missing binary path invariant: {required_text}",
            errors,
        )
    insert_source = (ROOT / "03_insert" / "mod_usr.t").read_text(encoding="utf-8")
    for required_text in (
        "liakh_mass_loading",
        "call prepare_liakh_mass_source_geometry()",
        "call add_liakh_filament_mass_source",
        "rho_source=liakh_density_contrast*rho_background/t_load*weight",
    ):
        require(
            required_text in insert_source,
            f"03_insert/mod_usr.t: missing Liakh source invariant: {required_text}",
            errors,
        )

    pars = {stage: parse_assignments(ROOT / stage / "amrvac.par") for stage in STAGES}
    profiles = {
        path.name: parse_assignments(path)
        for path in sorted((ROOT / "04_fullmhd" / "profiles").glob("*.par"))
    }
    global_par = parse_assignments(ROOT / GLOBAL_STAGE / "amrvac.par")
    atmosphere_pars = {**pars, GLOBAL_STAGE: global_par, **profiles}
    canonical_atmosphere = None
    for name, values in atmosphere_pars.items():
        try:
            atmosphere = atmosphere_from_namelist(values)
        except ValueError as exc:
            errors.append(f"{name}: {exc}")
            continue
        if canonical_atmosphere is None:
            canonical_atmosphere = atmosphere
        require(
            atmosphere == canonical_atmosphere,
            f"{name}: atmosphere differs from other stages/profiles; synchronize before starting or restarting",
            errors,
        )
    validate_atmosphere_sources(errors)
    validate_refinement_sources(errors)
    for stage, values in pars.items():
        require(values.get("meshlist.itfixgrid") == "0", f"{stage}: itfixgrid must be 0", errors)
        require(values.get("meshlist.ditregrid") == "3", f"{stage}: ditregrid must be 3", errors)
        require(values.get("meshlist.refine_criterion") == "3", f"{stage}: refine_criterion must be 3", errors)
        require(values.get("mhd_list.typedivbfix", "").strip("'\"") == "ct", f"{stage}: CT is required", errors)

        job = (ROOT / stage / "job").read_text(encoding="utf-8")
        rank_match = re.search(r"(?m)^#BSUB -n\s+(\d+)$", job)
        require(rank_match is not None, f"{stage}: missing #BSUB -n", errors)
        if rank_match:
            require(int(rank_match.group(1)) % 64 == 0, f"{stage}: ranks must be a multiple of 64", errors)
        require('#BSUB -R "span[ptile=64]"' in job, f"{stage}: missing span[ptile=64]", errors)

    require("filelist.restart_from_file" not in pars["01_pfss"], "01_pfss must not restart", errors)
    require(truth(pars["01_pfss"]["mag_list.do_rbsl"]), "01_pfss must subtract RBSL Br", errors)
    require(not truth(pars["02_relax"]["filelist.firstprocess"]), "02_relax must not insert RBSL", errors)
    require(truth(pars["03_insert"]["filelist.firstprocess"]), "03_insert must insert RBSL once at startup", errors)
    require(truth(pars["03_insert"]["mag_list.do_rbsl"]), "03_insert do_rbsl must be true", errors)
    require(pars["03_insert"].get("meshlist.iprob") == "2", "03_insert iprob must be 2", errors)
    require(not truth(pars["04_fullmhd"]["filelist.firstprocess"]), "04_fullmhd must be continuation-only", errors)
    require(not truth(pars["04_fullmhd"]["mag_list.do_rbsl"]), "04_fullmhd do_rbsl must be false", errors)
    instant_loading = truth(pars["03_insert"]["mag_list.load_filament_density"])
    liakh_loading = truth(pars["03_insert"]["mag_list.liakh_mass_loading"])
    require(
        not (instant_loading and liakh_loading),
        "03_insert instant density and Liakh source cannot both be enabled",
        errors,
    )
    if not instant_loading and not liakh_loading:
        warnings.append("both filament-loading methods are disabled")
    if not (ROOT / "03_insert" / "path.dat").is_file():
        warnings.append("03_insert/path.dat is absent; run configure_case.py before building")
    if liakh_loading:
        require(
            float(pars["03_insert"]["mag_list.liakh_density_contrast"].lower().replace("d", "e")) > 0.0,
            "Liakh density contrast must be positive",
            errors,
        )
        require(
            float(pars["03_insert"]["mag_list.liakh_t_start_s"].lower().replace("d", "e")) >= 0.0,
            "Liakh start time must be non-negative",
            errors,
        )
        for key in ("liakh_t_load_s", "liakh_sigma_horizontal_mm", "liakh_sigma_radial_mm"):
            require(
                float(pars["03_insert"][f"mag_list.{key}"].lower().replace("d", "e")) > 0.0,
                f"03_insert {key} must be positive",
                errors,
            )
    active_radiation = truth(pars["04_fullmhd"]["mhd_list.mhd_radiative_cooling"])
    active_heating = float(
        pars["04_fullmhd"]["mag_list.background_heating_cgs"].lower().replace("d", "e")
    )
    if not active_radiation and active_heating > 0.0:
        warnings.append("04_fullmhd has background heating but no radiative sink; expect secular heating")
    if active_radiation and active_heating == 0.0:
        warnings.append("04_fullmhd has radiation but no background heating; check thermal balance")

    require(
        "filelist.restart_from_file" not in global_par,
        f"{GLOBAL_STAGE} must not restart",
        errors,
    )
    require(
        truth(global_par.get("mag_list.do_rbsl", ".false.")),
        f"{GLOBAL_STAGE} must subtract RBSL Br",
        errors,
    )
    require(
        global_par.get("mhd_list.typedivbfix", "").strip("'\"") == "ct",
        f"{GLOBAL_STAGE}: CT is required",
        errors,
    )
    require(
        global_par.get("filelist.autoconvert", ".true.").lower() == ".false.",
        f"{GLOBAL_STAGE}: the current global-PFSS template disables online conversion",
        errors,
    )
    require(
        float(global_par["stoplist.time_max"].lower().replace("d", "e")) == 0.0,
        f"{GLOBAL_STAGE}: initial validation must use time_max=0",
        errors,
    )
    require(
        int(global_par["meshlist.max_blocks"]) >= 60000,
        f"{GLOBAL_STAGE}: max_blocks is below the retained template capacity (60000 per rank)",
        errors,
    )
    require(
        int(global_par["globalpfss_list.global_patch_x0"]) >= 0
        and int(global_par["globalpfss_list.global_patch_y0"]) >= 0,
        f"{GLOBAL_STAGE}: patch indices must be zero-based non-negative values",
        errors,
    )
    global_source = (ROOT / GLOBAL_STAGE / "mod_usr.t").read_text(encoding="utf-8")
    for required_text in (
        "use mod_pfss",
        "call pfss_vector_potential",
        "call prepare_global_boundary_with_rbsl_subtraction()",
        "br_rbsl(m+1-itheta,iphi)",
        "access='stream',form='unformatted',convert='little_endian'",
        "expected_path_bytes=24_8*int(np,kind=8)",
        "usr_before_main_loop      => initial_divb_diagnostics",
    ):
        require(
            required_text in global_source,
            f"{GLOBAL_STAGE}: missing source invariant: {required_text}",
            errors,
        )
    global_job = (ROOT / GLOBAL_STAGE / "job").read_text(encoding="utf-8")
    global_rank = re.search(r"(?m)^#BSUB -n\s+(\d+)$", global_job)
    require(global_rank is not None, f"{GLOBAL_STAGE}: missing #BSUB -n", errors)
    if global_rank:
        require(
            int(global_rank.group(1)) % 64 == 0,
            f"{GLOBAL_STAGE}: ranks must be a multiple of 64",
            errors,
        )
    require(
        '#BSUB -R "span[ptile=64]"' in global_job,
        f"{GLOBAL_STAGE}: missing span[ptile=64]",
        errors,
    )
    full_thermo = profiles.get("amrvac_full_thermo.par", {})
    require(
        truth(full_thermo.get("mhd_list.mhd_radiative_cooling", ".false.")),
        "full-thermo profile must enable radiative cooling",
        errors,
    )
    for name, values in profiles.items():
        require(
            not truth(values.get("filelist.firstprocess", ".true.")),
            f"{name}: firstprocess must be false for continuation",
            errors,
        )
        require(
            not truth(values.get("mag_list.do_rbsl", ".true.")),
            f"{name}: do_rbsl must be false",
            errors,
        )

    if errors:
        print("workflow validation failed:")
        for error in errors:
            print(f"  - {error}")
        return 1

    print("workflow validation passed")
    for warning in warnings:
        print(f"warning: {warning}")
    print(f"shared 01/02/04 mod_usr.t sha256: {next(iter(background_hashes.values()))}")
    print(f"03_insert mod_usr.t sha256: {hashes['03_insert']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

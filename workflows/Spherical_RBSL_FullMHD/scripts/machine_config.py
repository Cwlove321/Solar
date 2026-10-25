"""Machine-local build and launch settings; no compilation during configuration."""
from __future__ import annotations

import json
import os
import re
import shlex
from pathlib import Path

STAGES = ("01_pfss", "01_pfss_global_ct", "02_relax", "03_insert", "04_fullmhd")


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def load_machine(path: Path) -> dict:
    path = Path(path).expanduser().resolve()
    config = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(config, dict) or config.get("schema_version", 1) != 1:
        raise ValueError("machine configuration must use schema_version=1")
    defaults = {"schema_version": 1, "arch": "default", "make_jobs": 4,
                "module_setup": [], "environment": {}, "library_dirs": [],
                "libraries": ["pdfi_ss", "fishpack", "cfitsio"], "stage_libraries": {},
                "launcher": ["mpirun", "-np", "{ranks}"], "ranks": {},
                "scheduler": {"kind": "local"}, "compiler_commands": ["mpif90"],
                "executable": "./amrvac", "make_command": "make", "setup_command": None}
    unknown = set(config) - set(defaults) - {"amrvac_dir"}
    if unknown:
        raise ValueError(f"unknown machine settings: {sorted(unknown)}")
    result = defaults | config
    if "amrvac_dir" not in result:
        raise ValueError("machine.amrvac_dir is required")
    def resolve(value):
        value = Path(os.path.expandvars(str(value))).expanduser()
        return str((path.parent / value).resolve() if not value.is_absolute() else value.resolve())
    result["amrvac_dir"] = resolve(result["amrvac_dir"])
    if any(char.isspace() for char in result["amrvac_dir"]):
        raise ValueError("upstream AMRVAC build source path cannot contain whitespace; case directories may contain spaces")
    result["library_dirs"] = [resolve(v) for v in result["library_dirs"]]
    if not re.fullmatch(r"[A-Za-z0-9_.-]+", result["arch"]):
        raise ValueError("arch must be a Makefile architecture name")
    if not isinstance(result["make_jobs"], int) or isinstance(result["make_jobs"], bool) or result["make_jobs"] < 1:
        raise ValueError("make_jobs must be a positive integer")
    for key in ("module_setup", "library_dirs", "libraries", "launcher", "compiler_commands"):
        if not isinstance(result[key], list) or not all(isinstance(v, str) and "\n" not in v and "\x00" not in v for v in result[key]):
            raise ValueError(f"machine.{key} must be an array of single-line strings")
    if not result["launcher"]:
        raise ValueError("launcher must not be empty")
    if not isinstance(result["environment"], dict) or not all(re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", k) and isinstance(v, str) and "\x00" not in v for k, v in result["environment"].items()):
        raise ValueError("environment must map valid environment names to strings")
    if not isinstance(result["ranks"], dict) or set(result["ranks"]) - set(STAGES):
        raise ValueError("ranks keys must be workflow stage directory names")
    for stage in STAGES:
        ranks = result["ranks"].setdefault(stage, 1)
        if not isinstance(ranks, int) or isinstance(ranks, bool) or ranks <= 0:
            raise ValueError(f"ranks[{stage}] must be a positive integer")
    for key in ("libraries",):
        if not all(re.fullmatch(r"[A-Za-z0-9_+.-]+", v) for v in result[key]):
            raise ValueError(f"invalid library name in {key}")
    if not isinstance(result["stage_libraries"], dict) or set(result["stage_libraries"]) - set(STAGES):
        raise ValueError("stage_libraries keys must be workflow stages")
    for libraries in result["stage_libraries"].values():
        if not isinstance(libraries, list) or not all(isinstance(v, str) and re.fullmatch(r"[A-Za-z0-9_+.-]+", v) for v in libraries):
            raise ValueError("stage_libraries must contain arrays of library names")
    scheduler = result["scheduler"]
    if not isinstance(scheduler, dict) or scheduler.get("kind") not in {"local", "lsf"}:
        raise ValueError("scheduler.kind must be local or lsf (Slurm is not supported)")
    if set(scheduler) - {"kind", "queue", "ptile", "submit_command", "cancel_command"}:
        raise ValueError("unknown scheduler setting")
    if scheduler.get("queue") is not None and not re.fullmatch(r"[A-Za-z0-9_.-]+", scheduler["queue"]):
        raise ValueError("invalid LSF queue name")
    if scheduler.get("ptile") is not None and (not isinstance(scheduler["ptile"], int) or scheduler["ptile"] <= 0):
        raise ValueError("scheduler.ptile must be a positive integer")
    for token in result["launcher"]:
        try:
            token.format(ranks=1, executable="./amrvac", par="amrvac.par")
        except (KeyError, ValueError) as exc:
            raise ValueError(f"unsupported launcher placeholder: {token}") from exc
    if any(c in result["amrvac_dir"] for c in "\n\x00") or any(any(c in v for c in "\n\x00#$") for v in result["library_dirs"]):
        raise ValueError("invalid newline or Makefile metacharacter in source/library path")
    return result


def shell_environment(machine: dict) -> list[str]:
    lines = list(machine["module_setup"])
    values = dict(machine["environment"])
    values["AMRVAC_DIR"] = machine["amrvac_dir"]
    lines.extend(f"export {key}={shlex.quote(value)}" for key, value in values.items())
    if machine["library_dirs"]:
        joined = shlex.quote(":".join(machine["library_dirs"]))
        lines.append(f'export LD_LIBRARY_PATH={joined}:"${{LD_LIBRARY_PATH:-}}"')
        lines.append(f'export DYLD_LIBRARY_PATH={joined}:"${{DYLD_LIBRARY_PATH:-}}"')
    return lines


def launch_argv(machine, stage):
    values = {"ranks": machine["ranks"][stage], "executable": machine["executable"], "par": "amrvac.par"}
    argv = [token.format(**values) for token in machine["launcher"]]
    joined = " ".join(machine["launcher"])
    if "{executable}" not in joined:
        argv.append(machine["executable"])
    if "{par}" not in joined:
        argv.extend(["-i", "amrvac.par"])
    return argv


def active_stages(case_root):
    case = json.loads((Path(case_root) / "case.resolved.json").read_text(encoding="utf-8"))
    stages = case.get("stages")
    if stages is None:
        stages = [case.get("pfss_stage", "01_pfss"), "02_relax"]
        if case.get("fullmhd_source", "insert") == "insert":
            stages.append("03_insert")
        stages.append("04_fullmhd")
    if not isinstance(stages, list) or stages not in (["01_pfss", "02_relax", "03_insert", "04_fullmhd"], ["01_pfss_global_ct", "02_relax", "03_insert", "04_fullmhd"], ["01_pfss", "02_relax", "04_fullmhd"], ["01_pfss_global_ct", "02_relax", "04_fullmhd"]):
        raise ValueError("case.resolved stages are not a supported workflow chain")
    return stages


def configure_machine(case_root: Path, machine_path: Path) -> dict:
    root = Path(case_root).resolve()
    stages = active_stages(root)
    machine = load_machine(machine_path)
    if (root / "pipeline.state.json").exists() or (root / ".pipeline.lock").exists():
        raise ValueError("machine configuration is frozen after pipeline execution; create a new case for changes")
    for stage in stages:
        directory = root / stage
        if not (directory / "amrvac.par").is_file():
            raise FileNotFoundError(directory / "amrvac.par")
        libraries = machine["stage_libraries"].get(stage, machine["libraries"])
        # GNU Make's lists do not reliably preserve whitespace in library paths.
        # Case and executable paths may contain spaces; library directories may not.
        if any(any(c.isspace() for c in value) for value in machine["library_dirs"]):
            raise ValueError("AMRVAC LIB_DIRS cannot contain whitespace; use a no-space library installation path")
        make = ["# Generated from machine.resolved.json; machine-specific dependencies."]
        make += [f"LIBS += {library}" for library in libraries]
        make += [f"LIB_DIRS += {value}" for value in machine["library_dirs"]]
        (directory / "local.make").write_text("\n".join(make) + "\n", encoding="utf-8")
        lines = ["#!/usr/bin/env bash"]
        if machine["scheduler"]["kind"] == "lsf":
            scheduler = machine["scheduler"]
            lines += [f"#BSUB -n {machine['ranks'][stage]}", f"#BSUB -J rbsl_{stage}"]
            if scheduler.get("queue"):
                lines.append(f"#BSUB -q {scheduler['queue']}")
            if scheduler.get("ptile"):
                lines.append(f'#BSUB -R "span[ptile={scheduler["ptile"]}]"')
        lines += ["set -euo pipefail", 'cd -- "$(dirname -- "${BASH_SOURCE[0]}")"', "mkdir -p log data", "exec > log/runtime.stdout.log 2> log/runtime.stderr.log"]
        # LSF executes a spooled script: BASH_SOURCE then points outside the case.
        # An explicit quoted working directory is therefore required for LSF.
        if machine["scheduler"]["kind"] == "lsf":
            lines[lines.index('cd -- "$(dirname -- "${BASH_SOURCE[0]}")"')] = "cd -- " + shlex.quote(str(directory))
        lines += shell_environment(machine)
        lines.append("exec " + shlex.join(launch_argv(machine, stage)))
        job = directory / "job"
        job.write_text("\n".join(lines) + "\n", encoding="utf-8")
        job.chmod(0o755)
    write_json(root / "machine.resolved.json", machine)
    return {"configured": True, "case_root": str(root), "stages": stages, "scheduler": machine["scheduler"]["kind"], "machine": machine}

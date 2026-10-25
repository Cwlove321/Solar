"""Build, preflight, and sequentially run a case with validated DAT handoffs.

Local execution and blocking LSF submission are supported. No stage is
considered successful from scheduler status alone. Never re-run insertion on an
unvalidated partial stage: create a new case or inspect its failure first.
"""
from __future__ import annotations

import hashlib
import json
import math
import os
import re
import shlex
import signal
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path

from checkpoint_io import CheckpointError, select_checkpoint, sha256_file, validate_checkpoint
from case_inputs import audit_case_inputs
from machine_config import active_stages, launch_argv, shell_environment, write_json
from validate_workflow import parse_assignments

ERROR_PATTERN = re.compile(r"(?i)\b(?:nan|[+-]?inf(?:inity)?)\b|FATAL_NONFINITE_STATE|FATAL_INVALID_TIMESTEP|Time step too small|small value encountered|SIGSEGV|MPI_ABORT|Fortran runtime error")


def _now():
    return datetime.now(timezone.utc).isoformat()


def _load(root):
    machine_path = Path(root) / "machine.resolved.json"
    if not machine_path.is_file():
        raise ValueError("configure this case with a machine JSON before build/run")
    return json.loads(machine_path.read_text(encoding="utf-8"))


def _number(value):
    result = float(value.lower().replace("d", "e"))
    if not math.isfinite(result):
        raise ValueError("nonfinite namelist number")
    return result


def _text(value):
    if len(value) < 2 or value[0] not in "\"'" or value[-1] != value[0]:
        raise ValueError(f"expected a quoted Fortran string, got {value}")
    return value[1:-1].replace(value[0] * 2, value[0])


def stage_info(root, stage):
    directory = Path(root) / stage
    par = parse_assignments(directory / "amrvac.par")
    prefix = _text(par["filelist.base_filename"])
    prefix_path = directory / prefix
    if not prefix_path.resolve().is_relative_to(Path(root).resolve()):
        raise ValueError("output base_filename must stay inside the case")
    # AMRVAC parameters use fractions of a turn; DAT coordinates use radians.
    mesh = {"domain_nx": [int(par[f"meshlist.domain_nx{i}"]) for i in (1, 2, 3)],
            "block_nx": [int(par[f"meshlist.block_nx{i}"]) for i in (1, 2, 3)],
            "refine_max_level": int(par["meshlist.refine_max_level"]),
            "qstretch_baselevel1": _number(par["meshlist.qstretch_baselevel(1)"])}
    for name, stem in (("xmin", "xprobmin"), ("xmax", "xprobmax")):
        values = [_number(par[f"meshlist.{stem}{i}"]) for i in (1, 2, 3)]
        mesh[name] = [values[0], values[1] * 2 * math.pi, values[2] * 2 * math.pi]
    return {"stage": stage, "directory": directory, "par": par, "prefix": prefix_path,
            "target": _number(par["stoplist.time_max"]), "mesh": mesh}


def _fingerprint(root, stages, machine):
    root = Path(root)
    digest = hashlib.sha256()
    files = [root / "case.resolved.json", root / "machine.resolved.json"]
    for stage in stages:
        directory = root / stage
        files += [directory / name for name in ("amrvac.par", "mod_usr.t", "local.make", "job")]
        executable = Path(machine["executable"])
        executable = executable if executable.is_absolute() else directory / executable
        if not executable.is_file():
            raise FileNotFoundError(f"missing compiled stage executable: {executable}")
        digest.update(stage.encode())
        digest.update(sha256_file(executable).encode())
        # Hash actual copied input bytes, not only original source path strings.
        files += [path for name in ("path.dat", "Bphot/Brll0000.fits", "potential_boundary/pfss_boundary.dat") if (path := directory / name).is_file()]
    for path in files:
        digest.update(str(path.relative_to(root)).encode())
        digest.update(sha256_file(path).encode())
    return digest.hexdigest()


def _shell(machine, argv):
    return ["bash", "-lc", "\n".join(["set -euo pipefail", *shell_environment(machine), "exec " + shlex.join(argv)])]


def preflight(case_root: Path, runtime: bool = True) -> dict:
    root = Path(case_root).resolve()
    errors, checks = [], []
    try:
        stages = active_stages(root)
        machine = _load(root) if runtime or (root / "machine.resolved.json").is_file() else None
        case = json.loads((root / "case.resolved.json").read_text(encoding="utf-8"))
    except (ValueError, OSError) as exc:
        return {"ok": False, "runtime_checked": runtime, "errors": [str(exc)], "checks": checks}
    derived_grid = {}
    try:
        derived_grid = audit_case_inputs(root).get("grid", {})
        checks.append("preprocessed FITS, full-sphere boundary, axis, times and derived case grid revalidated")
    except (ValueError, KeyError, OSError) as exc:
        errors.append(f"scientific input audit: {exc}")
    previous_mesh = None
    previous_time = 0.
    shared_science = None
    for stage in stages:
        try:
            info = stage_info(root, stage)
            par = info["par"]
            if previous_mesh and info["mesh"] != previous_mesh:
                errors.append(f"{stage}: domain/block dimensions or bounds differ across restart stages")
            previous_mesh = info["mesh"]
            for key, expected in derived_grid.items():
                if key not in {"refine_max_level", "qstretch_baselevel1"} and not re.fullmatch(r"(?:domain_nx|xprobmin|xprobmax)[123]", key):
                    continue
                namelist_key = "qstretch_baselevel(1)" if key == "qstretch_baselevel1" else key
                actual = _number(par[f"meshlist.{namelist_key}"])
                if not math.isclose(actual, expected, rel_tol=1e-12, abs_tol=1e-12):
                    errors.append(f"{stage}: actual namelist {key}={actual} differs from input-derived grid value {expected}")
            shared = {key: value.lower() for key, value in par.items() if key.startswith("atmosphere_list.") or key in {"mag_list.bflux", "mag_list.rmfr", "mag_list.i_helix", "mag_list.npoints", "mag_list.flip_axis"}}
            if shared_science is not None and shared != shared_science:
                errors.append(f"{stage}: atmosphere or RBSL geometry/flux/helicity settings differ across stages")
            shared_science = shared
            if any(d <= 0 or b <= 0 or d % b for d, b in zip(info["mesh"]["domain_nx"], info["mesh"]["block_nx"])):
                errors.append(f"{stage}: positive domain dimensions must be divisible by positive block dimensions")
            if info["target"] < 0:
                errors.append(f"{stage}: time_max must be nonnegative")
            start_time = 0. if stage == stages[0] or par.get("stoplist.reset_time", ".false.").lower() == ".true." else previous_time
            if info["target"] < start_time - 1e-12:
                errors.append(f"{stage}: time_max={info['target']} precedes inherited restart time={start_time}")
            if stage == "03_insert" and par.get("mag_list.liakh_mass_loading", ".false.").lower() == ".true.":
                start_seconds = _number(par["mag_list.liakh_t_start_s"])
                end_seconds = start_seconds + _number(par["mag_list.liakh_t_load_s"])
                if start_seconds < start_time * 5972.47716155599 - 1e-4 or end_seconds > info["target"] * 5972.47716155599 + 1e-4:
                    errors.append("03_insert: actual code-time interval does not contain the complete Liakh loading window")
            previous_time = info["target"]
            for key in ("mhd_list.mhd_energy", "mhd_list.mhd_internal_e"):
                if par.get(key, "").lower() != ".true.":
                    errors.append(f"{stage}: {key} must be true for the validated internal-energy workflow")
            if _text(par.get("mhd_list.typedivbfix", "''")).lower() != "ct":
                errors.append(f"{stage}: constrained transport is required")
            for filename in (("mod_usr.t", "local.make", "job") if runtime else ("mod_usr.t",)):
                if not (info["directory"] / filename).is_file():
                    errors.append(f"{stage}: missing {filename}")
            if stage != stages[0] and _text(par.get("filelist.restart_from_file", "''")) != "./restart.dat":
                errors.append(f"{stage}: restart_from_file must be './restart.dat'")
            if stage == "03_insert":
                if par.get("filelist.firstprocess", "").lower() != ".true." or int(par.get("meshlist.iprob", "0")) != 2:
                    errors.append("03_insert: firstprocess=true and iprob=2 are required exactly once")
            elif stage in {"02_relax", "04_fullmhd"} and par.get("filelist.firstprocess", "").lower() != ".false.":
                errors.append(f"{stage}: continuation must keep firstprocess=false")
            checks.append(f"{stage}: namelist and handoff settings checked")
        except (ValueError, OSError, KeyError) as exc:
            errors.append(f"{stage}: {exc}")
    # create_case validates FITS/WCS/axis deeply; preserve its hash records and
    # check the copied inputs remain accessible even if original files moved.
    for key in ("input_checksums", "inputs"):
        inputs = case.get(key)
        if isinstance(inputs, dict):
            for name, record in inputs.items():
                if not isinstance(record, dict):
                    continue
                copied = record.get("copied_path") or record.get("case_path")
                if copied:
                    path = root / copied
                    if not path.is_file():
                        errors.append(f"missing copied input {name}: {path}")
                    elif record.get("sha256") and sha256_file(path) != record["sha256"]:
                        errors.append(f"copied input hash changed: {name}")
    first = root / stages[0]
    needed = [first / "Bphot/Brll0000.fits"] if stages[0] == "01_pfss" else [first / "potential_boundary/pfss_boundary.dat"]
    if "03_insert" in stages:
        needed += [first / "path.dat", root / "03_insert/path.dat"]
    errors += [f"missing initial input: {path}" for path in needed if not path.is_file()]
    manifest_path = root / case.get("input_manifest", "input_manifest.json")
    if manifest_path.is_file():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        for relative, expected_hash in manifest.get("files", {}).items():
            path = root / relative
            if not path.resolve().is_relative_to(root) or not path.is_file() or sha256_file(path) != expected_hash:
                errors.append(f"copied input missing or changed: {relative}")
        for stage in stages:
            for stage_relative, input_relative in (("path.dat", "inputs/axis/path.dat"), ("Bphot/Brll0000.fits", "inputs/boundary/Brll0000.fits"), ("potential_boundary/pfss_boundary.dat", "inputs/boundary/pfss_boundary.dat")):
                original = root / input_relative
                copied = root / stage / stage_relative
                if original.is_file() and (not copied.is_file() or sha256_file(original) != sha256_file(copied)):
                    errors.append(f"{stage}: {stage_relative} differs from immutable case input")
    if runtime:
        source = Path(machine["amrvac_dir"])
        requirements = {
            "src/mhd/mod_mhd_phys.t": ("mhd_internal_e", "mhd_trac_type", "'ct'"),
            "src/physics/mod_radiative_cooling.t": ("tcool_cutoff",),
            "src/io/mod_input_output.t": ("snapshot_read_header", "stagger_mark_dat"),
        }
        for relative, tokens in requirements.items():
            path = source / relative
            if not path.is_file():
                errors.append(f"missing AMRVAC source: {path}")
                continue
            text = path.read_text(encoding="utf-8").lower()
            for token in tokens:
                if token not in text:
                    errors.append(f"AMRVAC source lacks required {token}: {path}")
        setup = Path(machine["setup_command"] or source / "setup.pl")
        if not setup.is_file():
            errors.append(f"missing AMRVAC setup tool: {setup}")
        if not (source / "src/vacpp.pl").is_file():
            errors.append(f"missing AMRVAC Perl preprocessor: {source / 'src/vacpp.pl'}")
        if not (source / "arch" / f"{machine['arch']}.defs").is_file():
            errors.append(f"missing AMRVAC architecture definition: {machine['arch']}")
        for directory in machine["library_dirs"]:
            if not Path(directory).is_dir():
                errors.append(f"missing library directory: {directory}")
        libraries = set(machine["libraries"])
        for stage in stages:
            libraries.update(machine["stage_libraries"].get(stage, []))
        for library in sorted(libraries):
            found = any(any(Path(directory).glob(f"lib{library}.{suffix}")) for directory in machine["library_dirs"] for suffix in ("a", "so", "so.*", "dylib"))
            if not found:
                errors.append(f"library lib{library} not found in machine.library_dirs")
        commands = [machine["make_command"], "perl", "rsync", *machine["compiler_commands"], launch_argv(machine, stages[0])[0]]
        if machine["scheduler"]["kind"] == "lsf":
            commands += [machine["scheduler"].get("submit_command", "bsub"), machine["scheduler"].get("cancel_command", "bkill")]
        commands = list(dict.fromkeys(commands))
        try:
            # Load the compiler environment once, rather than once per command.
            lines = ["set -euo pipefail", *shell_environment(machine)]
            lines += [f"if ! command -v -- {shlex.quote(command)} >/dev/null 2>&1; then printf '%s\\n' {shlex.quote('RBSL_MISSING_COMMAND:' + command)}; fi" for command in commands]
            result = subprocess.run(["bash", "-lc", "\n".join(lines)], capture_output=True, text=True, timeout=30)
            if result.returncode:
                errors.append("machine module/environment setup failed: " + result.stderr.strip())
            errors += ["missing command after module setup: " + line.split(":", 1)[1] for line in result.stdout.splitlines() if line.startswith("RBSL_MISSING_COMMAND:")]
        except (OSError, subprocess.TimeoutExpired) as exc:
            errors.append(f"machine tool check: {exc}")
        checks.append("AMRVAC feature tokens, compiler/launcher commands and link library presence checked; actual build still required")
    return {"ok": not errors, "runtime_checked": runtime, "machine_configured": machine is not None, "errors": errors, "checks": checks, "stages": stages}


def build_case(case_root: Path, dry_run: bool = False) -> dict:
    root = Path(case_root).resolve()
    machine = _load(root)
    stages = active_stages(root)
    setup = machine["setup_command"] or str(Path(machine["amrvac_dir"]) / "setup.pl")
    setup_argv = [setup, "-d=3"] if machine["setup_command"] else ["perl", setup, "-d=3"]
    vacpp = "VACPP=perl " + str(Path(machine["amrvac_dir"]) / "src/vacpp.pl")
    plan = [{"stage": stage, "cwd": str(root / stage), "commands": [setup_argv, [machine["make_command"], f"ARCH={machine['arch']}", vacpp, "-j", str(machine["make_jobs"])]]} for stage in stages]
    if dry_run:
        return {"ok": True, "dry_run": True, "plan": plan}
    report = preflight(root, runtime=True)
    if not report["ok"]:
        raise ValueError("preflight failed: " + "; ".join(report["errors"]))
    if (root / "pipeline.state.json").exists() or (root / ".pipeline.lock").exists():
        raise ValueError("build is frozen after pipeline execution; create a new case to rebuild")
    for item in plan:
        log = root / item["stage"] / "log"
        log.mkdir(exist_ok=True)
        with (log / "build.stdout.log").open("w") as stdout, (log / "build.stderr.log").open("w") as stderr:
            for argv in item["commands"]:
                result = subprocess.run(_shell(machine, argv), cwd=item["cwd"], stdout=stdout, stderr=stderr)
                if result.returncode:
                    raise RuntimeError(f"{item['stage']} build failed (exit {result.returncode}); see {log}")
        executable = Path(machine["executable"])
        executable = executable if executable.is_absolute() else root / item["stage"] / executable
        if not executable.is_file() or not os.access(executable, os.X_OK):
            raise RuntimeError(f"build produced no executable: {executable}")
    return {"ok": True, "stages": stages, "plan": plan}


def _errors_in(paths):
    found = []
    for path in paths:
        if not path.is_file():
            continue
        with path.open(encoding="utf-8", errors="replace") as stream:
            for number, line in enumerate(stream, 1):
                if ERROR_PATTERN.search(line):
                    found.append({"file": str(path), "line": number, "text": line.strip()[:2000]})
                    if len(found) >= 20:
                        return found
    return found


def _new_errors(paths, cursors):
    """Read appended bytes only; production iteration logs can be very large."""
    found = []
    for path in paths:
        if not path.is_file():
            continue
        cursor, line_number, carry = cursors.get(str(path), (0, 0, ""))
        if path.stat().st_size < cursor:
            cursor, line_number, carry = 0, 0, ""
        with path.open("rb") as stream:
            stream.seek(cursor)
            raw = stream.read()
            cursor = stream.tell()
        if not raw:
            continue
        lines = (carry + raw.decode("utf-8", errors="replace")).splitlines(keepends=True)
        carry = ""
        for line in lines:
            if ERROR_PATTERN.search(line):
                found.append({"file": str(path), "line": line_number + 1, "text": line.strip()[:2000]})
            if line.endswith(("\n", "\r")):
                line_number += 1
            else:
                carry = line
        cursors[str(path)] = (cursor, line_number, carry)
    return found


def _stop(process, machine, log_directory):
    if machine["scheduler"]["kind"] == "lsf":
        submission = (log_directory / "submission.stdout.log").read_text(errors="replace")
        matches = re.findall(r"Job\s*<(\d+)>", submission)
        if not matches:
            raise RuntimeError("numerical failure detected but LSF job ID was not available; inspect submission log and cancel the job")
        result = subprocess.run(_shell(machine, [machine["scheduler"].get("cancel_command", "bkill"), matches[-1]]), capture_output=True, text=True)
        (log_directory / "cancel.log").write_text(result.stdout + result.stderr)
        if result.returncode:
            raise RuntimeError("LSF cancellation failed; see cancel.log")
    else:
        try:
            os.killpg(process.pid, signal.SIGTERM)
        except ProcessLookupError:
            pass
    try:
        process.wait(timeout=15)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.wait()


def _execute(directory, machine):
    log = directory / "log"
    log.mkdir(exist_ok=True)
    (directory / "data").mkdir(exist_ok=True)
    local = machine["scheduler"]["kind"] == "local"
    stdout_path = log / ("launcher.stdout.log" if local else "submission.stdout.log")
    stderr_path = log / ("launcher.stderr.log" if local else "submission.stderr.log")
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr, (directory / "job").open("r") as job:
        argv = ["bash", "./job"] if local else _shell(machine, [machine["scheduler"].get("submit_command", "bsub"), "-K"])
        process = subprocess.Popen(argv, cwd=directory, stdin=subprocess.DEVNULL if local else job, stdout=stdout, stderr=stderr, start_new_session=True)
        faults = []
        cursors = {}
        try:
            while process.poll() is None:
                faults = _new_errors([log / name for name in ("runtime.stdout.log", "runtime.stderr.log", stdout_path.name, stderr_path.name)] + list((directory / "data").glob("*.log")), cursors)
                if faults:
                    _stop(process, machine, log)
                    break
                time.sleep(0.2)
        except BaseException:
            if process.poll() is None:
                _stop(process, machine, log)
            raise
        faults = faults or _errors_in([log / name for name in ("runtime.stdout.log", "runtime.stderr.log", stdout_path.name, stderr_path.name)] + list((directory / "data").glob("*.log")))
        return process.returncode, faults


def _output_paths(info):
    prefix = info["prefix"]
    return sorted(prefix.parent.glob(prefix.name + "*.dat"))


def run_pipeline(case_root: Path, dry_run: bool = False) -> dict:
    root = Path(case_root).resolve()
    machine = _load(root)
    stages = active_stages(root)
    plan = [{"stage": stage, "ranks": machine["ranks"][stage], "target_time": stage_info(root, stage)["target"], "launch": launch_argv(machine, stage)} for stage in stages]
    if dry_run:
        return {"ok": True, "dry_run": True, "plan": plan, "note": "No job submitted and no restart links or state files created."}
    report = preflight(root, runtime=True)
    if not report["ok"]:
        raise ValueError("preflight failed: " + "; ".join(report["errors"]))
    fingerprint = _fingerprint(root, stages, machine)
    state_path = root / "pipeline.state.json"
    if state_path.exists():
        state = json.loads(state_path.read_text(encoding="utf-8"))
        if state.get("configuration_sha256") != fingerprint:
            raise ValueError("case/machine/source/input configuration changed since execution; create a new case")
    else:
        state = {"schema_version": 1, "configuration_sha256": fingerprint, "started": _now(), "status": "running", "stages": {}}
    lock = root / ".pipeline.lock"
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY, 0o600)
    except FileExistsError as exc:
        raise RuntimeError("another pipeline holds .pipeline.lock; inspect that process before removing a stale lock") from exc
    os.write(fd, f"pid={os.getpid()}\n".encode())
    os.close(fd)
    previous = None
    try:
        for stage in stages:
            info = stage_info(root, stage)
            directory = info["directory"]
            record = state["stages"].get(stage)
            if record:
                if record.get("status") != "validated":
                    raise RuntimeError(f"{stage} has an incomplete/failed prior attempt; create a new case to avoid reapplying insertion/loading")
                checkpoint = Path(record["checkpoint"]["path"])
                if not checkpoint.is_file() or sha256_file(checkpoint) != record["checkpoint"]["sha256"]:
                    raise RuntimeError(f"validated checkpoint changed: {checkpoint}")
                previous = checkpoint
                continue
            if _output_paths(info):
                raise RuntimeError(f"{stage} already contains unregistered DAT output; use a clean case to prevent overwrite/reinsertion")
            executable = Path(machine["executable"])
            executable = executable if executable.is_absolute() else directory / executable
            if not executable.is_file() or not os.access(executable, os.X_OK):
                raise RuntimeError(f"missing compiled executable: {executable}; build the case first")
            if previous is not None:
                restart = directory / "restart.dat"
                relative = os.path.relpath(previous, directory)
                if restart.is_symlink():
                    if restart.resolve() != previous.resolve():
                        raise RuntimeError(f"existing restart link targets another file: {restart}")
                elif restart.exists():
                    raise RuntimeError(f"refusing to overwrite existing restart: {restart}")
                else:
                    restart.symlink_to(relative)
            record = {"status": "running", "started": _now(), "target_time": info["target"], "ranks": machine["ranks"][stage], "executable_sha256": sha256_file(executable)}
            state["stages"][stage] = record
            state["status"] = "running"
            write_json(state_path, state)
            code, faults = _execute(directory, machine)
            record["exit_code"] = code
            record["faults"] = faults
            if faults:
                raise RuntimeError(f"{stage}: numerical/runtime fault detected: {faults[0]['text']}")
            if code:
                raise RuntimeError(f"{stage}: process or scheduler exited {code}")
            checkpoint = select_checkpoint(_output_paths(info), info["target"], info["mesh"])
            record.update(status="validated", finished=_now(), checkpoint=checkpoint)
            write_json(state_path, state)
            previous = Path(checkpoint["path"])
        state.update(status="complete", finished=_now())
        write_json(state_path, state)
        return state
    except BaseException as exc:
        state.update(status="failed", failed=_now(), error=str(exc))
        for record in state["stages"].values():
            if record.get("status") == "running":
                record.update(status="failed", error=str(exc), finished=_now())
        write_json(state_path, state)
        raise
    finally:
        lock.unlink()


def status(case_root: Path) -> dict:
    root = Path(case_root).resolve()
    path = root / "pipeline.state.json"
    if not path.is_file():
        return {"status": "not_started", "case_root": str(root), "machine_configured": (root / "machine.resolved.json").is_file()}
    state = json.loads(path.read_text(encoding="utf-8"))
    state["case_root"] = str(root)
    state["lock_present"] = (root / ".pipeline.lock").exists()
    return state

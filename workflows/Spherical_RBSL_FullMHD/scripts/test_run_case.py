"""Small v5 DAT/fake solver tests. No MPI or real scheduler is executed."""
from __future__ import annotations

import json
import math
import os
import struct
import sys
import tempfile
import textwrap
import unittest
from unittest.mock import patch
from pathlib import Path

import numpy as np

from checkpoint_io import CheckpointError, read_header, select_checkpoint, validate_checkpoint
from machine_config import configure_machine, load_machine
from run_case import build_case, preflight, run_pipeline, status

PYTHON = sys.executable


def make_dat(path, code_time=0., bad=None, staggered=True, geometry="spherical_3D", names=None):
    path = Path(path)
    names = names or ["rho", "m1", "m2", "m3", "e", "b1", "b2", "b3"]
    def name(value):
        return value.encode().ljust(16, b" ")
    header = bytearray(struct.pack("=i9id", 5, 0, 0, len(names), 3, 3, 1, 1, 0, int(code_time * 1000), code_time))
    header += struct.pack("=6d", 1., .1*2*math.pi, .2*2*math.pi, 1.5, .3*2*math.pi, .4*2*math.pi)
    header += struct.pack("=9i", 2, 2, 2, 2, 2, 2, 0, 0, 0)
    header += name(geometry) + struct.pack("=i", int(staggered))
    header += b"".join(name(v) for v in names) + name("mhd")
    header += struct.pack("=i", 1) + struct.pack("=d", 5/3) + name("gamma") + struct.pack("=3i", 1, 0, 0)
    tree_start = len(header)
    blocks_start = tree_start + 28
    struct.pack_into("=2i", header, 4, tree_start, blocks_start)
    tree = struct.pack("=5iq", 1, 1, 1, 1, 1, blocks_start)
    fields = np.zeros((len(names), 8), dtype="=f8")
    fields[names.index("rho")] = 2.
    fields[names.index("e")] = 3.
    if bad == "rho":
        fields[0, 3] = -1.
    elif bad == "energy":
        fields[4, 1] = math.nan
    elif bad == "field":
        fields[7, 4] = math.inf
    ct = np.zeros(3 * 27, dtype="=f8") if staggered else np.zeros(0, dtype="=f8")
    if bad == "ct":
        ct[50] = math.nan
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(header + tree + struct.pack("=6i", 0, 0, 0, 0, 0, 0) + fields.tobytes() + ct.tobytes())
    return path


def make_case(root, *, scheduler="local", failure=None):
    root.mkdir()
    stages = ["01_pfss", "02_relax", "03_insert", "04_fullmhd"]
    (root / "case.resolved.json").write_text(json.dumps({"stages": stages, "pfss_stage": "01_pfss"}))
    source = root.parent / ("fake_source_" + root.name.replace(" ", "_"))
    (source / "src/mhd").mkdir(parents=True)
    (source / "src/physics").mkdir()
    (source / "src/io").mkdir()
    (source / "arch").mkdir()
    (source / "src/mhd/mod_mhd_phys.t").write_text("mhd_internal_e mhd_trac_type 'ct'")
    (source / "src/physics/mod_radiative_cooling.t").write_text("tcool_cutoff")
    (source / "src/io/mod_input_output.t").write_text("snapshot_read_header stagger_mark_dat")
    (source / "arch/default.defs").write_text("F90=mpif90")
    (source / "setup.pl").write_text("#!/usr/bin/env perl\nexit 0;\n")
    (source / "src/vacpp.pl").write_text('open(my $fh, ">", "vacpp-ran.txt") or die $!; print $fh "ran\\n"; close($fh);\n')
    (source / "setup.pl").chmod(0o644)
    (source / "src/vacpp.pl").chmod(0o644)
    for index, stage in enumerate(stages):
        directory = root / stage
        directory.mkdir()
        (directory / "mod_usr.t").write_text("! fake fixture, not a scientific solver\n")
        (directory / "path.dat").write_bytes(b"axis fixture")
        if index == 0:
            (directory / "Bphot").mkdir()
            (directory / "Bphot/Brll0000.fits").write_bytes(b"magnetogram fixture")
        target = index / 10
        par = f"""&filelist
 base_filename='./data/test_'
 restart_from_file='./restart.dat'
 firstprocess={'.true.' if stage=='03_insert' else '.false.'}
/
&stoplist
 time_max={target}
/
&meshlist
 domain_nx1=2
 domain_nx2=2
 domain_nx3=2
 block_nx1=2
 block_nx2=2
 block_nx3=2
 refine_max_level=1
 qstretch_baselevel(1)=1.04
 xprobmin1=1.
 xprobmax1=1.5
 xprobmin2=.1
 xprobmax2=.3
 xprobmin3=.2
 xprobmax3=.4
 iprob={2 if stage=='03_insert' else 1}
/
&mhd_list
 mhd_energy=.true.
 mhd_internal_e=.true.
 typedivbfix='ct'
/
"""
        (directory / "amrvac.par").write_text(par)
        script = f"""#!{PYTHON}
import sys
from pathlib import Path
sys.path.insert(0, {str(Path(__file__).parent)!r})
from test_run_case import make_dat
stage = Path.cwd().name
with Path('ran.txt').open('a') as stream:
 stream.write('ran\\n')
if stage != '01_pfss':
 assert Path('restart.dat').is_file()
target = {target!r}
failure = {failure!r}
if stage == '02_relax' and failure == 'nanlog':
 print('FATAL_NONFINITE_STATE e NaN', flush=True)
 import time; time.sleep(20)
if stage == '02_relax' and failure == 'exit':
 sys.exit(7)
if stage == '02_relax' and failure == 'short':
 target = .01
bad = 'ct' if stage == '02_relax' and failure == 'baddat' else None
make_dat('data/test_0042.dat', target, bad=bad)
"""
        (directory / "amrvac").write_text(textwrap.dedent(script))
        (directory / "amrvac").chmod(0o755)
    machine = {"amrvac_dir": str(source), "arch": "default", "library_dirs": [], "libraries": [], "compiler_commands": ["true"], "launcher": [PYTHON, "{executable}", "-i", "{par}"], "scheduler": {"kind": scheduler}}
    if scheduler == "lsf":
        submit = root / "fake bsub"
        submit.write_text("#!/usr/bin/env bash\nprintf 'Job <123> is submitted\\n'\nbash\n")
        submit.chmod(0o755)
        cancel = root / "fake bkill"
        cancel.write_text("#!/usr/bin/env bash\nexit 0\n")
        cancel.chmod(0o755)
        machine["scheduler"] |= {"submit_command": str(submit), "cancel_command": str(cancel), "queue": "test", "ptile": 2}
    machine_path = root / "machine.json"
    machine_path.write_text(json.dumps(machine))
    configure_machine(root, machine_path)
    return root


class DatTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="checkpoint fixtures ")
        self.root = Path(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def test_all_states_and_ct_validated(self):
        path = make_dat(self.root / "ok.dat", .3)
        report = validate_checkpoint(path, target_time=.3)
        self.assertEqual(report["ranges"]["rho"], [2., 2.])
        self.assertEqual(report["levels"], {"1": 1})
        self.assertEqual(report["ct_range"], [0., 0.])

    def test_bad_cc_or_ct_blocks_handoff(self):
        for bad in ("rho", "energy", "field", "ct"):
            with self.subTest(bad=bad), self.assertRaises(CheckpointError):
                validate_checkpoint(make_dat(self.root / (bad + ".dat"), bad=bad))

    def test_missing_ct_wrong_geometry_or_incomplete_mhd_rejected(self):
        for settings in ({"staggered": False}, {"geometry": "Cartesian_3D"}, {"names": ["rho", "e"]}):
            with self.subTest(settings=settings), self.assertRaises(CheckpointError):
                validate_checkpoint(make_dat(self.root / "incompatible.dat", **settings))

    def test_truncation_and_offsets_rejected(self):
        path = make_dat(self.root / "bad.dat")
        path.write_bytes(path.read_bytes()[:-1])
        with self.assertRaises(CheckpointError):
            validate_checkpoint(path)
        make_dat(path)
        blob = bytearray(path.read_bytes())
        struct.pack_into("=i", blob, 8, len(blob) + 100)
        path.write_bytes(blob)
        with self.assertRaises(CheckpointError):
            read_header(path)

    def test_actual_time_selection_and_target(self):
        old = make_dat(self.root / "output_9999.dat", .1)
        new = make_dat(self.root / "output_0000.dat", .3)
        self.assertEqual(select_checkpoint([old, new], .3)["time"], .3)
        with self.assertRaises(CheckpointError):
            select_checkpoint([old], .3)
        make_dat(new, .3, bad="ct")
        with self.assertRaises(CheckpointError):
            select_checkpoint([old, new], .3)
        make_dat(new, .3)
        make_dat(old, .1, bad="energy")
        with self.assertRaises(CheckpointError):
            select_checkpoint([old, new], .3)


class PipelineTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="runtime_fixtures_")
        self.root = Path(self.tmp.name) / "case with spaces"
        # These are intentionally tiny runtime fixtures, not real FITS/axis
        # inputs. Scientific FITS/axis validation has its own independent tests.
        self.input_audit = patch("run_case.audit_case_inputs", return_value={})
        self.input_audit_mock = self.input_audit.start()

    def tearDown(self):
        self.input_audit.stop()
        self.tmp.cleanup()

    def test_local_complete_relative_links_and_safe_resume(self):
        make_case(self.root)
        self.assertTrue(preflight(self.root)["ok"], preflight(self.root))
        result = run_pipeline(self.root)
        self.assertEqual(result["status"], "complete")
        link = self.root / "03_insert/restart.dat"
        self.assertTrue(link.is_symlink())
        self.assertFalse(os.path.isabs(os.readlink(link)))
        run_pipeline(self.root)
        self.assertEqual((self.root / "03_insert/ran.txt").read_text().count("ran"), 1)

    def test_lsf_blocking_submission_fixture(self):
        make_case(self.root, scheduler="lsf")
        result = run_pipeline(self.root)
        self.assertEqual(result["status"], "complete")
        self.assertIn("Job <123>", (self.root / "01_pfss/log/submission.stdout.log").read_text())

    def test_dry_run_creates_no_state_output_or_restart(self):
        make_case(self.root)
        self.assertTrue(run_pipeline(self.root, dry_run=True)["dry_run"])
        self.assertTrue(build_case(self.root, dry_run=True)["dry_run"])
        self.assertFalse((self.root / "pipeline.state.json").exists())
        self.assertFalse((self.root / "02_relax/restart.dat").exists())
        self.assertFalse((self.root / "01_pfss/data").exists())

    def test_failure_does_not_run_insertion(self):
        for failure in ("exit", "short", "baddat", "nanlog"):
            with self.subTest(failure=failure):
                root = Path(self.tmp.name) / failure
                make_case(root, failure=failure)
                with self.assertRaises((RuntimeError, CheckpointError)):
                    run_pipeline(root)
                self.assertFalse((root / "03_insert/ran.txt").exists())
                self.assertEqual(status(root)["status"], "failed")
                with self.assertRaises(RuntimeError):
                    run_pipeline(root)

    def test_configuration_change_blocks_resume(self):
        make_case(self.root)
        run_pipeline(self.root)
        par = self.root / "03_insert/amrvac.par"
        par.write_text(par.read_text() + "! changed configuration\n")
        with self.assertRaises(ValueError):
            run_pipeline(self.root)

    def test_binary_change_blocks_resume(self):
        make_case(self.root)
        run_pipeline(self.root)
        executable = self.root / "04_fullmhd/amrvac"
        executable.write_text(executable.read_text() + "\n# another build\n")
        with self.assertRaises(ValueError):
            run_pipeline(self.root)

    def test_mutated_checkpoint_blocks_resume(self):
        make_case(self.root)
        run_pipeline(self.root)
        make_dat(self.root / "02_relax/data/test_0042.dat", .1, bad="ct")
        with self.assertRaises(RuntimeError):
            run_pipeline(self.root)

    def test_runtime_false_does_not_need_compiler(self):
        make_case(self.root)
        machine = self.root / "machine.resolved.json"
        value = json.loads(machine.read_text())
        value["compiler_commands"] = ["deliberately_missing_compiler_12345"]
        machine.write_text(json.dumps(value))
        self.assertTrue(preflight(self.root, runtime=False)["ok"])
        self.assertFalse(preflight(self.root, runtime=True)["ok"])

    def test_config_only_before_machine_configuration(self):
        make_case(self.root)
        (self.root / "machine.resolved.json").unlink()
        report = preflight(self.root, runtime=False)
        self.assertTrue(report["ok"], report)
        self.assertFalse(report["machine_configured"])
        self.assertFalse(preflight(self.root, runtime=True)["ok"])

    def test_build_commands_are_machine_configured(self):
        make_case(self.root)
        make = self.root / "fake make"
        make.write_text(f"#!{PYTHON}\nimport sys, shlex, subprocess\nfrom pathlib import Path\nPath('built.txt').write_text('\\n'.join(sys.argv[1:]))\nfor arg in sys.argv[1:]:\n if arg.startswith('VACPP='):\n  subprocess.run(shlex.split(arg[6:]), check=True)\n")
        make.chmod(0o755)
        machine_path = self.root / "machine.json"
        machine = json.loads(machine_path.read_text())
        machine["make_command"] = str(make)
        machine["make_jobs"] = 3
        machine_path.write_text(json.dumps(machine))
        configure_machine(self.root, machine_path)
        source = Path(machine["amrvac_dir"])
        self.assertEqual((source / "setup.pl").stat().st_mode & 0o111, 0)
        self.assertEqual((source / "src/vacpp.pl").stat().st_mode & 0o111, 0)
        self.assertTrue(build_case(self.root)["ok"])
        built = (self.root / "04_fullmhd/built.txt").read_text()
        self.assertIn("ARCH=default", built)
        self.assertIn("3", built)
        self.assertIn("VACPP=perl", built)
        self.assertEqual((self.root / "04_fullmhd/vacpp-ran.txt").read_text(), "ran\n")

    def test_invalid_machine_scheduler_rejected(self):
        make_case(self.root)
        path = self.root / "machine.json"
        machine = json.loads(path.read_text())
        machine["scheduler"]["kind"] = "slurm"
        path.write_text(json.dumps(machine))
        with self.assertRaises(ValueError):
            load_machine(path)

    def test_actual_refinement_and_stretch_must_match_inputs(self):
        self.input_audit_mock.return_value = {"grid": {"refine_max_level": 1, "qstretch_baselevel1": 1.04}}
        for old, new in (("refine_max_level=1", "refine_max_level=2"), ("qstretch_baselevel(1)=1.04", "qstretch_baselevel(1)=1.05")):
            root = Path(self.tmp.name) / ("changed_level" if "refine" in old else "changed_stretch")
            make_case(root)
            for stage in ("01_pfss", "02_relax", "03_insert", "04_fullmhd"):
                par = root / stage / "amrvac.par"
                par.write_text(par.read_text().replace(old, new))
            report = preflight(root, runtime=False)
            self.assertFalse(report["ok"], report)
            self.assertTrue(any("input-derived grid" in error for error in report["errors"]), report)


if __name__ == "__main__":
    unittest.main()

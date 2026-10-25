#!/usr/bin/env python3
"""Configuration tests using disposable copies; never modify the workflow."""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path
import re
import shutil
import tempfile
import unittest
from unittest.mock import patch

import configure_case
from atmosphere_config import DEFAULT_ATMOSPHERE, atmosphere_from_namelist, normalize_atmosphere
from validate_workflow import normalized_subroutine, parse_assignments


SOURCE_ROOT = Path(__file__).resolve().parents[1]
PAR_PATHS = [Path(stage) / "amrvac.par" for stage in (*configure_case.STAGES, configure_case.GLOBAL_STAGE)]
PAR_PATHS += [path.relative_to(SOURCE_ROOT) for path in sorted((SOURCE_ROOT / "04_fullmhd/profiles").glob("*.par"))]


def without_atmosphere(text: str) -> str:
    return re.sub(r"(?ims)^\s*&atmosphere_list\s*\n.*?^\s*/\s*\n", "", text)


class AtmosphereConfigurationTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp = tempfile.TemporaryDirectory(prefix="rbsl_atmosphere_config_test_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        for relative in PAR_PATHS:
            target = self.root / relative
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(SOURCE_ROOT / relative, target)
        self.root_patch = patch.object(configure_case, "ROOT", self.root)
        self.root_patch.start()
        self.addCleanup(self.root_patch.stop)

    def snapshot(self) -> dict[str, bytes]:
        return {str(path.relative_to(self.root)): path.read_bytes() for path in self.root.rglob("*") if path.is_file()}

    def configure(self, payload: dict, atmosphere_only: bool = True) -> None:
        config = self.root / "test_case.json"
        config.write_text(json.dumps(payload), encoding="utf-8")
        with contextlib.redirect_stdout(io.StringIO()):
            configure_case.configure(config, atmosphere_only=atmosphere_only)

    def assert_all_atmospheres(self, expected: dict) -> None:
        for relative in PAR_PATHS:
            self.assertEqual(atmosphere_from_namelist(parse_assignments(self.root / relative)), expected, str(relative))

    def test_defaults_and_roundtrip_all_seven_namelists(self) -> None:
        self.assertEqual(len(PAR_PATHS), 7)
        self.assert_all_atmospheres(normalize_atmosphere())
        self.configure({})
        self.assert_all_atmospheres(DEFAULT_ATMOSPHERE)
        self.assertFalse((self.root / "case.resolved.json").exists())
        self.assertFalse((self.root / "01_pfss/Bphot").exists())

    def test_atmosphere_only_preserves_every_other_parameter(self) -> None:
        before = {relative: without_atmosphere((self.root / relative).read_text()) for relative in PAR_PATHS}
        overrides = {"model": "isothermal", "isothermal_k": 1.2e6, "reference_nh_cm3": 8e8}
        self.configure({"atmosphere": overrides})
        self.assert_all_atmospheres(normalize_atmosphere(overrides))
        for relative, original in before.items():
            self.assertEqual(without_atmosphere((self.root / relative).read_text()), original, str(relative))

    def test_pressure_requires_explicit_anchor_and_accepts_ghost_height(self) -> None:
        for incomplete in ({"normalization": "pressure"}, {"normalization": "pressure", "reference_height_mm": 0}, {"normalization": "pressure", "reference_pressure_cgs": 138.365}):
            with self.assertRaisesRegex(ValueError, "explicit"):
                normalize_atmosphere(incomplete)
        overrides = {"normalization": "pressure", "reference_height_mm": -0.1, "reference_pressure_cgs": 138.365}
        self.configure({"atmosphere": overrides})
        self.assert_all_atmospheres(normalize_atmosphere(overrides))

    def test_zero_conduction_flux_and_legal_model_geometry(self) -> None:
        values = normalize_atmosphere({"conductive_flux_cgs": 0.0, "floor_height_mm": 0.5, "transition_height_mm": 2.0})
        self.assertEqual(values["conductive_flux_cgs"], 0.0)
        with self.assertRaisesRegex(ValueError, "floor_height"):
            normalize_atmosphere({"floor_height_mm": 2.0})

    def test_invalid_configuration_does_not_write_any_file(self) -> None:
        invalid = [None, [], {"unknown": 1}, {"model": "VAL_C"}, {"normalization": "unknown"}, {"normalization": "pressure"}, {"reference_nh_cm3": True}, {"reference_nh_cm3": "1e9"}, {"reference_nh_cm3": float("nan")}, {"isothermal_k": float("inf")}, {"reference_pressure_cgs": 0}, {"transition_width_mm": 0}, {"transition_k": 8000}, {"transition_k": 2e6}, {"conductive_flux_cgs": -1}, {"floor_height_mm": 3}, {"model": "analytic_tr'\n&evil"}]
        for values in invalid:
            with self.subTest(values=values):
                config = self.root / "test_case.json"
                config.write_text(json.dumps({"atmosphere": values}), encoding="utf-8")
                before = self.snapshot()
                with self.assertRaises(ValueError):
                    configure_case.configure(config, atmosphere_only=True)
                self.assertEqual(self.snapshot(), before)

    def test_missing_target_key_preflight_does_not_partially_write(self) -> None:
        target = self.root / configure_case.GLOBAL_STAGE / "amrvac.par"
        target.write_text(target.read_text().replace("atmosphere_isothermal_K", "missing_key"))
        config = self.root / "test_case.json"
        config.write_text(json.dumps({"atmosphere": {"model": "isothermal"}}))
        before = self.snapshot()
        with self.assertRaisesRegex(ValueError, "atmosphere_isothermal_K"):
            configure_case.configure(config, atmosphere_only=True)
        self.assertEqual(self.snapshot(), before)

    def test_full_configuration_does_not_apply_local_grid_to_global_stage(self) -> None:
        payload = json.loads((SOURCE_ROOT / "case.example.json").read_text())
        inputs = self.root / "test_inputs"
        inputs.mkdir()
        (inputs / "Brll0000.fits").write_bytes(b"placeholder copied without FITS parsing")
        (inputs / "path.dat").write_bytes(bytes(payload["rbsl"]["npoints"] * 3 * 8))
        payload["boundary_dir"] = str(inputs)
        payload["path_dat"] = str(inputs / "path.dat")
        payload["grid"]["domain_nx1"] = 80
        payload["atmosphere"] = {"model": "isothermal", "isothermal_k": 1.3e6}
        global_file = self.root / configure_case.GLOBAL_STAGE / "amrvac.par"
        global_original = without_atmosphere(global_file.read_text())
        preserved_physics = {relative: {key: val for key, val in parse_assignments(self.root / relative).items() if key.startswith(("mhd_list.", "rc_list."))} for relative in PAR_PATHS}
        self.configure(payload, atmosphere_only=False)
        self.assert_all_atmospheres(normalize_atmosphere(payload["atmosphere"]))
        self.assertEqual(without_atmosphere(global_file.read_text()), global_original)
        self.assertEqual(parse_assignments(self.root / "01_pfss/amrvac.par")["meshlist.domain_nx1"], "80")
        for relative, before in preserved_physics.items():
            after = parse_assignments(self.root / relative)
            self.assertEqual({key: after[key] for key in before}, before)
        resolved = json.loads((self.root / "case.resolved.json").read_text())
        self.assertEqual(resolved["atmosphere"], normalize_atmosphere(payload["atmosphere"]))

    def test_semantic_namelist_comparison_and_source_layout(self) -> None:
        values = parse_assignments(self.root / "04_fullmhd/amrvac.par")
        values["atmosphere_list.atmosphere_reference_height_mm"] = "1.0D1"
        self.assertEqual(atmosphere_from_namelist(values), DEFAULT_ATMOSPHERE)
        values["atmosphere_list.extra"] = "1"
        with self.assertRaisesRegex(ValueError, "unknown"):
            atmosphere_from_namelist(values)
        a = "SUBROUTINE inithdstatic()\n x = 1 ! comment\nEND SUBROUTINE inithdstatic"
        b = "subroutine inithdstatic()\nx=1\nend subroutine inithdstatic"
        self.assertEqual(normalized_subroutine(a, "inithdstatic"), normalized_subroutine(b, "inithdstatic"))

    def test_static_parser_rejects_hidden_overrides_and_incomplete_groups(self) -> None:
        original = (self.root / "04_fullmhd/amrvac.par").read_text()
        invalid = [
            original + "\n$atmosphere_list\n atmosphere_model='isothermal'\n$end\n",
            original + "\n&atmosphere_list\n atmosphere_model='isothermal'\n/\n",
            original.rstrip().removesuffix("/"),
            original.replace("/\n", "", 1),
            "&atmosphere_list atmosphere_model='isothermal'\n/\n",
            "&atmosphere_list\n atmosphere_model='analytic_tr', atmosphere_normalization='pressure'\n/\n",
            "&meshlist\n domain_nx1=64, &\n 32\n/\n",
            "&meshlist\n domain_nx1=\n 64\n/\n",
            "&meshlist\n domain_nx1=64 /\n",
            "&meshlist\n domain_nx1=64\n continuation without assignment\n/\n",
        ]
        target = self.root / "parser_test.par"
        for text in invalid:
            with self.subTest(text=text[-120:]):
                target.write_text(text)
                with self.assertRaises(ValueError):
                    parse_assignments(target)
        target.write_text("&meshlist\n refine_threshold = 0.2d0, 0.4d0, 16*0.2d0\n/\n&filelist\n name='a!=b/path' ! comment\n/\n")
        values = parse_assignments(target)
        self.assertEqual(values["meshlist.refine_threshold"], "0.2d0, 0.4d0, 16*0.2d0")
        self.assertEqual(values["filelist.name"], "'a!=b/path'")


if __name__ == "__main__":
    unittest.main()

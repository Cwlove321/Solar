"""Exercise the portable case boundary/axis/science contract with actual arrays."""
from __future__ import annotations

import copy
import hashlib
import json
import shutil
import struct
import tempfile
import unittest
from pathlib import Path

import numpy as np
from astropy.io import fits

from create_case import create_case, DYNAMIC_ITFIXGRID, _science
from case_inputs import audit_case_inputs
from verify_dependencies import verify_sources
from validate_workflow import parse_assignments


class CreateCaseTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix="fullmhd_case_fixture_")
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.boundary = self.root/"boundary"
        self.boundary.mkdir()
        self.axisdir = self.root/"axis"
        self.axisdir.mkdir()
        self.axis = self.axisdir/"path.dat"
        self.settings = self.root/"settings.json"
        self.output = self.root/"new_case"
        self.config = {
            "grid":{"domain_nx1":24,"xprobmin1":1.,"xprobmax1":1.5,"refine_max_level":2,"qstretch_baselevel1":1.02},
            "rbsl":{"bflux_mx":2e20,"radius_mm":20.,"positive_helicity":False,"flip_axis":False},
            "filament":{"method":"none"},
            "global_pfss":{"lmax":12},
            "stage_times":{"02_relax":{"time_max_code":.2,"save_interval_code":.01},"03_insert":{"time_max_code":0.},"04_fullmhd":{"time_max_code":.1,"save_interval_code":.005}},
        }
        grid={"nlat":24,"nlon":24,"resolution_degree":5.,"effective_nghost":0,"lon_range":[-60.,60.],"lat_range":[-60.,60.],"x0":24,"y0":6}
        (self.boundary/"parameters.json").write_text(json.dumps({"grid":grid}))
        self.br = np.arange(24*24,dtype="f8").reshape(24,24)/10
        global_br=np.zeros((36,72))
        global_br[6:30,24:48] = self.br.T[::-1,:]
        with (self.boundary/"pfss_boundary.dat").open("wb") as f:
            f.write(struct.pack("<ii",72,36))
            np.radians(180-(np.arange(36)+.5)*5).astype("<f8").tofile(f)
            np.radians((np.arange(72)+.5)*5).astype("<f8").tofile(f)
            global_br.astype("<f8").tofile(f)
        with (self.boundary/"cutout.dat").open("wb") as f:
            f.write(struct.pack("<ii",24,24))
            np.stack([self.br,np.zeros_like(self.br),np.zeros_like(self.br)]).astype("<f8").tofile(f)
        h=fits.Header({"BUNIT":"G","NGHOST":0,"LONMIN":-60.,"LONMAX":60.,"LATMIN":-60.,"LATMAX":60.})
        fits.PrimaryHDU(self.br.astype("f4"),header=h).writeto(self.boundary/"Brll0000.fits")
        np.asarray([[170.,0.,1.],[180.,5.,1.1],[190.,0.,1.]],dtype="<f8").tofile(self.axis)
        (self.axisdir/"path_summary.json").write_text(json.dumps({"npoints":3}))

    def run_case(self):
        self.settings.write_text(json.dumps(self.config))
        return create_case(self.boundary,self.axis,self.settings,self.output)

    def test_generated_case_syncs_every_stage_and_has_auditable_inputs(self):
        self.config.update(pfss_stage="01_pfss_global_ct",thermal_profile="full_thermo",amr_mode="dynamic")
        self.assertEqual(self.run_case(),self.output.resolve())
        resolved=json.loads((self.output/"case.resolved.json").read_text())
        self.assertEqual(resolved["grid"]["domain_nx2"],12)
        self.assertEqual(resolved["rbsl"]["npoints"],3)
        for path in self.output.glob("*/amrvac.par"):
            v=parse_assignments(path)
            self.assertEqual(v["meshlist.domain_nx2"],"12")
            self.assertEqual(v["mag_list.npoints"],"3")
            itfixgrid=int(v["meshlist.itfixgrid"])
            self.assertEqual(itfixgrid,DYNAMIC_ITFIXGRID)
            # AMRVAC's actual freeze predicate is global_time>=tfixgrid OR it>=itfixgrid.
            self.assertFalse(0>=itfixgrid)
            self.assertFalse(600000>=itfixgrid)
            self.assertAlmostEqual(float(v["mag_list.bflux"].replace("d","e")),2e20)
        global_values=parse_assignments(self.output/"01_pfss_global_ct/amrvac.par")
        self.assertEqual(global_values["globalpfss_list.global_patch_x0"],"24")
        self.assertEqual(global_values["globalpfss_list.global_expected_ntheta"],"36")
        self.assertEqual(parse_assignments(self.output/"04_fullmhd/amrvac.par")["mhd_list.mhd_radiative_cooling"],".true.")
        manifest=json.loads((self.output/"input_manifest.json").read_text())["files"]
        for path,expected in manifest.items():
            self.assertEqual(hashlib.sha256((self.output/path).read_bytes()).hexdigest(),expected)
        self.assertEqual((self.output/"03_insert/path.dat").read_bytes(),self.axis.read_bytes())
        par=self.output/"04_fullmhd/amrvac.par"
        digest=hashlib.sha256(par.read_bytes()).hexdigest()
        self.assertIn(digest,par.with_name(par.name+".readme").read_text())
        self.assertEqual(audit_case_inputs(self.output)["npoints"],3)

    def test_no_rope_chain_disables_pfss_subtraction_and_skips_insert(self):
        self.config["fullmhd_source"]="relax"
        self.run_case()
        r=json.loads((self.output/"case.resolved.json").read_text())
        self.assertEqual(r["stages"],["01_pfss","02_relax","04_fullmhd"])
        for stage in ("01_pfss","01_pfss_global_ct","03_insert"):
            self.assertEqual(parse_assignments(self.output/stage/"amrvac.par")["mag_list.do_rbsl"],".false.")

    def test_isolated_case_contains_locked_sources_and_requirement_includes(self):
        self.run_case()
        source_report=verify_sources(self.output)
        self.assertTrue(source_report["ok"],source_report["errors"])
        self.assertTrue((self.output/"vendor/SOURCE_PROVENANCE.md").is_file())
        self.assertTrue((self.output/"tested-environment.txt").is_file())
        self.assertFalse((self.output/"vendor/amrvac/lib").exists())
        self.assertFalse(any((self.output/"vendor").rglob("__pycache__")))
        for requirements in self.output.glob("requirements*.txt"):
            for line in requirements.read_text().splitlines():
                if line.strip().startswith("-r "):
                    self.assertTrue((requirements.parent/line.strip()[3:].strip()).is_file(),f"missing include from {requirements}")

    def test_missing_brll_generated_only_in_new_snapshot(self):
        (self.boundary/"Brll0000.fits").unlink()
        self.run_case()
        self.assertFalse((self.boundary/"Brll0000.fits").exists())
        self.assertTrue((self.output/"inputs/boundary/Brll0000.fits").is_file())

    def test_reject_inconsistent_npoints_and_leave_no_partial_case(self):
        self.config["rbsl"]["npoints"]=6000
        with self.assertRaisesRegex(ValueError,"npoints disagrees"):
            self.run_case()
        self.assertFalse(self.output.exists())

    def test_reject_old_angular_and_preprocessor_radial_overrides(self):
        for key,value,expected in (("xprobmin2",.16,"conflicts"),("xprobmin1",1.001,"solar surface"),("domain_nx2",36,"conflicts"),("domain_nx1",25,"divisible"),("qstretch_baselevel1",1.,"singular"),("xprobmax1",1.7,"radial sampling")):
            with self.subTest(key=key):
                original=dict(self.config["grid"])
                self.config["grid"][key]=value
                with self.assertRaisesRegex(ValueError,expected):
                    self.run_case()
                self.config["grid"]=original
                self.assertFalse(self.output.exists())

    def test_reject_nan_axis_wrong_longitude_and_incorrect_binary_count(self):
        for points in (np.array([[170.,0.,1.],[180.,np.nan,1.1]]),np.array([[-10.,0.,1.],[0.,5.,1.1]])):
            with self.subTest(points=points.tolist()):
                points.astype("<f8").tofile(self.axis)
                (self.axisdir/"path_summary.json").write_text(json.dumps({"npoints":2}))
                with self.assertRaises(ValueError):
                    self.run_case()
        self.axis.write_bytes(b"wrong format")
        with self.assertRaisesRegex(ValueError,"24 bytes"):
            self.run_case()

    def test_reject_fits_ghost_and_global_patch_mismatch(self):
        with fits.open(self.boundary/"Brll0000.fits",mode="update") as f:
            f[0].header["NGHOST"]=2
        with self.assertRaisesRegex(ValueError,"no ghosts"):
            self.run_case()
        with fits.open(self.boundary/"Brll0000.fits",mode="update") as f:
            f[0].header["NGHOST"]=0
            f[0].data[0,0]+=10
        with self.assertRaisesRegex(ValueError,"orientation/values"):
            self.run_case()

    def test_refuses_existing_output_and_invalid_no_rope_load(self):
        self.config["fullmhd_source"]="relax"
        self.config["filament"]["method"]="instant"
        with self.assertRaisesRegex(ValueError,"requires filament.method=none"):
            self.run_case()
        self.output.mkdir()
        marker=self.output/"keep.txt"
        marker.write_text("existing")
        with self.assertRaises(FileExistsError):
            self.run_case()
        self.assertEqual(marker.read_text(),"existing")

    def test_reference_time_mismatch_requires_explicit_override(self):
        metadata=json.loads((self.boundary/"parameters.json").read_text())
        metadata["reference_time"]="2024-10-07T22:00:00"
        (self.boundary/"parameters.json").write_text(json.dumps(metadata))
        (self.axisdir/"path_summary.json").write_text(json.dumps({"npoints":3,"reference_time":"2024-10-07T23:00:00"}))
        with self.assertRaisesRegex(ValueError,"reference_time differ"):
            self.run_case()
        self.config["allow_reference_time_mismatch"]=True
        self.run_case()

    def test_unknown_scientific_settings_are_not_silently_ignored(self):
        self.config["amr_mdoe"]="dynamic"
        with self.assertRaisesRegex(ValueError,"unknown settings"):
            self.run_case()

    def test_restart_targets_follow_actual_reset_clocks(self):
        self.config["stage_times"]["03_insert"].update(reset_time=False,time_max_code=.1)
        with self.assertRaisesRegex(ValueError,"below its restart start time"):
            _science(self.config)
        self.config["fullmhd_source"]="relax"
        # The skipped insertion target is irrelevant to this chain.
        self.config["stage_times"]["04_fullmhd"].update(reset_time=False,time_max_code=.3)
        self.assertEqual(_science(self.config)["stages"],["01_pfss","02_relax","04_fullmhd"])
        self.config["stage_times"]["04_fullmhd"]["time_max_code"]=.1
        with self.assertRaisesRegex(ValueError,"04_fullmhd.*below its restart"):
            _science(self.config)

    def test_liakh_window_must_fit_the_preserved_restart_clock(self):
        self.config["filament"]={"method":"liakh_source","liakh_mass_loading":{"start_time_s":0.,"load_time_s":100.}}
        self.config["stage_times"]["03_insert"].update(reset_time=False,time_max_code=.3)
        with self.assertRaisesRegex(ValueError,"loading window.*actual stage"):
            _science(self.config)
        self.config["filament"]["liakh_mass_loading"]["start_time_s"]=1200.
        _science(self.config)
        self.config["filament"]["liakh_mass_loading"]["load_time_s"]=1000.
        with self.assertRaisesRegex(ValueError,"loading window"):
            _science(self.config)
        self.config["stage_times"]["03_insert"].update(reset_time=True,time_max_code=.018417818440236957)
        self.config["filament"]["liakh_mass_loading"].update(start_time_s=0.,load_time_s=100.)
        _science(self.config)


if __name__=="__main__":
    unittest.main()

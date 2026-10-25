#!/usr/bin/env python3
"""Compile the actual atmosphere routines in a small, MPI-free test harness.

This exercises the shared Fortran physics, not a replacement implementation.
It is not a substitute for compiling and running the complete AMRVAC solver.
Run with: python -m unittest discover -s scripts -p 'test_atmosphere.py' -v
"""

from __future__ import annotations

import re
import shutil
import subprocess
import tempfile
import unittest
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
STAGES = ("01_pfss", "01_pfss_global_ct", "02_relax", "03_insert", "04_fullmhd")
LENGTH = 6.955e10
TEMPERATURE = 1.0e6
NUMBER_DENSITY = 1.0e9
VELOCITY = LENGTH / 5972.47716155599
PRESSURE = 1.99760887990772**2 / (4.0 * np.pi)
DENSITY = PRESSURE / VELOCITY**2
JMAX = 8000
NGHOST = 2


def extract_subroutine(source: str, name: str) -> str:
    match = re.search(
        rf"(?ims)^\s*subroutine\s+{name}\b.*?^\s*end\s+subroutine\s+{name}\b[^\n]*",
        source,
    )
    if match is None:
        raise ValueError(f"missing subroutine {name}")
    return match.group(0)


def driver_source(source: str, *, legacy: bool = False) -> str:
    if legacy:
        declarations = ""
        routines = extract_subroutine(source, "inithdstatic")
        read_call = ""
    else:
        declarations = source.split("! BEGIN ATMOSPHERE PARAMETERS", 1)[1].split(
            "! END ATMOSPHERE PARAMETERS", 1
        )[0]
        routines = "\n".join(
            extract_subroutine(source, name)
            for name in ("read_atmosphere_params", "validate_atmosphere_parameters", "inithdstatic")
        )
        read_call = "call read_atmosphere_params(files)"
    return f"""module atmosphere_harness
  implicit none
  integer, parameter :: dp=kind(1.d0), jmax={JMAX}, nghostcells={NGHOST}
  integer :: mype=0, unitpar=19, refine_max_level=4
  double precision, parameter :: zero=0.d0
  double precision :: unit_length={LENGTH:.16e}_dp
  double precision :: unit_temperature={TEMPERATURE:.16e}_dp
  double precision :: unit_numberdensity={NUMBER_DENSITY:.16e}_dp
  double precision :: unit_pressure={PRESSURE:.16e}_dp
  double precision :: unit_velocity={VELOCITY:.16e}_dp
  double precision :: unit_density={DENSITY:.16e}_dp
  double precision :: ya(jmax),ra(jmax),pa(jmax),dx(3,4)
  double precision :: usr_grav,SRadius,rhob,dr,gzone,rmax
  double precision, allocatable :: rbc(:),pbc(:)
  {declarations}
contains
  subroutine mpistop(message)
    character(len=*), intent(in) :: message
    print *, 'STOP: ', message
    stop 1
  end subroutine mpistop
  {routines}
end module atmosphere_harness
program test_driver
  use atmosphere_harness
  implicit none
  character(len=512), allocatable :: files(:)
  integer :: i, outunit, narg
  narg=command_argument_count()
  allocate(files(narg))
  do i=1,narg
    call get_command_argument(i,files(i))
  end do
  rmax=0.5d0
  gzone=dble(nghostcells)*rmax/(64.d0*2.d0**(refine_max_level-1))
  dr=(2.d0*gzone+rmax)/dble(jmax)
  usr_grav=-2.74d4*unit_length/unit_velocity**2
  SRadius=6.955d10/unit_length
  dx=0.d0
  do i=1,4
    dx(1,i)=rmax/(64.d0*2.d0**(i-1))
  end do
  {read_call}
  call inithdstatic
  open(newunit=outunit,file='profile.bin',access='stream',form='unformatted',status='replace')
  write(outunit) ya,ra,pa,rbc,pbc
  close(outunit)
end program test_driver
"""


def compile_driver(source: str, directory: Path, *, legacy: bool = False) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    driver = directory / "atmosphere_test.f90"
    driver.write_text(driver_source(source, legacy=legacy), encoding="utf-8")
    executable = directory / "atmosphere_test"
    result = subprocess.run(
        ["gfortran", "-std=f2008", "-Wall", "-Wextra", "-fcheck=all", "-ffpe-trap=invalid,zero,overflow",
         "-O0", str(driver), "-o", str(executable)],
        cwd=directory, text=True, capture_output=True,
    )
    if result.returncode:
        raise AssertionError(result.stdout + result.stderr)
    return executable


def run_driver(executable: Path, directory: Path, texts: list[str]):
    directory.mkdir(parents=True, exist_ok=True)
    paths = []
    for index, value in enumerate(texts):
        path = directory / f"input_{index}.par"
        path.write_text(value, encoding="utf-8")
        paths.append(str(path))
    result = subprocess.run([str(executable), *paths], cwd=directory, text=True, capture_output=True)
    data = None
    if result.returncode == 0:
        data = np.fromfile(directory / "profile.bin", dtype=np.float64)
        if data.size != 3 * JMAX + 2 * NGHOST:
            raise AssertionError(f"unexpected profile size {data.size}")
    return result, data


def reference_success_profile(height: np.ndarray):
    """Independent evaluation of the successful RBSL15 analytic atmosphere."""
    tpho, ttop, ttr = 8000 / TEMPERATURE, 1.5e6 / TEMPERATURE, 1.6e5 / TEMPERATURE
    htr, wtr, floor = 2e8 / LENGTH, .25e8 / LENGTH, -.2e8 / LENGTH
    fc = 2e5 / PRESSURE / VELOCITY
    conductivity = 8e-7 * TEMPERATURE**3.5 / LENGTH / DENSITY / VELOCITY**3
    shift = wtr * np.arctanh(2 * (ttr - tpho) / (ttop - tpho) - 1)
    temperature = tpho + .5 * (ttop - tpho) * (np.tanh((np.maximum(height, floor) - htr + shift) / wtr) + 1)
    upper = height > htr
    temperature[upper] = (3.5 * fc / conductivity * (height[upper] - htr) + ttr**3.5)**(2 / 7)
    gravity = (-2.74e4 * LENGTH / VELOCITY**2) / (1 + height)**2
    integral = np.concatenate(([0.0], np.cumsum(np.diff(height) * (gravity[1:] / temperature[1:] + gravity[:-1] / temperature[:-1]) * .5)))
    pressure = 1e5 * temperature[0] * np.exp(integral)
    density = pressure / temperature
    ratio = (6.2e8 / NUMBER_DENSITY) / np.interp(1e9 / LENGTH, height, density)
    return density * ratio, pressure * ratio


@unittest.skipUnless(shutil.which("gfortran"), "gfortran is required for actual Fortran numerical tests")
class AtmosphereFortranTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.temp = tempfile.TemporaryDirectory(prefix="rbsl_atmosphere_test_")
        cls.directory = Path(cls.temp.name)
        cls.executables = {
            stage: compile_driver((ROOT / stage / "mod_usr.t").read_text(), cls.directory / stage)
            for stage in STAGES
        }

    @classmethod
    def tearDownClass(cls):
        cls.temp.cleanup()

    def run_case(self, text: str, stage: str = "04_fullmhd", extra: list[str] | None = None):
        with tempfile.TemporaryDirectory(dir=self.directory) as directory:
            return run_driver(self.executables[stage], Path(directory), [text, *(extra or [])])

    def require_profile(self, text: str, **kwargs):
        result, data = self.run_case(text, **kwargs)
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
        self.assertTrue(np.isfinite(data).all())
        self.assertTrue((data[JMAX:] > 0).all())
        return data

    def test_defaults_match_successful_analytic_model(self):
        data = self.require_profile("&atmosphere_list\n/\n")
        height, density, pressure = data[:3 * JMAX].reshape(3, JMAX)
        expected_density, expected_pressure = reference_success_profile(height)
        np.testing.assert_allclose(density, expected_density, rtol=3e-12)
        np.testing.assert_allclose(pressure, expected_pressure, rtol=3e-12)
        self.assertAlmostEqual(np.interp(1e9 / LENGTH, height, density), .62, places=13)

    def test_all_stages_and_profiles_share_actual_parameters(self):
        baseline = None
        for stage in STAGES:
            data = self.require_profile((ROOT / stage / "amrvac.par").read_text(), stage=stage)
            if baseline is None:
                baseline = data
            else:
                np.testing.assert_array_equal(data, baseline)
        for path in (ROOT / "04_fullmhd" / "profiles").glob("*.par"):
            data = self.require_profile(path.read_text())
            np.testing.assert_array_equal(data, baseline)

    def test_missing_namelist_uses_defaults(self):
        data = self.require_profile("&unrelated\nvalue=1\n/\n")
        expected = self.require_profile("&atmosphere_list\n/\n")
        np.testing.assert_array_equal(data, expected)

    def test_pressure_normalization_uses_position_interpolation(self):
        data = self.require_profile("""&atmosphere_list
 atmosphere_normalization='pressure'
 atmosphere_reference_height_Mm=0.d0
 atmosphere_reference_pressure_cgs=138.365d0
/
""")
        height, density, pressure = data[:3 * JMAX].reshape(3, JMAX)
        self.assertAlmostEqual(np.interp(0, height, pressure) * PRESSURE, 138.365, places=10)
        default = self.require_profile("&atmosphere_list\n/\n")
        np.testing.assert_allclose(pressure / density, default[2 * JMAX:3 * JMAX] / default[JMAX:2 * JMAX], rtol=1e-14)

    def test_isothermal_preserves_hydrostatic_energy_model(self):
        data = self.require_profile("""&atmosphere_list
 atmosphere_model='isothermal'
 atmosphere_isothermal_K=1.2d6
/
""")
        height, density, pressure = data[:3 * JMAX].reshape(3, JMAX)
        np.testing.assert_allclose(pressure / density * TEMPERATURE, 1.2e6, rtol=1e-14)
        self.assertTrue((np.diff(pressure) < 0).all())
        self.assertAlmostEqual(np.interp(1e9 / LENGTH, height, density), .62, places=13)

    def test_case_insensitive_namelist_and_ordered_overrides(self):
        data = self.require_profile("&ATMOSPHERE_LIST\n atmosphere_reference_nH_cgs=1.d9\n/\n",
                                    extra=["&atmosphere_list\n atmosphere_reference_nH_cgs=3.d8\n/\n"])
        height, density, _ = data[:3 * JMAX].reshape(3, JMAX)
        self.assertAlmostEqual(np.interp(1e9 / LENGTH, height, density), .3, places=13)

    def test_invalid_inputs_fail_before_writing_profile(self):
        bad_values = [
            "atmosphere_model='typo'", "atmosphere_normalization='typo'",
            "atmosphere_reference_nH_cgs=-1.d0", "atmosphere_reference_pressure_cgs=0.d0",
            "atmosphere_transition_width_Mm=0.d0", "atmosphere_transition_K=2.d6",
            "atmosphere_reference_height_Mm=1.d6", "atmosphere_conductive_flux_cgs=-1.d0",
            "atmosphere_isothermal_K=-1.d0", "atmosphere_reference_nH_cgs=NaN",
            "unknown_atmosphere_option=1",
        ]
        for assignment in bad_values:
            with self.subTest(assignment=assignment):
                result, data = self.run_case(f"&atmosphere_list\n {assignment}\n/\n")
                self.assertNotEqual(result.returncode, 0)
                self.assertIsNone(data)
        result, _ = self.run_case("&atmosphere_list\n atmosphere_model='isothermal'\n")
        self.assertNotEqual(result.returncode, 0)


if __name__ == "__main__":
    unittest.main()

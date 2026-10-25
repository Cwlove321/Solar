#!/usr/bin/env python3
"""Regress the bottom-only user callback, without compiling MPI-AMRVAC.

The executable test compiles the actual callback statement body with a tiny
Fortran state stub. It tests the user policy, not the full AMR/CT solver.
"""
from __future__ import annotations

from pathlib import Path
import re
import shutil
import subprocess
import tempfile
import unittest

from validate_workflow import normalized_subroutine


ROOT = Path(__file__).resolve().parents[1]
STAGES = ('01_pfss', '01_pfss_global_ct', '02_relax', '03_insert', '04_fullmhd')


def callback(source: str) -> str:
    matches = re.findall(
        r'(?ims)^\s*subroutine special_refine_grid\(.*?^\s*end subroutine special_refine_grid',
        source,
    )
    if len(matches) != 1:
        raise AssertionError('Expected exactly one special_refine_grid')
    return matches[0]


def executable_body(source: str) -> str:
    routine = callback(source)
    # This declaration is the last declaration in all five callback copies.
    parts = re.split(r'(?im)^\s*integer, intent\(inout\)\s*::\s*refine,\s*coarsen\s*$', routine)
    if len(parts) != 2:
        raise AssertionError('Review the callback declaration layout before extracting it')
    return parts[1].rsplit('end subroutine', 1)[0]


def driver_source(body: str) -> str:
    return '''module policy_test
  implicit none
  integer, parameter :: refine_max_level=4
  type :: grid_stub
    logical :: is_physical_boundary(6)=.false.
  end type
  type(grid_stub) :: ps(1)
contains
  subroutine apply_policy(igrid,level,refine,coarsen)
    integer, intent(in) :: igrid,level
    integer, intent(inout) :: refine,coarsen
''' + body + '''
  end subroutine apply_policy
end module policy_test

program check_policy
  use policy_test
  implicit none
  integer :: boundary,level,seed_refine,seed_coarsen,refine,coarsen,checks
  checks=0
  do boundary=0,1
    ps(1)%is_physical_boundary=.false.
    ps(1)%is_physical_boundary(1)=boundary==1
    do level=1,refine_max_level
      do seed_refine=-1,1
        do seed_coarsen=-1,1
          refine=seed_refine
          coarsen=seed_coarsen
          call apply_policy(1,level,refine,coarsen)
          if(boundary==1) then
            if(refine/=1 .or. coarsen/=-1) error stop 1
          else
            if(refine/=seed_refine .or. coarsen/=seed_coarsen) error stop 2
          end if
          checks=checks+1
        end do
      end do
    end do
  end do
  print *, checks
end program check_policy
'''


class RefinementPolicyTests(unittest.TestCase):
    def test_all_five_callbacks_match(self) -> None:
        reference = normalized_subroutine(
            (ROOT / STAGES[0] / 'mod_usr.t').read_text(), 'special_refine_grid'
        )
        for stage in STAGES:
            with self.subTest(stage=stage):
                self.assertEqual(reference, normalized_subroutine(
                    (ROOT / stage / 'mod_usr.t').read_text(), 'special_refine_grid'))

    def test_no_coordinate_box_or_fixed_interior_level(self) -> None:
        expected = 'if(ps(igrid)%is_physical_boundary(1))thenrefine=1coarsen=-1endif'
        for stage in STAGES:
            with self.subTest(stage=stage):
                body = executable_body((ROOT / stage / 'mod_usr.t').read_text())
                uncommented = '\n'.join(line.split('!', 1)[0] for line in body.splitlines())
                self.assertEqual(re.sub(r'\s+', '', uncommented).lower(), expected)

    def compile_and_run(self, body: str, directory: Path) -> subprocess.CompletedProcess:
        compiler = shutil.which('gfortran')
        if not compiler:
            self.skipTest('gfortran unavailable; executable policy test not run')
        source = directory / 'policy.f90'
        source.write_text(driver_source(body), encoding='utf-8')
        executable = directory / 'policy_test'
        build = subprocess.run(
            [compiler, '-std=f2008', '-fcheck=all', '-O0', str(source), '-o', str(executable)],
            cwd=directory, text=True, capture_output=True,
        )
        self.assertEqual(build.returncode, 0, build.stderr)
        return subprocess.run([str(executable)], cwd=directory, text=True, capture_output=True)

    def test_actual_fortran_all_levels_and_request_combinations(self) -> None:
        for stage in STAGES:
            with self.subTest(stage=stage), tempfile.TemporaryDirectory(prefix='rbsl_refine_') as name:
                body = executable_body((ROOT / stage / 'mod_usr.t').read_text())
                result = self.compile_and_run(body, Path(name))
                self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
                self.assertEqual(result.stdout.strip(), '72')

    def test_old_fixed_l3_rule_fails_behavior_check(self) -> None:
        old_body = '''
    if(ps(igrid)%is_physical_boundary(1)) then
      coarsen=-1
      if(level<refine_max_level) refine=1
    else
      if(level<refine_max_level-1) refine=1
      if(level>refine_max_level-1) coarsen=1
    end if
'''
        with tempfile.TemporaryDirectory(prefix='rbsl_refine_old_') as name:
            result = self.compile_and_run(old_body, Path(name))
        self.assertNotEqual(result.returncode, 0)


if __name__ == '__main__':
    unittest.main()

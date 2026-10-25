#!/usr/bin/env python3
"""Create, check, build and run an isolated spherical RBSL case."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys


ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT / 'scripts'))


def emit(value) -> None:
    print(json.dumps(value, ensure_ascii=False, indent=2, default=str))


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest='command', required=True)
    new = commands.add_parser('new', help='Generate an isolated case from boundary and axis products')
    new.add_argument('--boundary', required=True, type=Path, help='Boundary product directory containing parameters.json')
    new.add_argument('--axis', required=True, type=Path, help='Binary path.dat from the path picker')
    new.add_argument('--settings', required=True, type=Path, help='Scientific settings JSON; see case.delivery.example.json')
    new.add_argument('--machine', type=Path, help='Optional machine profile; can also be configured later')
    new.add_argument('--output', required=True, type=Path, help='New case directory; must not already exist')
    machine = commands.add_parser('machine', help='Apply a machine profile to an existing case')
    machine.add_argument('case', type=Path)
    machine.add_argument('--profile', required=True, type=Path)
    check = commands.add_parser('check', help='Check inputs and, by default, the actual build/run environment')
    check.add_argument('case', type=Path)
    check.add_argument('--config-only', action='store_true', help='Check configuration without requiring installed compilers/libraries')
    for name in ('build', 'run'):
        sub = commands.add_parser(name, help=f'{name.capitalize()} the selected case stages')
        sub.add_argument('case', type=Path)
        sub.add_argument('--dry-run', action='store_true', help='Print planned operations without compiling or submitting')
    stat = commands.add_parser('status', help='Show recorded pipeline state and case location')
    stat.add_argument('case', type=Path)
    refresh = commands.add_parser('refresh-docs', help='Regenerate all .par.readme snapshots from actual case values')
    refresh.add_argument('case', type=Path)
    verify = commands.add_parser('verify-package', help='Verify a packaged source manifest')
    verify.add_argument('--root', type=Path, default=ROOT)
    package = commands.add_parser('package', help='Create a clean ZIP source delivery archive')
    package.add_argument('--output', type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        if args.command == 'new':
            from create_case import create_case
            if args.machine:
                from machine_config import load_machine
                load_machine(args.machine.resolve())
            case = create_case(args.boundary.resolve(), args.axis.resolve(), args.settings.resolve(), args.output.resolve())
            if args.machine:
                from machine_config import configure_machine
                configure_machine(case, args.machine.resolve())
            docs = list(case.glob('*/amrvac.par.readme')) + list(case.glob('04_fullmhd/profiles/*.par.readme'))
            emit({'case': str(case), 'parameter_readmes': len(docs),
                  'next': ['python', str(ROOT / 'workflow.py'), 'check', str(case), *([] if args.machine else ['--config-only'])]})
        elif args.command == 'machine':
            from machine_config import configure_machine
            emit(configure_machine(args.case.resolve(), args.profile.resolve()))
        elif args.command == 'check':
            from run_case import preflight
            report = preflight(args.case.resolve(), runtime=not args.config_only)
            emit(report)
            if report.get('errors') or report.get('ok') is False:
                return 2
        elif args.command == 'build':
            from run_case import build_case
            emit(build_case(args.case.resolve(), dry_run=args.dry_run))
        elif args.command == 'run':
            from run_case import run_pipeline
            emit(run_pipeline(args.case.resolve(), dry_run=args.dry_run))
        elif args.command == 'status':
            from run_case import status
            emit(status(args.case.resolve()))
        elif args.command == 'refresh-docs':
            from refresh_par_readmes import refresh_all
            emit({'updated': refresh_all(args.case.resolve())})
        elif args.command == 'verify-package':
            from package_workflow import verify_manifest
            emit(verify_manifest(args.root.resolve()))
        elif args.command == 'package':
            from package_workflow import create_archive
            emit(create_archive(ROOT, args.output.resolve()))
    except ImportError as exc:
        print(f'error: missing Python dependency ({exc}); install with this interpreter: python -m pip install -r requirements.txt', file=sys.stderr)
        return 2
    except (ValueError, OSError, RuntimeError, KeyError) as exc:
        print(f'error: {exc}', file=sys.stderr)
        return 2
    return 0


if __name__ == '__main__':
    raise SystemExit(main())

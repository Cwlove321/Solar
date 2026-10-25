#!/usr/bin/env python3
"""Build and verify a source-only, self-described workflow archive."""
from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import tempfile
import zipfile


MANIFEST = 'DELIVERY_MANIFEST.json'
STAGES = {'01_pfss', '01_pfss_global_ct', '02_relax', '03_insert', '04_fullmhd'}
EXCLUDED_NAMES = {'.git', '__pycache__', '.pytest_cache', '.venv', '.DS_Store', 'runs', 'releases'}


def payload_paths(root: Path) -> list[Path]:
    locked_vendor = None
    source_lock = root / 'vendor/dependencies.lock.json'
    if source_lock.is_file():
        lock = json.loads(source_lock.read_text(encoding='utf-8'))
        locked_vendor = {'vendor/dependencies.lock.json', 'vendor/SOURCE_PROVENANCE.md'}
        for package in lock['packages'].values():
            locked_vendor.update(str(Path(package['path']) / relative) for relative in package['file_sha256'])
    paths = []
    for path in sorted(root.rglob('*')):
        rel = path.relative_to(root)
        if any(part in EXCLUDED_NAMES for part in rel.parts) or path.name == MANIFEST:
            continue
        if rel.parts[0] == 'inputs':
            continue
        if rel.parts[0] == 'vendor' and locked_vendor is not None and rel.as_posix() not in locked_vendor:
            continue
        if path.name.endswith(('.pyc', '.o', '.mod', '.a', '.so', '.dylib')):
            continue
        if path.name in {'case.json', 'case.resolved.json', 'machine.json', 'machine.resolved.json', 'pipeline.state.json'}:
            continue
        if rel.parts[0] in STAGES:
            if any(part in {'data', 'log', 'Bphot'} for part in rel.parts[1:]):
                continue
            if path.name in {'path.dat', 'restart.dat', 'amrvac', 'amrvac.h', 'makefile', 'mod_usr.f', 'output.dat'}:
                continue
            if path.suffix in {'.fits', '.coef', '.dat', '.vtu'}:
                continue
        if path.is_symlink():
            raise ValueError(f'Delivery payload must not contain symbolic links: {rel}')
        if not path.is_file():
            continue
        paths.append(path)
    return paths


def sha256_file(path: Path) -> str:
    result = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(chunk)
    return result.hexdigest()


def build_manifest(root: Path) -> dict:
    if (root / 'vendor/dependencies.lock.json').is_file():
        from verify_dependencies import verify_sources
        source_report = verify_sources(root)
        if not source_report['ok']:
            raise ValueError('Dependency source verification failed: ' + '; '.join(source_report['errors']))
    paths = payload_paths(root)
    if not paths:
        raise ValueError('Empty delivery payload')
    return {
        'schema_version': 1,
        'package': 'Spherical_RBSL_FullMHD',
        'created_utc': datetime.now(timezone.utc).isoformat(),
        'description': 'Source workflow, exact solver source snapshot where bundled, and isolated-case tools. Observational inputs and platform libraries are supplied by the receiver.',
        'file_count': len(paths),
        'files': {path.relative_to(root).as_posix(): {'size': path.stat().st_size, 'sha256': sha256_file(path)} for path in paths},
    }


def verify_manifest(root: Path) -> dict:
    manifest_path = root / MANIFEST
    if not manifest_path.is_file():
        raise ValueError(f'Missing {MANIFEST}; verify an extracted delivery package')
    manifest = json.loads(manifest_path.read_text(encoding='utf-8'))
    if manifest.get('schema_version') != 1 or not isinstance(manifest.get('files'), dict):
        raise ValueError('Unsupported or malformed delivery manifest')
    errors = []
    for name, record in manifest['files'].items():
        rel = Path(name)
        if rel.is_absolute() or '..' in rel.parts:
            raise ValueError(f'Unsafe path in manifest: {name}')
        path = root / rel
        if path.is_symlink() or not path.is_file():
            errors.append(f'missing/nonregular: {name}')
        elif path.stat().st_size != record['size'] or sha256_file(path) != record['sha256']:
            errors.append(f'changed: {name}')
    if errors:
        raise ValueError('Package verification failed: ' + '; '.join(errors[:12]))
    return {'ok': True, 'checked_files': len(manifest['files']), 'package': str(root)}


def create_archive(root: Path, output: Path) -> dict:
    root, output = root.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite delivery archive: {output}')
    if root == output or root in output.parents:
        raise ValueError('Place the archive outside the workflow source tree')
    manifest = build_manifest(root)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix='rbsl_package_', dir=output.parent) as name:
        archive = Path(name) / output.name
        with zipfile.ZipFile(archive, 'w', compression=zipfile.ZIP_DEFLATED, compresslevel=6) as target:
            for relative in manifest['files']:
                target.write(root / relative, 'Spherical_RBSL_FullMHD/' + relative)
            target.writestr('Spherical_RBSL_FullMHD/' + MANIFEST, json.dumps(manifest, ensure_ascii=False, indent=2) + '\n')
        # The destination was checked above; never replace an existing archive.
        with output.open('xb') as destination, archive.open('rb') as source:
            for chunk in iter(lambda: source.read(1024 * 1024), b''):
                destination.write(chunk)
    return {'archive': str(output), 'bytes': output.stat().st_size,
            'sha256': sha256_file(output), 'source_files': manifest['file_count']}

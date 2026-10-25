"""Test distribution integrity and a new-user case creation through the CLI."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
import zipfile

from make_demo_inputs import make_demo
from package_workflow import create_archive, payload_paths, verify_manifest


ROOT = Path(__file__).resolve().parents[1]


class DeliveryTests(unittest.TestCase):
    def test_source_package_roundtrip_and_corruption_detection(self):
        with tempfile.TemporaryDirectory(prefix='rbsl_delivery_package_') as name:
            root = Path(name)
            source = root / 'source'
            (source / '01_pfss/data').mkdir(parents=True)
            (source / '01_pfss/amrvac.par').write_text('&stoplist\n time_max=0.d0\n/\n')
            (source / 'README.md').write_text('source archive fixture\n')
            (source / '01_pfss/data/private.dat').write_bytes(b'must not be distributed')
            (source / '01_pfss/restart.dat').symlink_to('/some/external/checkpoint')
            archive = root / 'delivery.zip'
            result = create_archive(source, archive)
            self.assertEqual(result['source_files'], 2)
            with zipfile.ZipFile(archive) as zipped:
                zipped.extractall(root / 'unpacked')
                self.assertFalse(any('private.dat' in item for item in zipped.namelist()))
            unpacked = root / 'unpacked/Spherical_RBSL_FullMHD'
            self.assertEqual(verify_manifest(unpacked)['checked_files'], 2)
            (unpacked / 'README.md').write_text('corrupt\n')
            with self.assertRaisesRegex(ValueError, 'verification failed'):
                verify_manifest(unpacked)
            with self.assertRaises(FileExistsError):
                create_archive(source, archive)

    def test_package_rejects_unexpected_external_symlink(self):
        with tempfile.TemporaryDirectory(prefix='rbsl_delivery_link_') as name:
            root = Path(name)
            (root / 'external.py').symlink_to('/tmp/not-a-distribution-file')
            with self.assertRaisesRegex(ValueError, 'symbolic links'):
                payload_paths(root)

    def test_fresh_case_cli_and_parameter_docs(self):
        with tempfile.TemporaryDirectory(prefix='rbsl_new_user_') as name:
            root = Path(name)
            inputs = make_demo(root / 'products')
            original = {str(path): hashlib.sha256(path.read_bytes()).hexdigest()
                        for path in (root / 'products').rglob('*') if path.is_file()}
            case = root / 'case with spaces'
            result = subprocess.run([sys.executable, str(ROOT / 'workflow.py'), 'new',
                '--boundary', inputs['boundary'], '--axis', inputs['axis'],
                '--settings', inputs['settings'], '--output', str(case)],
                text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertTrue((case / 'workflow.py').is_file())
            self.assertEqual((case / 'VALIDATION.md').read_bytes(), (ROOT / 'VALIDATION.md').read_bytes())
            resolved = json.loads((case / 'case.resolved.json').read_text())
            self.assertEqual(resolved['grid']['domain_nx2'], 12)
            self.assertEqual(resolved['grid']['domain_nx3'], 12)
            docs = list(case.glob('*/amrvac.par.readme')) + list(case.glob('04_fullmhd/profiles/*.par.readme'))
            self.assertEqual(len(docs), 7)
            for doc in docs:
                par = doc.with_name(doc.name.removesuffix('.readme'))
                self.assertIn(hashlib.sha256(par.read_bytes()).hexdigest(), doc.read_text())
            for path, digest in original.items():
                self.assertEqual(hashlib.sha256(Path(path).read_bytes()).hexdigest(), digest)
            repeated = subprocess.run([sys.executable, str(ROOT / 'workflow.py'), 'new',
                '--boundary', inputs['boundary'], '--axis', inputs['axis'],
                '--settings', inputs['settings'], '--output', str(case)],
                text=True, capture_output=True)
            self.assertNotEqual(repeated.returncode, 0)


if __name__ == '__main__':
    unittest.main()

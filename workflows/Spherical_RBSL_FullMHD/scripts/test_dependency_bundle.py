"""Verify dependency-integrity failure handling and standalone analysis imports."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from verify_dependencies import SOURCE_FEATURES, verify_sources


ROOT = Path(__file__).resolve().parents[1]


class DependencyBundleTests(unittest.TestCase):
    def test_locked_bundle_matches(self):
        report = verify_sources(ROOT)
        self.assertTrue(report["ok"], report["errors"])

    def test_damaged_source_is_rejected(self):
        with tempfile.TemporaryDirectory(prefix="rbsl_source_lock_") as temporary:
            root = Path(temporary)
            solver = root / "vendor/amrvac"
            hashes = {}
            for relative, tokens in SOURCE_FEATURES.items():
                path = solver / relative
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_text("\n".join(tokens))
                hashes[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
            (root / "vendor/dependencies.lock.json").write_text(json.dumps({
                "packages": {"amrvac": {"path": "vendor/amrvac", "file_sha256": hashes}}
            }))
            self.assertTrue(verify_sources(root)["ok"])
            first = solver / next(iter(hashes))
            first.write_text("wrong source")
            result = verify_sources(root)
            self.assertFalse(result["ok"])
            self.assertTrue(result["packages"]["amrvac"]["changed"])

    def test_balance_entry_help_from_unrelated_directory(self):
        with tempfile.TemporaryDirectory(prefix="rbsl_balance_help_") as temporary:
            result = subprocess.run([
                sys.executable, str(ROOT / "scripts/calculate_rbsl_balance_from_stretched_dat.py"), "--help"
            ], cwd=temporary, text=True, capture_output=True)
            self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
            self.assertIn("--grid-config", result.stdout)


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Verify delivered source files without importing optional scientific packages."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys


SOURCE_FEATURES = {
    "src/physics/mod_radiative_cooling.t": ["tcool_cutoff", "rc_split"],
    "src/mhd/mod_mhd_phys.t": ["mhd_internal_e", "mhd_trac_type", "mhd_semirelativistic", "case ('ct')"],
    "src/io/mod_input_output_helper.t": ["version_number = 5", "snapshot_write_header1"],
}


def verify_sources(root: Path) -> dict:
    lock_path = root / "vendor" / "dependencies.lock.json"
    lock = json.loads(lock_path.read_text(encoding="utf-8"))
    report = {"ok": True, "packages": {}, "errors": []}
    for name, package in lock["packages"].items():
        directory = root / package["path"]
        missing, changed = [], []
        for relative, expected in package["file_sha256"].items():
            file = directory / relative
            if not file.is_file() or file.is_symlink():
                missing.append(relative)
            elif hashlib.sha256(file.read_bytes()).hexdigest() != expected:
                changed.append(relative)
        report["packages"][name] = {"files": len(package["file_sha256"]),
                                     "missing": missing, "changed": changed}
        if missing or changed:
            report["errors"].append(f"{name}: {len(missing)} missing, {len(changed)} changed source files")
    solver = root / lock["packages"]["amrvac"]["path"]
    for relative, tokens in SOURCE_FEATURES.items():
        path = solver / relative
        text = path.read_text(encoding="utf-8", errors="replace") if path.is_file() else ""
        for token in tokens:
            if token not in text:
                report["errors"].append(f"amrvac missing required feature {token!r} in {relative}")
    report["ok"] = not report["errors"]
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--json", action="store_true")
    args = parser.parse_args()
    try:
        result = verify_sources(args.root.resolve())
    except (OSError, ValueError, KeyError) as error:
        result = {"ok": False, "packages": {}, "errors": [str(error)]}
    if args.json:
        print(json.dumps(result, ensure_ascii=False, indent=2))
    else:
        for name, item in result["packages"].items():
            print(f"{name}: {item['files']} source files; missing={len(item['missing'])}, changed={len(item['changed'])}")
        for error in result["errors"]:
            print(error, file=sys.stderr)
        print("SOURCE CHECK PASSED" if result["ok"] else "SOURCE CHECK FAILED")
    raise SystemExit(0 if result["ok"] else 1)


if __name__ == "__main__":
    main()

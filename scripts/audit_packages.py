#!/usr/bin/env python3
"""Verify pinned commits and deterministic source-tree hashes in Packages/."""

from __future__ import annotations

import argparse
import fnmatch
import hashlib
import json
import subprocess
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_LOCK = REPOSITORY_ROOT / "Packages" / "dependencies.lock.json"
IGNORED_DIRECTORIES = {
    ".git",
    ".pytest_cache",
    ".ruff_cache",
    ".mypy_cache",
    "__pycache__",
    "build",
    "dist",
}
IGNORED_SUFFIXES = {".a", ".dylib", ".mod", ".o", ".pyc", ".so", ".x"}


def source_tree_files(root: Path, extra_excludes: list[str] | None = None) -> list[Path]:
    """Return deployable source files after excluding generated artifacts."""
    excludes = extra_excludes or []
    files = []
    for path in sorted(root.rglob("*")):
        if not path.is_file():
            continue
        relative = path.relative_to(root)
        relative_text = relative.as_posix()
        if any(part in IGNORED_DIRECTORIES or part.endswith(".egg-info") for part in relative.parts):
            continue
        if path.suffix in IGNORED_SUFFIXES or path.name == ".DS_Store":
            continue
        if any(fnmatch.fnmatch(relative_text, pattern) for pattern in excludes):
            continue
        files.append(path)
    return files


def source_tree_hash(root: Path, extra_excludes: list[str] | None = None) -> tuple[str, int]:
    """Return a stable hash over deployable relative paths and file contents."""

    files = source_tree_files(root, extra_excludes)
    digest = hashlib.sha256()
    for path in files:
        relative_text = path.relative_to(root).as_posix()
        file_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        digest.update(relative_text.encode("utf-8"))
        digest.update(b"\0")
        digest.update(file_hash.encode("ascii"))
        digest.update(b"\n")
    return digest.hexdigest(), len(files)


def git_output(path: Path, *arguments: str) -> str:
    return subprocess.check_output(
        ["git", "-C", str(path), *arguments], text=True, stderr=subprocess.DEVNULL
    ).strip()


def marker_contains_commit(path: Path, marker: str, commit: str) -> bool:
    marker_path = path / marker
    return marker_path.is_file() and commit in marker_path.read_text(errors="replace")


def audit(lock_path: Path, selected: set[str] | None = None, *, show_hashes: bool = False) -> int:
    lock = json.loads(lock_path.read_text())
    failures: list[str] = []
    for name, metadata in lock["packages"].items():
        if selected and name not in selected:
            continue
        package_path = REPOSITORY_ROOT / metadata["path"]
        if not package_path.exists():
            if metadata.get("optional", False):
                print(f"[skip] {name}: optional package is absent")
                continue
            failures.append(f"{name}: missing {package_path}")
            continue

        expected_commit = metadata.get("upstream_commit")
        if expected_commit and metadata.get("management", "").startswith("external-git"):
            git_metadata = package_path / ".git"
            if git_metadata.exists():
                actual_commit = git_output(package_path, "rev-parse", "HEAD")
                if actual_commit != expected_commit:
                    failures.append(
                        f"{name}: commit {actual_commit} != locked {expected_commit}"
                    )
                expected_origin = metadata.get("upstream_url")
                if expected_origin:
                    actual_origin = git_output(package_path, "remote", "get-url", "origin")
                    normalized_actual = actual_origin.removesuffix(".git")
                    normalized_expected = expected_origin.removesuffix(".git")
                    if normalized_actual != normalized_expected:
                        failures.append(
                            f"{name}: origin {actual_origin} != locked {expected_origin}"
                        )
            elif not marker_contains_commit(
                package_path, metadata.get("commit_marker", "UPSTREAM_COMMIT"), expected_commit
            ):
                failures.append(f"{name}: no Git metadata or matching commit marker")

        expected_hash = metadata.get("tree_sha256")
        if expected_hash is not None:
            actual_hash, file_count = source_tree_hash(
                package_path, metadata.get("hash_exclude", [])
            )
            if show_hashes:
                print(f"{name}|files={file_count}|sha256={actual_hash}")
            elif actual_hash != expected_hash:
                failures.append(
                    f"{name}: source hash {actual_hash} != locked {expected_hash}"
                )
            else:
                print(f"[ok] {name}: {file_count} source files")

            if metadata.get("management") == "vendored" and (REPOSITORY_ROOT / ".git").exists():
                repository_paths = [
                    str(path.relative_to(REPOSITORY_ROOT))
                    for path in source_tree_files(package_path, metadata.get("hash_exclude", []))
                ]
                ignored = subprocess.run(
                    ["git", "-C", str(REPOSITORY_ROOT), "check-ignore", "--stdin"],
                    input="\n".join(repository_paths) + "\n",
                    capture_output=True,
                    text=True,
                ).stdout.splitlines()
                if ignored:
                    preview = ", ".join(ignored[:5])
                    failures.append(f"{name}: locked source files are Git-ignored: {preview}")
        elif not show_hashes:
            print(f"[ok] {name}: locked commit available")

        for package_relative, source_relative in metadata.get("overlay_files", {}).items():
            package_file = package_path / package_relative
            source_file = REPOSITORY_ROOT / source_relative
            if not package_file.is_file() or not source_file.is_file():
                failures.append(f"{name}: missing overlay file {package_relative}")
            elif package_file.read_bytes() != source_file.read_bytes():
                failures.append(f"{name}: overlay drift in {package_relative}")

    if failures:
        for failure in failures:
            print(f"[fail] {failure}", file=sys.stderr)
        return 1
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lock", type=Path, default=DEFAULT_LOCK)
    parser.add_argument("--package", action="append", dest="packages")
    parser.add_argument("--show-hashes", action="store_true")
    args = parser.parse_args()
    return audit(args.lock, set(args.packages or []), show_hashes=args.show_hashes)


if __name__ == "__main__":
    raise SystemExit(main())

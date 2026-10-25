#!/usr/bin/env bash
set -euo pipefail

repository_root=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
remote_host=${SOLAR_CODE_REMOTE_HOST:-guoguo}
remote_dir=${SOLAR_CODE_REMOTE_DIR:-/share/home/guoyang1/zhouchang/solar_code}
remote_python=${SOLAR_CODE_REMOTE_PYTHON:-/share/home/guoyang1/apps/miniconda3/envs/zcc/bin/python}
mode=dry-run
run_full_tests=false

usage() {
    echo "Usage: $0 [--dry-run|--apply] [--full-test]"
    echo
    echo "Environment overrides:"
    echo "  SOLAR_CODE_REMOTE_HOST"
    echo "  SOLAR_CODE_REMOTE_DIR"
    echo "  SOLAR_CODE_REMOTE_PYTHON"
}

while (($#)); do
    case "$1" in
        --dry-run) mode=dry-run ;;
        --apply) mode=apply ;;
        --full-test) run_full_tests=true ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
    esac
    shift
done

if [[ ! "$remote_dir" =~ ^/[A-Za-z0-9._/-]+/solar_code$ ]]; then
    echo "Refusing unexpected remote directory: $remote_dir" >&2
    exit 2
fi
if [[ ! "$remote_python" =~ ^/[A-Za-z0-9._/-]+/python[0-9.]*$ ]]; then
    echo "Refusing unexpected remote Python path: $remote_python" >&2
    exit 2
fi

common_excludes=(
    --exclude=.git/
    --exclude=.DS_Store
    --exclude=.pytest_cache/
    --exclude=.ruff_cache/
    --exclude=.mypy_cache/
    --exclude=__pycache__/
    --exclude='*.pyc'
    --exclude='*.egg-info/'
    --exclude=build/
    --exclude='*.o'
    --exclude='*.mod'
    --exclude='*.so'
    --exclude='*.dylib'
    --exclude='*.a'
    --exclude='*.x'
)
core_sources=(
    README.md
    pyproject.toml
    bootstrap_dependencies.sh
    .gitignore
    integrations
    src
    tests
    workflows
    scripts
)
package_sources=(
    Packages/README.md
    Packages/dependencies.lock.json
    Packages/simesh
    Packages/radsyn_tools
    Packages/qsl
    Packages/FastQSL2
)

sync_core() {
    local -a extra_options=("$@")
    rsync -ric "${extra_options[@]}" "${common_excludes[@]}" \
        --exclude=helicity \
        "${core_sources[@]}" \
        "${remote_host}:${remote_dir}/"
}

sync_packages() {
    local -a extra_options=("$@")
    rsync -ric "${extra_options[@]}" "${common_excludes[@]}" \
        --exclude=bin/ \
        "${package_sources[@]}" \
        "${remote_host}:${remote_dir}/Packages/"
}

cd "$repository_root"

if [[ "$mode" == dry-run ]]; then
    echo "Dry run: $repository_root -> ${remote_host}:${remote_dir}"
    sync_core --dry-run --itemize-changes
    sync_packages --dry-run --itemize-changes
    echo "No files were changed. Re-run with --apply to synchronize."
    exit 0
fi

echo "Creating a recoverable remote source backup..."
backup_result=$(
    ssh "$remote_host" bash -s -- "$remote_dir" <<'REMOTE_BACKUP'
set -euo pipefail
sync_remote_dir=$1
cd "$sync_remote_dir"
sync_stamp=$(date +%Y%m%d_%H%M%S)
sync_backup="backups/solar_code_pre_sync_${sync_stamp}.tar.gz"
mkdir -p backups
tar --exclude=.pytest_cache --exclude=__pycache__ --exclude='*.pyc' \
    -czf "$sync_backup" \
    .gitignore README.md bootstrap_dependencies.sh pyproject.toml \
    src tests workflows scripts integrations Packages 2>/dev/null || {
        tar --exclude=.pytest_cache --exclude=__pycache__ --exclude='*.pyc' \
            -czf "$sync_backup" \
            .gitignore README.md bootstrap_dependencies.sh pyproject.toml \
            src tests workflows Packages
    }
tar -tzf "$sync_backup" >/dev/null
sync_bytes=$(stat -c %s "$sync_backup")
sync_hash=$(sha256sum "$sync_backup" | awk '{print $1}')
printf 'backup=%s\nbytes=%s\nsha256=%s\n' "$sync_backup" "$sync_bytes" "$sync_hash"
REMOTE_BACKUP
)
printf '%s\n' "$backup_result"

echo "Synchronizing source files without deleting remote-only files..."
sync_core --itemize-changes
sync_packages --itemize-changes

verification_output=$(mktemp -t solar-code-sync-verify.XXXXXX)
trap 'rm -f "$verification_output"' EXIT
{
    sync_core --dry-run --itemize-changes
    sync_packages --dry-run --itemize-changes
} >"$verification_output"
if [[ -s "$verification_output" ]]; then
    echo "Post-sync verification found remaining differences:" >&2
    sed -n '1,200p' "$verification_output" >&2
    exit 1
fi
echo "Checksum verification passed."

echo "Running remote syntax and core-consistency checks..."
ssh "$remote_host" \
    "cd '$remote_dir' && '$remote_python' -m compileall -q src tests workflows scripts && PYTHONPATH=src '$remote_python' -m pytest -q tests/test_preprocess_core_sync.py"

if [[ "$run_full_tests" == true ]]; then
    echo "Running the complete remote test suite..."
    ssh "$remote_host" \
        "cd '$remote_dir' && PYTHONPATH=src '$remote_python' -m pytest -q"
fi

echo "Synchronization completed successfully."

#!/usr/bin/env bash
set -euo pipefail

repository_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
packages_dir=${PACKAGES_DIR:-"${repository_root}/Packages"}
python_bin=${PYTHON:-/opt/anaconda3/bin/python}
lock_file="${repository_root}/Packages/dependencies.lock.json"
overlay_dir="${repository_root}/integrations/FastQSL2"

if [[ ! -x "$python_bin" ]]; then
    echo "错误：Python 不可执行：$python_bin" >&2
    echo "可通过 PYTHON=/absolute/path/to/python 覆盖。" >&2
    exit 1
fi

for vendored_package in simesh radsyn_tools qsl; do
    if [[ ! -d "${packages_dir}/${vendored_package}" ]]; then
        echo "错误：缺少主仓库应携带的 Packages/${vendored_package}。" >&2
        echo "请先恢复完整的 solar_code checkout，不要从未锁定的上游最新版替代。" >&2
        exit 1
    fi
done

fastqsl_metadata=$(
    "$python_bin" -c 'import json,sys; d=json.load(open(sys.argv[1]))["packages"]["FastQSL2"]; print(d["upstream_url"], d["upstream_commit"], sep="|")' "$lock_file"
)
IFS='|' read -r fastqsl_url fastqsl_commit <<<"$fastqsl_metadata"
fastqsl_dir="${packages_dir}/FastQSL2"

if [[ ! -e "$fastqsl_dir" ]]; then
    echo "[clone] FastQSL2 @ ${fastqsl_commit}"
    git clone "$fastqsl_url" "$fastqsl_dir"
    git -C "$fastqsl_dir" checkout --detach "$fastqsl_commit"
fi

if [[ -d "${fastqsl_dir}/.git" ]]; then
    actual_commit=$(git -C "$fastqsl_dir" rev-parse HEAD)
    if [[ "$actual_commit" != "$fastqsl_commit" ]]; then
        echo "错误：FastQSL2 当前提交 $actual_commit 与锁定提交 $fastqsl_commit 不符。" >&2
        exit 1
    fi
    if git -C "$fastqsl_dir" apply --reverse --check "$overlay_dir/fastqsl-local.patch" >/dev/null 2>&1; then
        echo "[skip] FastQSL2 Python 补丁已应用"
    elif git -C "$fastqsl_dir" apply --check "$overlay_dir/fastqsl-local.patch"; then
        git -C "$fastqsl_dir" apply "$overlay_dir/fastqsl-local.patch"
        echo "[patch] FastQSL2 Python 入口"
    else
        echo "错误：FastQSL2 本地状态无法安全应用固定补丁。" >&2
        exit 1
    fi
fi

mkdir -p "$fastqsl_dir/scripts"
install -m 0644 "$overlay_dir/LOCAL_README.md" "$fastqsl_dir/LOCAL_README.md"
install -m 0644 "$overlay_dir/package.gitignore" "$fastqsl_dir/.gitignore"
install -m 0755 "$overlay_dir/build_fastqsl.sh" "$fastqsl_dir/scripts/build_fastqsl.sh"
printf 'repository=%s\nbranch=main\ncommit=%s\n' "$fastqsl_url" "$fastqsl_commit" >"$fastqsl_dir/UPSTREAM_COMMIT"

"$python_bin" "$repository_root/scripts/audit_packages.py"
echo "依赖源码与锁文件一致：${packages_dir}"

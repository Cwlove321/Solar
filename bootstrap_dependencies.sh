#!/usr/bin/env bash
set -euo pipefail

repository_root=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
packages_dir=${PACKAGES_DIR:-"${repository_root}/Packages"}

clone_if_missing() {
    local url=$1
    local destination=$2
    local name=$3

    if [[ -d "${destination}/.git" || -f "${destination}/pyproject.toml" || -d "${destination}/DDA_method1" ]]; then
        echo "[skip] ${name}: ${destination}"
        return
    fi
    if [[ -e "${destination}" ]]; then
        echo "错误：${destination} 已存在，但不像完整的 ${name} 仓库。" >&2
        exit 1
    fi

    echo "[clone] ${name}: ${url}"
    git clone "${url}" "${destination}"
}

mkdir -p "${packages_dir}"
clone_if_missing \
    https://github.com/Astery0502/simesh.git \
    "${packages_dir}/simesh" \
    simesh
clone_if_missing \
    https://github.com/fuwentai/radsyn_tools.git \
    "${packages_dir}/radsyn_tools" \
    radsyn_tools

echo "外部依赖已就绪：${packages_dir}"

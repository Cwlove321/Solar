#!/usr/bin/env bash
set -euo pipefail

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
repository_root=$(cd "${script_dir}/../.." && pwd)
python_bin=${PYTHON_BIN:-python3}
venv_dir=${VENV_DIR:-"${script_dir}/.venv"}
simesh_dir=${SIMESH_DIR:-"${repository_root}/Packages/simesh"}

if [[ ! -d "${simesh_dir}" ]]; then
    echo "未找到 simesh，从原作者仓库拉取到：${simesh_dir}"
    mkdir -p "$(dirname "${simesh_dir}")"
    git clone https://github.com/Astery0502/simesh.git "${simesh_dir}"
fi

echo "[1/3] 创建 Python 环境：${venv_dir}"
"$python_bin" -m venv "$venv_dir"
venv_python="${venv_dir}/bin/python"

echo "[2/3] 安装 Python 依赖和仓库内新版 simesh"
"$venv_python" -m pip install --upgrade pip setuptools wheel
"$venv_python" -m pip install -r "${script_dir}/requirements.txt"
"$venv_python" -m pip install -e "${simesh_dir}"

echo "[3/3] 验证公开 Python API"
"$venv_python" -c "from simesh.amrvac import open_dataset, read_uniform; print('simesh.amrvac: OK')"

if [[ ${BUILD_FORTRAN:-0} == 1 ]]; then
    echo "检测到 BUILD_FORTRAN=1，开始构建 Fortran 程序。"
    make -C "$script_dir" FC="${FC:-ifort}"
else
    echo "未构建 Fortran。加载 Intel MKL 后可运行："
    echo "  BUILD_FORTRAN=1 FC=ifort bash ${script_dir}/install.sh"
fi

echo "安装完成。Python: ${venv_python}"

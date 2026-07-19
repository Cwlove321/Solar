#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
用法:
  bash run_pipeline.sh INPUT_DIR OUTPUT_DIR [START] [END] [LEVEL]

START/END 是 Python 半开区间 [START, END)，默认处理全部 .dat 文件。
可选环境变量：DAT_PATTERN=*.dat、ALLOW_ANISOTROPIC=0、INTERPOLATION=zero、ERROR_ESTIMATION=0
EOF
}

if (( $# < 2 )); then
    usage
    exit 1
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
input_dir=${1%/}
output_dir=${2%/}
start=${3:-0}
level=${5:-2}
interpolation=${INTERPOLATION:-zero}
error_estimation=${ERROR_ESTIMATION:-0}
pattern=${DAT_PATTERN:-*.dat}
allow_anisotropic=${ALLOW_ANISOTROPIC:-0}

total=$(find "$input_dir" -maxdepth 1 -type f -name "$pattern" | wc -l | tr -d ' ')
if (( total == 0 )); then
    echo "错误：${input_dir} 中没有 .dat 文件。" >&2
    exit 1
fi
end=${4:-$total}

if ! [[ $start =~ ^[0-9]+$ && $end =~ ^[0-9]+$ && $level =~ ^[1-9][0-9]*$ ]]; then
    echo "错误：START、END 必须是非负整数，LEVEL 必须是正整数。" >&2
    exit 1
fi
if (( end <= start || end > total )); then
    echo "错误：范围应满足 0 <= START < END <= ${total}。" >&2
    exit 1
fi
if [[ $error_estimation != 0 && $error_estimation != 1 ]]; then
    echo "错误：ERROR_ESTIMATION 只能是 0 或 1。" >&2
    exit 1
fi

if [[ -n ${PYTHON_BIN:-} ]]; then
    python_bin=$PYTHON_BIN
elif [[ -x "${script_dir}/.venv/bin/python" ]]; then
    python_bin="${script_dir}/.venv/bin/python"
else
    python_bin=python3
fi

if [[ ! -x "${script_dir}/helicity" ]]; then
    echo "错误：找不到 ${script_dir}/helicity。请先加载 MKL 并运行 make。" >&2
    exit 1
fi

echo "[1/2] AMRVAC 预处理：快照 [${start}, ${end})，level=${level}"
anisotropic_option=()
if [[ $allow_anisotropic == 1 ]]; then anisotropic_option=(--allow-anisotropic); fi
"$python_bin" "${script_dir}/process_dat.py" \
    "$input_dir" "$output_dir" "$start" "$end" \
    --level "$level" --interpolation "$interpolation" --pattern "$pattern" "${anisotropic_option[@]}"

echo "[2/2] 磁能与磁螺度计算：快照 [${start}, $((end-1))]"
"${script_dir}/helicity" "$output_dir" "$start" "$((end-1))" "$error_estimation"
echo "完成：${output_dir}/out.binary"

#!/usr/bin/env bash
set -euo pipefail

usage() {
    cat <<'EOF'
用法:
  bash py_jobs.sh GROUPS INPUT_DIR [OUTPUT_DIR] [MODE]

MODE:
  helicity  生成 nlfff/ 和 bd/（默认）
  vtk       只生成 VTK
  all       生成 nlfff、bd、VTK 和 NPZ

可选环境变量:
  DAT_PATTERN=*.dat ALLOW_ANISOTROPIC=0 LEVEL=2 INTERPOLATION=zero PYTHON_BIN=/path/to/python
  BSUB_WALLTIME=24:00 BSUB_CORES=1 BSUB_QUEUE=queue_name
EOF
}

if (( $# < 2 )); then
    usage
    exit 1
fi

groups=$1
input_dir=${2%/}
output_dir=${3:-"${input_dir}/helicity_input"}
mode=${4:-helicity}
level=${LEVEL:-2}
interpolation=${INTERPOLATION:-zero}
pattern=${DAT_PATTERN:-*.dat}
allow_anisotropic=${ALLOW_ANISOTROPIC:-0}
walltime=${BSUB_WALLTIME:-24:00}
cores=${BSUB_CORES:-1}
queue=${BSUB_QUEUE:-}

if ! [[ $groups =~ ^[1-9][0-9]*$ ]]; then
    echo "错误：GROUPS 必须是正整数。" >&2
    exit 1
fi
if [[ $interpolation != zero && $interpolation != linear ]]; then
    echo "错误：INTERPOLATION 只能是 zero 或 linear。" >&2
    exit 1
fi

script_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
script_path="${script_dir}/process_dat.py"
if [[ -n ${PYTHON_BIN:-} ]]; then
    python_bin=$PYTHON_BIN
elif [[ -x "${script_dir}/.venv/bin/python" ]]; then
    python_bin="${script_dir}/.venv/bin/python"
else
    python_bin=python3
fi

case "$mode" in
    helicity) output_options=(--save-nlfff --save-bd) ;;
    vtk)      output_options=(--save-vtk) ;;
    all)      output_options=(--save-nlfff --save-bd --save-vtk --save-npz) ;;
    *)
        echo "错误：MODE 必须是 helicity、vtk 或 all。" >&2
        exit 1
        ;;
esac
if [[ $allow_anisotropic == 1 ]]; then output_options+=(--allow-anisotropic); fi

total=$(find "$input_dir" -maxdepth 1 -type f -name "$pattern" | wc -l | tr -d ' ')
if (( total == 0 )); then
    echo "错误：${input_dir} 中没有 .dat 文件。" >&2
    exit 1
fi
if (( groups > total )); then
    echo "任务组数大于文件数，自动从 ${groups} 调整为 ${total}。"
    groups=$total
fi

log_dir="${output_dir}/logs/preprocess"
job_dir="${log_dir}/jobs"
mkdir -p "$output_dir" "$log_dir" "$job_dir"

size=$((total / groups))
remainder=$((total % groups))
start=0

echo "输入文件: ${total}；LSF 作业: ${groups}；模式: ${mode}；level: ${level}"
for ((group=0; group<groups; group++)); do
    group_size=$size
    if (( group < remainder )); then
        group_size=$((group_size + 1))
    fi
    end=$((start + group_size))

    command=(
        "$python_bin" "$script_path" "$input_dir" "$output_dir" "$start" "$end"
        --level "$level" --interpolation "$interpolation" --pattern "$pattern" "${output_options[@]}"
    )
    printf -v quoted_command '%q ' "${command[@]}"
    job_script="${job_dir}/preprocess_${group}.bsub"

    {
        echo '#!/usr/bin/env bash'
        echo "#BSUB -J helicity_pre_${group}"
        echo "#BSUB -o ${log_dir}/pre_${group}_%J.out"
        echo "#BSUB -e ${log_dir}/pre_${group}_%J.err"
        echo "#BSUB -n ${cores}"
        echo "#BSUB -W ${walltime}"
        if [[ -n $queue ]]; then echo "#BSUB -q ${queue}"; fi
        echo 'set -euo pipefail'
        printf '%s\n' "$quoted_command"
    } > "$job_script"

    bsub < "$job_script"
    start=$end
done

echo "作业脚本保存在 ${job_dir}，日志将写入 ${log_dir}。"

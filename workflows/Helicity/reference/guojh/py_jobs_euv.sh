#!/bin/bash

# ================== 用户配置区域 ==================
# 在此处修改以下变量以控制输出行为（无需命令行传入）
level=4                 # 分辨率倍数 (nx = level * domain_nx)
save_nlfff=false       # 是否保存 NLFFF 二进制文件
save_bd=false         # 是否保存边界数据二进制文件
save_vtk=true          # 是否保存 VTK 文件
save_npz=false          # 是否保存其他变量为 NPZ 文件
# =================================================

# 检查参数：需要分组数 n 和 base_dir
if [ $# -lt 2 ]; then
    echo "Usage: $0 <number_of_groups> <base_directory>"
    echo "  (输出选项请直接编辑脚本开头的变量)"
    exit 1
fi

n=$1
base_dir="$2"

# 确保 base_dir 以斜杠结尾
[[ "${base_dir}" != */ ]] && base_dir="${base_dir}/"

out_dir="${base_dir}out_path/"

script_path="/share/home/guojh/zc_codes/euv_wave/process_dat_new.py"
python_bin="/share/home/guojh/zc_codes/zc/bin/python"
log_dir="/share/home/guojh/zc_codes/euv_wave/logs/py_jobs"

mkdir -p "${out_dir}" "${log_dir}"

# 获取文件列表
files=($(ls -1v "${base_dir}"*.dat 2>/dev/null))
total=${#files[@]}
if [ $total -eq 0 ]; then
    echo "Error: No .dat files found in ${base_dir}"
    exit 1
fi
echo "Found $total .dat files. Splitting into $n groups."

# 计算每组大小
size=$((total / n))
remainder=$((total % n))

# 构建 Python 选项字符串（根据内部变量）
python_opts="--level ${level}"
$save_nlfff && python_opts="${python_opts} --save_nlfff"
$save_bd    && python_opts="${python_opts} --save_bd"
$save_vtk   && python_opts="${python_opts} --save_vtk"
$save_npz   && python_opts="${python_opts} --save_npz"

start=0
for ((i=0; i<n; i++)); do
    if [ $i -lt $remainder ]; then
        group_size=$((size + 1))
    else
        group_size=$size
    fi
    end=$((start + group_size))

    # 完整 Python 命令
    cmd="${python_bin} ${script_path} ${base_dir} ${out_dir} ${start} ${end} ${python_opts}"

    # 创建临时作业脚本
    job_script="pyjob_${i}.bsub"
    cat > ${job_script} <<EOF
#!/bin/bash
#BSUB -J py_${i}
#BSUB -o ${log_dir}/out_%J
#BSUB -e ${log_dir}/err_%J
#BSUB -n 1
#BSUB -W 24:00

source /share/home/guojh/zc_codes/zc/bin/activate

${cmd}
EOF

    bsub < ${job_script}
    rm ${job_script}

    start=$end
done

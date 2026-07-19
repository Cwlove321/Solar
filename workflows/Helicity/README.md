# Cartesian 磁能与磁螺度工作流

本项目把 Cartesian 3D AMRVAC `.dat` 快照重采样到规则网格，生成有限体积磁螺度
程序的输入，并计算总磁能、自由磁能、相对磁螺度 `H_V`、电流承载磁螺度 `H_J`
和底边界无符号磁通量。

项目使用原作者维护的新版 `simesh`：
<https://github.com/Astery0502/simesh>。`Packages/` 不上传到本仓库，
安装脚本会在本地缺失时自动拉取到 `Packages/simesh`。

## 目录结构

```text
solar_code/
├── Packages/
│   └── simesh/                    # 新版 AMRVAC Python/Cython 工具
└── workflows/
    └── Helicity/
        ├── process_dat.py         # 唯一 Python 预处理入口
        ├── run_pipeline.sh        # 单机完整流程
        ├── py_jobs.sh             # LSF 并行预处理
        ├── read_results.py        # out.binary → 屏幕/CSV
        ├── install.sh             # 服务器 Python 环境安装
        ├── Makefile               # Intel Fortran + MKL 构建
        ├── helicity.f90           # 精简后的主程序
        ├── helicity_numerics.f90  # 梯度、DeVore 矢势、Poisson、积分
        ├── helicity_diagnostics.f90
        ├── requirements.txt
        └── reference/guojh/       # 远端已运行版本，只用于追溯
```

正式运行只使用顶层文件。`reference/guojh` 原样保存从服务器下载的测试版本，其中
包含个人绝对路径和 EUV/VTK 默认参数，不应作为新任务入口。

## 最快使用方法

在新的 Linux 服务器上克隆本仓库后：

```bash
cd /path/to/solar_code
bash workflows/Helicity/install.sh
```

安装脚本会在缺失时克隆原作者 `simesh`，然后在 `workflows/Helicity/.venv`
创建独立 Python 环境并以 editable 方式安装。安装完成后先加载服务器的 Intel/MKL 环境，
再编译 Fortran：

```bash
module load intel                    # 模块名以服务器实际配置为准
cd /path/to/solar_code/workflows/Helicity
make FC=ifort
```

然后一条命令处理全部 `.dat` 并计算磁螺度：

```bash
bash run_pipeline.sh /data/amrvac_run /data/amrvac_run/helicity_input
```

输出结果在：

```text
/data/amrvac_run/helicity_input/out.binary
```

转换成 CSV：

```bash
.venv/bin/python read_results.py \
  /data/amrvac_run/helicity_input/out.binary \
  --csv /data/amrvac_run/helicity_input/helicity.csv
```

## 1. Python 安装

要求：

- Python 3.11 或更高版本；
- NumPy 1.23.5–2.x；
- C/C++ 编译器，用来构建 `simesh` 的 Cython 扩展；
- `pyevtk`，仅 VTK 输出使用，但完整安装默认一并安装。

默认安装：

```bash
cd /path/to/solar_code
bash workflows/Helicity/install.sh
```

指定 Python 或虚拟环境位置：

```bash
PYTHON_BIN=/path/to/python3 \
VENV_DIR=/data/venvs/helicity \
bash workflows/Helicity/install.sh
```

也可以手动安装：

```bash
git clone https://github.com/Astery0502/simesh.git Packages/simesh
python3 -m venv workflows/Helicity/.venv
workflows/Helicity/.venv/bin/python -m pip install \
  -r workflows/Helicity/requirements.txt
workflows/Helicity/.venv/bin/python -m pip install -e Packages/simesh
```

把项目复制到另一种操作系统或换 Python 小版本后，应重新安装 `simesh`，不要直接
复用其他平台生成的 `.so` 文件。

## 2. Fortran 构建

程序使用 Intel MKL Poisson/Helmholtz 接口，需要 Intel Fortran 编译器以及完整的
MKL include 和 library。普通 `gfortran + LAPACK` 不能直接构建。

### 推荐：使用 Makefile

加载 Intel 环境后，确认 `MKLROOT` 已设置：

```bash
echo "$MKLROOT"
cd /path/to/solar_code/workflows/Helicity
make FC=ifort
```

使用新版 oneAPI 编译器时可改为：

```bash
make FC=ifx
```

Makefile 会自动完成：

1. 编译 `$MKLROOT/include/mkl_dfti.f90`；
2. 编译 `$MKLROOT/include/mkl_poisson.f90`；
3. 按依赖顺序编译三个 Helicity Fortran 文件；
4. 链接 sequential MKL，生成 `workflows/Helicity/helicity`。

### 在 guoguo 服务器生成 `.mod` 和 `.o`

`zcc` 环境只负责 Python 预处理；Fortran 模块由 Intel 编译器生成，不需要进入
Python/Conda 环境。在 `guoguo` 上已经验证可用的命令是：

```bash
ssh guoguo
module load intel/oneapi_2023

cd /share/home/guoyang1/zhouchang/solar_code/workflows/Helicity
echo "$MKLROOT"
make FC=ifort
```

正常情况下 `MKLROOT` 类似：

```text
/share/apps/intel/oneapi_2023.1.0/mkl/2023.1.0
```

构建完成后检查：

```bash
find build -maxdepth 1 -type f -printf '%f\n' | sort
ls -lh helicity
if ldd helicity | grep -q 'not found'; then
  echo '存在缺失动态库'
else
  echo '动态库解析正常'
fi
```

应该生成：

```text
build/mkl_dft_type.mod
build/mkl_dfti.mod
build/mkl_dfti.o
build/mkl_poisson.mod
build/mkl_poisson.o
build/helicity_numerics.mod
build/helicity_numerics.o
build/helicity_diagnostics.mod
build/helicity_diagnostics.o
build/helicity.o
helicity
```

`.mod` 是编译器读取的模块接口，`.o` 是链接时需要的目标文件。两类文件都与编译器
及版本有关，不应从其他服务器直接复制；更换 `ifort/ifx` 或 oneMKL 版本后应重新
生成。

### 手动逐步编译 `.mod` 文件

一般直接使用 Makefile 即可。如果需要排查某一步，可以严格按照依赖顺序手动执行：

```bash
cd /share/home/guoyang1/zhouchang/solar_code/workflows/Helicity
mkdir -p build

MKL_INCLUDE="$MKLROOT/include"
MKL_LIBDIR="$MKLROOT/lib/intel64"

# 1. 生成 mkl_dft_type.mod、mkl_dfti.mod 和 mkl_dfti.o
ifort -O3 -c "$MKL_INCLUDE/mkl_dfti.f90" \
  -module build -Ibuild -I"$MKL_INCLUDE" \
  -o build/mkl_dfti.o

# 2. 读取 mkl_dfti.mod，生成 mkl_poisson.mod 和 mkl_poisson.o
ifort -O3 -c "$MKL_INCLUDE/mkl_poisson.f90" \
  -module build -Ibuild -I"$MKL_INCLUDE" \
  -o build/mkl_poisson.o

# 3. 读取 mkl_poisson.mod，生成 helicity_numerics.mod
ifort -O3 -c helicity_numerics.f90 \
  -module build -Ibuild -I"$MKL_INCLUDE" \
  -o build/helicity_numerics.o

# 4. 生成 helicity_diagnostics.mod
ifort -O3 -c helicity_diagnostics.f90 \
  -module build -Ibuild -I"$MKL_INCLUDE" \
  -o build/helicity_diagnostics.o

# 5. 编译主程序
ifort -O3 -c helicity.f90 \
  -module build -Ibuild -I"$MKL_INCLUDE" \
  -o build/helicity.o

# 6. 链接可执行文件
ifort -O3 \
  build/mkl_dfti.o \
  build/mkl_poisson.o \
  build/helicity_numerics.o \
  build/helicity_diagnostics.o \
  build/helicity.o \
  -L"$MKL_LIBDIR" \
  -lmkl_intel_lp64 -lmkl_sequential -lmkl_core \
  -lpthread -lm -ldl \
  -o helicity
```

不要改变上述顺序，因为后一个源文件会读取前一步产生的 `.mod`。如需彻底重新编译：

```bash
make clean
make FC=ifort
```

`make clean` 会删除当前目录的 `build/` 和 `helicity` 可执行文件，但不会删除源码、
Python 环境或计算数据。

也可以在安装 Python 时请求构建：

```bash
module load intel
BUILD_FORTRAN=1 FC=ifort bash workflows/Helicity/install.sh
```

如果提示找不到 `mkl_poisson.f90`，说明服务器提供的 MKL 安装不含旧 Poisson
Fortran 接口，需要管理员提供对应 Intel MKL 版本或模块。

## 3. Python 预处理

最简单的调用只需要输入和输出目录：

```bash
workflows/Helicity/.venv/bin/python workflows/Helicity/process_dat.py \
  /data/amrvac_run \
  /data/amrvac_run/helicity_input
```

没有指定任何 `--save-*` 开关时，默认生成 `nlfff/` 和 `bd/`，可以直接交给
`helicity`。默认处理全部按名称排序的 `.dat` 文件，默认 `level=2`。

如果目录还包含边界 `.dat`，可用 `--pattern` 只选择时间快照：

```bash
python process_dat.py INPUT OUTPUT \
  --pattern 'driven[0-9][0-9][0-9][0-9].dat'
```

指定范围和分辨率：

```bash
workflows/Helicity/.venv/bin/python workflows/Helicity/process_dat.py \
  /data/amrvac_run /data/amrvac_run/helicity_input \
  0 80 --level 2 --interpolation zero
```

Python 范围是半开区间 `[start,end)`；上例处理编号 0–79。

### 输出模式

```bash
# Helicity 输入
python process_dat.py INPUT OUTPUT --save-nlfff --save-bd

# VTK
python process_dat.py INPUT OUTPUT --save-vtk

# 全部输出
python process_dat.py INPUT OUTPUT \
  --save-nlfff --save-bd --save-vtk --save-npz
```

只生成 Helicity 输入时，程序只从 AMRVAC 文件读取 `b1,b2,b3`，避免为所有流体
变量分配规则网格，适合高分辨率数据。生成 VTK 时还需要 `rho,m1,m2,m3,e`。

### 分辨率和插值

- `--level 2` 表示每个方向的目标分辨率为 `2 * domain_nx`，三维总单元数约为
  level 1 的 8 倍；
- `--interpolation zero` 是推荐默认值，分片常数、内存低，与已使用流程一致；
- `--interpolation linear` 使用三线性插值和一层 ghost cells，更平滑但更耗内存，
  也会改变尖锐梯度附近的电流；
- Fortran 使用单一网格间距，因此预处理会拒绝 x、y、z 间距不同的网格。

若需要复现旧程序对轻微非等距数据的处理，可显式添加 `--allow-anisotropic`。程序会
打印实际 `(dx,dy,dz)`，并在 Fortran 输入头中统一写入 x 方向的 `dx`；这属于兼容
模式，不会静默启用。

VTK 速度使用 `v_i=m_i/rho*116.45`。模型归一化不同可通过
`--velocity-unit VALUE` 修改；该参数不影响 Helicity 二进制输入。

## 4. LSF 集群并行预处理

`py_jobs.sh` 不包含个人路径，会自动使用 `.venv/bin/python`。命令格式：

```bash
bash py_jobs.sh GROUPS INPUT_DIR [OUTPUT_DIR] [MODE]
```

例如提交 8 个 Helicity 预处理作业：

```bash
cd /path/to/solar_code/workflows/Helicity
bash py_jobs.sh 8 /data/amrvac_run /data/amrvac_run/helicity_input helicity
```

模式：

- `helicity`：生成 `nlfff+bd`，默认；
- `vtk`：只生成 VTK；
- `all`：生成 `nlfff+bd+VTK+NPZ`。

通过环境变量适配服务器：

```bash
LEVEL=4 \
INTERPOLATION=zero \
DAT_PATTERN='driven[0-9][0-9][0-9][0-9].dat' \
PYTHON_BIN=/data/venvs/helicity/bin/python \
BSUB_QUEUE=normal \
BSUB_WALLTIME=24:00 \
bash py_jobs.sh 16 /data/amrvac_run /data/amrvac_run/helicity_input helicity
```

生成的 `.bsub` 文件保存在 `OUTPUT/logs/preprocess/jobs`，便于复现；标准输出和错误
日志位于 `OUTPUT/logs/preprocess`。

## 5. 运行 Fortran

命令格式：

```text
./helicity [INPUT_DIR] [START] [END] [ERROR_ESTIMATION]
```

Fortran 的 `START` 和 `END` 两端都包含。例如：

```bash
./helicity /data/amrvac_run/helicity_input 0 79 0
```

默认值依次是 `./helicity_input`、`0`、`79`、`0`。路径和编号通过命令行传入，
不需要修改源码或重新编译。

`ERROR_ESTIMATION=1` 会额外计算 `divergence_error` 和 `curl_a_error`，并输出：

```text
INPUT_DIR/A/A_XXXX.binary
INPUT_DIR/A/Ap_XXXX.binary
INPUT_DIR/potential/potential_XXXX.binary
```

## 6. 一键完整流程

`run_pipeline.sh` 将 Python 的半开区间自动转换为 Fortran 的闭区间：

```bash
bash run_pipeline.sh INPUT_DIR OUTPUT_DIR [START] [END] [LEVEL]
```

例如：

```bash
INTERPOLATION=zero ERROR_ESTIMATION=0 \
bash run_pipeline.sh /data/amrvac_run /data/helicity_output 0 80 2
```

这会预处理第 0–79 个 `.dat`，然后运行 `helicity ... 0 79 0`。

## 7. 文件格式与结果

`nlfff/nlfff_XXXX.binary`：

1. 一个 `float64`：网格间距 `dx`；
2. 三个 `int32`：`nx,ny,nz`；
3. 三个连续 `float64` 数组：`Bx,By,Bz`。

`bd/bd_XXXX.binary` 保存六个外法向边界面：

```text
-Bx(xmin), +Bx(xmax), -By(ymin), +By(ymax), -Bz(zmin), +Bz(zmax)
```

新版 `simesh` 返回 `(x,y,z,variable)`；预处理会显式转换成 Fortran 文件接口使用的
`(z,y,x)` 内存顺序。

`out.binary` 默认每帧包含 5 个 `float64`：

```text
energy, free_energy, H_V, H_J, b_flux
```

误差估计开启时，每帧再增加：

```text
divergence_error, curl_a_error
```

读取或导出 CSV：

```bash
.venv/bin/python read_results.py OUTPUT/out.binary --start 0

.venv/bin/python read_results.py OUTPUT/out.binary \
  --start 0 --error-estimation --csv OUTPUT/helicity.csv
```

## 8. Fortran 优化说明

远端已经运行过的单文件版本完整保存在 `reference/guojh/helicity.f90`。正式版本没有
改变势场边界条件、DeVore gauge 积分、磁能/磁螺度公式和输出量，主要整理为：

- `helicity.f90`：参数、I/O 和逐帧流程；
- `helicity_numerics.f90`：可复用数值核；
- `helicity_diagnostics.f90`：可选误差估计；
- 旋度计算复用两个临时三维数组，替代原来的六个，减少该步骤约 2/3 scratch 内存；
- 每帧检查输入文件、维度和 MKL 返回状态；
- 输入路径、编号和误差开关改为命令行参数，不再为不同任务修改源码。

## 9. 运行前检查

1. `.dat` 必须是 Cartesian 3D AMRVAC 数据；
2. Helicity 至少需要 `b1,b2,b3`；
3. Python 日志中的分辨率应符合预期；
4. `nlfff_XXXX.binary` 与 `bd_XXXX.binary` 必须一一对应；
5. 注意 Python 的 `END` 不包含，Fortran 的 `END` 包含；
6. Fortran 磁能使用 `B²/(8π)`，磁场和长度必须采用一致的物理单位；
7. 高分辨率数据先用单帧和 `--level 1` 检查内存与输出目录。

# AMRVAC 辐射合成

本项目集中存放从 AMRVAC `.dat` 快照生成合成辐射图像的工作流入口。
通用的物理量转换、探测器平面和视线积分辅助函数仍位于
`src/solarutils/sim/radiation.py`，以便测试和其他项目复用。

Python 入口按模拟坐标系明确命名，避免把笛卡尔数据读取器用于球坐标拉伸网格：

| 坐标网格 | 推荐入口 | 输入处理 |
| --- | --- | --- |
| 笛卡尔 AMRVAC | `synthesize_radiation_cartesian.py` | 批量均匀笛卡尔 slab |
| 笛卡尔 AMRVAC | `synthesize_aia_optically_thick_cartesian.py` | 单帧笛卡尔 `read_uniform` |
| 球坐标拉伸 AMRVAC | `synthesize_aia_spherical_stretched.py` | 用 `amrvac.par` 重建真实 `(r,theta,phi)` 后重采样 |

原来的 `synthesize_radiation.py` 和 `synthesize_aia_optically_thick.py` 保留为笛卡尔
兼容入口；新脚本和文档应优先使用带坐标后缀的名称。

目前有两条互补的计算路径：

| 路径 | 适用场景 | 主要输出 |
| --- | --- | --- |
| AMRVAC 原生 `mod_thermal_emission` | AIA/EUV 等光学薄辐射，直接并行读取 AMR `.dat` | VTI/VTU |
| 本目录 Python + `radsyn_tools` | 任意视角、光学薄或包含吸收的光学厚辐射传输 | NPY/PNG/JSON |

AMRVAC 原生 EUV 模块计算的是光学薄视线积分，不处理冷物质吸收。需要光学厚度、
背景强度 `I0` 或任意吸收系数时，应使用后文的 Python/DDA 工作流。

## AMRVAC 原生辐射合成

### 1. 模块位置与使用条件

AMRVAC 3.3 的原生热辐射模块位于：

```text
src/physics/mod_thermal_emission.t
```

辐射合成通过 AMRVAC 的 `.dat` 转换模式调用，不需要单独编写 Fortran 主程序。
建议使用生成模拟快照时的同一个 `amrvac` 可执行文件，因为读取 `.dat` 时必须保持以下内容兼容：

- 三维设置、变量数目和变量排列；
- MHD/HD 物理模块及能量方程设置；
- `mod_usr.t` 中的单位和用户参数；
- AMRVAC `.dat` 文件格式版本。

若需要重新编译，应在原模拟目录设置 AMRVAC 源码路径后编译。例如：

```bash
export AMRVAC_DIR=/path/to/amrvac
cd /path/to/simulation
perl "$AMRVAC_DIR/setup.pl" -d=3
make -j 8
```

服务器使用 Intel MPI 编译的可执行文件时，运行阶段也必须加载相同的 Intel MPI
运行库。可以先用 `ldd ./amrvac` 检查实际链接的 MPI 库。

### 2. 单帧参数文件

下面是 AIA 193 Å、视线 `[0,1,0]` 的核心配置。实际使用时，应在原模拟
`amrvac.par` 的基础上修改 `filelist` 和添加 `emissionlist`，保留与模拟一致的
`meshlist`、`mhd_list` 和用户参数。

```fortran
&filelist
    convert = .true.
    convert_type = 'EIvtiCCmpi'
    autoconvert = .false.
    firstprocess = .false.
    restart_from_file = '/path/to/eru0001.dat'
/

&emissionlist
    filename_euv = '/path/to/output/euv193_y_eru0001_'
    wavelength = 193
    LOS_theta = 90.0d0
    LOS_phi = -90.0d0
    image_rotate = 0.0d0
    dat_resolution = .false.
    activate_unit_arcsec = .true.
/
```

参数含义：

| 参数 | 含义 |
| --- | --- |
| `convert` | 进入后处理转换模式，不继续时间演化。 |
| `convert_type='EIvtiCCmpi'` | 并行计算 EUV image，并输出 VTK ImageData (`.vti`)。 |
| `restart_from_file` | 要处理的 AMRVAC `.dat` 快照。 |
| `filename_euv` | 输出文件前缀；建议包含波段、视角和输入帧编号，防止覆盖。 |
| `wavelength` | AIA 波段，原生模块支持 94、131、171、193、211、304 和 335 Å。 |
| `LOS_theta/LOS_phi` | 视线的球坐标角，单位为度，定义见下一节。 |
| `image_rotate` | 探测器平面内的旋转角，单位为度。 |
| `dat_resolution` | 是否使用沿模拟坐标轴的原生数据分辨率快捷算法。 |
| `activate_unit_arcsec` | 将 VTI 图像坐标输出为 arcsec。 |
| `instrument_resolution_factor` | 调整仪器空间/光谱采样，默认值为 1。 |

不要把非官方参数 `ray_method` 写入 AMRVAC 3.3 的 `emissionlist`。标准 3.3
不认识该变量，会报 `invalid reference to variable in NAMELIST input`。

### 3. 视线方向的定义

笛卡尔网格中，AMRVAC 根据下面的关系生成实际视线向量：

```text
LOS = (-cos(phi) sin(theta),
       -sin(phi) sin(theta),
       -cos(theta))
```

其中 `theta=LOS_theta`、`phi=LOS_phi`，角度单位为度。常见视线对应关系如下：

| 目标视线 | `LOS_theta` | `LOS_phi` |
| --- | ---: | ---: |
| `[+1,0,0]` | 90 | 180 |
| `[-1,0,0]` | 90 | 0 |
| `[0,+1,0]` | 90 | -90 或 270 |
| `[0,-1,0]` | 90 | 90 |
| `[0,0,+1]` | 180 | 任意 |
| `[0,0,-1]` | 0 | 任意 |

程序启动后应在日志中核对，而不是只相信参数值。例如本项目 `[0,1,0]` 的正确日志为：

```text
LOS vector: [ 0.00  1.00  0.00]
xI1 vector: [ 1.00  0.00  0.00]
xI2 vector: [ 0.00  0.00  1.00]
```

这表示探测器横轴是 `+x`、纵轴是 `+z`。光学薄积分对视线正负方向的总强度通常相同，
但探测器横轴可能发生镜像；比较两套程序时仍应使用完全一致的 LOS 和图像基向量。

### 4. `dat_resolution` 的选择

AMRVAC 3.3 中 `dat_resolution=.true.` 的快捷路径只接受源码中明确列出的三个方向：

- `theta=90, phi=0`，即 `[-1,0,0]`；
- `theta=90, phi=90`，即 `[0,-1,0]`；
- `theta=0`，即 `[0,0,-1]`。

其他方向应设为：

```fortran
dat_resolution = .false.
```

这时模块会在探测器平面上按 AIA 仪器像元执行通用笛卡尔视线积分。比如本项目要求
精确使用 `[0,1,0]`，因此不能把 `phi=-90` 与 `dat_resolution=.true.` 组合，否则
会出现 `Wrong LOS for synthesizing emission`。

如果只关心光学薄总强度，也可以用 `[0,-1,0]` 的快捷路径计算，再在绘图时修正横轴镜像；
但为了与 Python/DDA 结果逐像素比较，推荐使用精确的 `[0,1,0]` 通用路径。

### 5. 运行单帧

参数文件假设命名为 `euv193_y.par`：

```bash
mpirun -np 128 ./amrvac -i euv193_y.par \
  > logs/eru0001.out 2> logs/eru0001.err
```

成功日志中应出现：

```text
convert specified: convert=T
Systhesizing EUV image
Wavelength: 193.500 Angstrom
Unit of EUV flux: DN s^-1 pixel^-1
LOS vector: [ 0.00  1.00  0.00]
```

`EIvtiCCmpi` 的输出是 `.vti`，辐射数组位于 CellData 中，名称类似 `AIA193`。
文件可直接用 ParaView 打开。

### 6. guoguo 集群上的 MPI 任务

当前 guoguo 的 `mpi` 队列必须按完整 64 核节点申请。单帧原生合成默认使用一整个
64 核节点：

```bash
#!/bin/bash
#BSUB -q mpi
#BSUB -J amrvac_aia193
#BSUB -n 64
#BSUB -R "span[ptile=64]"
#BSUB -W 7200
#BSUB -o logs/job_%J.out
#BSUB -e logs/job_%J.err

cd /path/to/workdir
mpirun -np 64 ./amrvac -i euv193.par
```

提交和检查：

```bash
bsub < submit_amrvac_193_y.lsf
bjobs JOB_ID
bjobs -l JOB_ID
bpeek JOB_ID
```

不要把一个节点拆成 32 核或其他部分，也不要申请 64 核后启动不同数量的 MPI rank。

### 7. 批量处理快照

批处理时，每个快照应使用独立的输出前缀。AMRVAC 输出文件还会自动追加四位编号，
所以不要假设输出一定叫 `0000.vti`。推荐每次运行后用输入标签查找对应文件：

```bash
for number in $(seq 0 30); do
    tag=$(printf 'eru%04d' "$number")
    input="/path/to/datamr_eru/${tag}.dat"
    par="params/${tag}.par"
    prefix="/path/to/vti/euv193_y_${tag}_"

    sed \
      -e "s|__INPUT_DAT__|$input|" \
      -e "s|__OUTPUT_PREFIX__|$prefix|" \
      euv193_y_template.par > "$par"

    mpirun -np 128 ./amrvac -i "$par" \
      > "logs/${tag}.out" 2> "logs/${tag}.err"
done
```

在一个 LSF 作业中顺序处理快照，可以重复使用同一批 128 MPI 资源，避免同时提交几十个
大型作业。应检查每个 `mpirun` 的退出状态和 VTI 是否存在，不要仅凭总作业状态判断成功。

### 8. VTI 转为 PNG

使用 Python 的 VTK 接口读取 `AIA193`：

```python
import numpy as np
import vtk
from vtk.util.numpy_support import vtk_to_numpy

reader = vtk.vtkXMLImageDataReader()
reader.SetFileName("euv193_y_eru0001_0001.vti")
reader.Update()

data = reader.GetOutput()
nx, ny, _ = data.GetDimensions()
intensity = vtk_to_numpy(
    data.GetCellData().GetArray("AIA193")
).reshape(ny - 1, nx - 1)

image_log = np.zeros_like(intensity)
mask = intensity > 0
image_log[mask] = np.log10(intensity[mask])
```

绘图时应使用 AIA193 色表，并在整个时间序列中固定 `vmin/vmax`，才能比较亮度随时间的变化。
初始帧可能物理上没有落入 193 Å 温度响应范围的辐射，例如本次 `eru0000` 在 AMRVAC
和 Python/DDA 两种方法中都为全零；这不代表程序出错。

### 9. 常见问题

| 现象 | 原因与处理 |
| --- | --- |
| `invalid reference to variable in NAMELIST` | 参数不属于当前 AMRVAC 版本，检查报错行并删除非官方变量。 |
| `Wrong LOS for synthesizing emission` | `dat_resolution=.true.` 与 LOS 不兼容；改为 `.false.` 或使用受支持的负坐标轴视线。 |
| VTI 全零 | 先检查相同快照的温度范围和其他方法；初始帧可能确实没有该波段辐射。 |
| 图像左右镜像 | LOS 正负方向或 `xI1` 基向量不一致；检查日志中的三个向量。 |
| 无法读取 `.dat` | 可执行文件与快照的 AMRVAC 版本、维数、变量排列或 `mod_usr` 不兼容。 |
| `MPI_Abort` 或内存不足 | 减少并发作业数、增加 MPI rank/内存，或先用较小快照验证参数。 |
| VTI 文件名编号与输入不同 | AMRVAC 使用快照内部编号；批处理时在前缀中加入输入标签。 |

### 10. 本次服务器案例

本次 193 Å、`[0,1,0]`、31 帧原生 AMRVAC 合成位于：

```text
/share/home/guojh/euv_wave/amrvac_radiation_193_thin_y/
├── euv193_y_template.par
├── submit_amrvac_193_y.lsf
├── render_vti_aia193.py
├── params/
├── vti/
├── png/
└── logs/
```

该目录保留了可直接复用的参数模板、128 MPI 进程作业脚本、VTI 绘图脚本和逐帧日志。

### 11. 球坐标拉伸网格双视角入口

球坐标拉伸 AMR 快照使用：

```text
synthesize_aia_amrvac_native_spherical_stretched.sh
render_amrvac_native_aia.py
job_amrvac_native_aia_spherical_stretched.lsf
```

运行脚本会顺序生成两个光学薄视角：

| 输出标签 | 物理视角 | `LOS_theta` | `LOS_phi` |
| --- | --- | ---: | ---: |
| `top` | rbsl14 活动区上方径向向内俯视 | 72.0 | 195.0 |
| `side` | 南向北穿过磁绳的切向视角 | 159.925625 | 204.957372 |

示例：

```bash
MPI_RANKS=64 \
INSTRUMENT_RESOLUTION_FACTOR=0.25 \
./synthesize_aia_amrvac_native_spherical_stretched.sh \
  /path/to/amrvac.par \
  /path/to/data/frame0012.dat \
  /path/to/matching/amrvac \
  /path/to/output \
  193
```

脚本会自动处理原生转换的几个易错点：

- 为过长的集群路径建立短符号链接，避免 AMRVAC 的 Fortran 命令行字符串截断；
- 保留输入 DAT 文件名末尾四位编号，供 AMRVAC 解析输出帧号；
- 链接匹配可执行文件需要的 `path.dat` 等运行时输入；
- 转换时关闭 `do_rbsl` 和 `liakh_mass_loading`，避免对既有快照重复构造磁绳或质量源；
- 使用 `dat_resolution=.false.`，因为 AMRVAC 的球坐标 EUV 合成不支持数据分辨率快捷路径；
- 只导入无界面的 VTK XML 模块，因此可在没有 `libGLX.so.0` 的计算节点读取 VTI。

已有两张匹配 VTI、只需重新检查或调整绘图时，可以显式设置
`REUSE_EXISTING_VTI=true`；默认值是 `false`，不会误用旧结果。

rbsl14 的 `r14_g14lc_0012.dat` 已用这套入口验证成功。测试采用 AIA 193、64 MPI rank、
`instrument_resolution_factor=0.25`，原生输出为两张 `1200 x 1200` VTI；俯视和侧视
积分分别耗时约 47.5 s 和 24.3 s，两个 AMRVAC stderr 均为空。渲染脚本会检查有限值、
正强度像素和 `AIA193` CellData，并生成两张单视角 PNG、双视角 PNG 与统计 JSON。

## 脚本说明

| 脚本 | 用途 |
| --- | --- |
| `synthesize_aia_amrvac_native_spherical_stretched.sh` | 调用 AMRVAC 原生模块，对球坐标拉伸 AMR 快照生成双视角光学薄 AIA VTI。 |
| `render_amrvac_native_aia.py` | 无界面读取原生 VTI，验证 `AIA*` 数组并生成单图、双视角图和 JSON。 |
| `job_amrvac_native_aia_spherical_stretched.lsf` | guoguo 上按完整 64 核节点运行原生双视角合成的作业模板。 |
| `synthesize_radiation_cartesian.py` | 批量处理多个笛卡尔 AMRVAC 快照，支持光学薄合成和带吸收的 DDA 合成。 |
| `synthesize_aia_optically_thick_cartesian.py` | 对单个笛卡尔 AMRVAC 快照进行 AIA 171/211/304 光学厚辐射传输。 |
| `synthesize_aia_spherical_stretched.py` | 对球坐标拉伸 AMRVAC 快照进行双视角 AIA 304 光学薄/暗云转移对比。 |

球坐标拉伸网格入口可用 `--volume-output <目录>` 同时保留均匀笛卡尔三维辐射体。
该目录包含三个坐标轴以及 `hydrogen_density_cm3.npy`、`temperature_k.npy`、
`aia304_emissivity.npy`、`opacity_304_cm-1.npy`、有效域和冷云掩膜，并用
`volume_metadata.json` 记录网格范围、时间和文件映射。浮点体数据使用 `float32`，
可在不重新读取大体积 AMRVAC DAT 的情况下更换视角或辐射转移设置。

这些脚本均依赖原作者的 `radsyn_tools`：
<https://github.com/fuwentai/radsyn_tools>。本仓库不上传 `Packages/`，
首次使用时在仓库根目录执行：

```bash
bash bootstrap_dependencies.sh
```

拉取后目录中应包含 `RST_cart/` 和 `AIA_response/`，默认使用：

```text
Packages/radsyn_tools/DDA_method1
```

整个 `solar_code` 移动到其他机器后不需要修改绝对路径。需要使用另一份辐射合成工具时，
仍可通过 `--dda-root` 覆盖默认目录。

## 批量辐射合成

`synthesize_radiation_cartesian.py` 通过 glob 读取一组笛卡尔 AMRVAC `.dat`，将密度和能量转换为
数密度、压强和温度，再针对指定 AIA 波段执行辐射合成。

```bash
cd /path/to/solar_code
PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_radiation_cartesian.py \
  --dat-glob '/Volumes/ZC/jet_mhd/*dat' \
  --output-dir '/Volumes/ZC/jet_mhd/radsyn' \
  --start 1 --stop 180 --step 3 \
  --wavelength 171 131 --save-physical --plot
```

可通过 `--normal` 设置视线方向；`--save-physical` 会保存中间物理量立方，
`--plot` 会生成 PNG 快速检查图。

## 单快照 AIA 光学厚合成

`synthesize_aia_optically_thick_cartesian.py` 直接使用当前 `simesh.amrvac` 接口读取一个
笛卡尔 AMRVAC `.dat`，计算发射率、吸收系数、总光学深度和输出强度：

```text
I = integral(exp(-tau_before) * emissivity * ds) + I0 * exp(-tau_total)
```

示例：

```bash
PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_aia_optically_thick_cartesian.py \
  /path/to/mhd0006.dat \
  --wavelength 171 --normal -1 0 0 \
  --vmin 3 --vmax 10 --label example_171
```

未设置 `--output-dir` 时，输出写入输入 `.dat` 所在目录下的 `radiation/`。
每次运行会保存线性强度、`log10` 强度、总光学深度、PNG 快速检查图和元数据
JSON。没有显式设置 `--resolution` 时，默认使用 AMRVAC 最高有效分辨率。

## 球坐标拉伸网格 AIA 304 合成

`synthesize_aia_spherical_stretched.py` 使用
`simesh.amrvac.open_dataset(dat, grid_config=amrvac.par)` 重建真实拉伸径向坐标，
再将 `(r,theta,phi)` AMR 数据线性采样到均匀笛卡尔体。辐射积分仍沿单一物理视线：

- 径向俯视为一组相互平行的局部 `-e_r` 射线；
- 切向侧视为穿过磁绳轴顶的南向北平行射线；
- 射线与笛卡尔包围盒采用厘米单位 slab 相交和逐体素 DDA，不能调用旧的
  `DDA.walk_volume()` 默认 `max_t=1.0` 入口判断。

输入 `path.dat` 的前三列依次为经度（度）、纬度（度）和半径（Rsun）。当前物理转换
要求 DAT 包含 `rho`、内能变量 `e` 和参数 `gamma`；数密度单位为每代码密度
`1e9 cm^-3`，温度单位为每代码温度 `1e6 K`。

```bash
/opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_aia_spherical_stretched.py \
  /path/to/frame_0000.dat \
  /path/to/amrvac.par \
  /path/to/path.dat \
  /path/to/output/frame_0000 \
  --frame 0000 --resolution 160,120,80 --image-size 160
```

每帧输出光学薄强度、暗云强度和光学深度 NPY，两张双视角 PNG，以及记录视线、
网格、物理假设和统计量的 JSON。暗云模式使用近似 H I/He I 连续吸收并把冷云源函数
设为零，只是吸收上限测试，不是 He II 304 非 LTE 辐射转移。

## 检查命令行参数

```bash
PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_radiation_cartesian.py --help

PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_aia_optically_thick_cartesian.py --help

/opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_aia_spherical_stretched.py --help
```

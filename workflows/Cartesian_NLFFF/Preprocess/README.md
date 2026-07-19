# NLFFF 底边界预处理

本目录存放笛卡尔势场与 NLFFF 外推所需的 Python 预处理程序。
`prepare_nlfff_boundary.py` 读取 SHARP CEA 矢量磁图，完成分量转换、区域裁剪、
multigrid 降采样、可选的 Wiegelmann 风格预处理，并生成 Guo 版
MPI-AMRVAC 可直接读取的底边界文件。

如果输入是全日面 HMI `hmi.B_720s` 矢量磁场分段，而不是已经生成的 SHARP
CEA 磁图，应先运行同目录的 `make_cea_patch.py`，把 `field`、`inclination`、
`azimuth` 和 `disambig` 重投影成统一网格的 `Br.fits`、`Bt.fits`、`Bp.fits`。

对应的 AMRVAC 势场与 NLFFF 运行模板位于上一级目录：

```text
Cartesian_NLFFF/
```

## 独立安装与分发

`Cartesian_NLFFF` 可以整体复制到另一台机器，不需要同时复制 `solar_code/src`
或设置 `PYTHONPATH`。本项目所需的数值处理代码已经放在
`Preprocess/nlfff_core/`。

接收方进入 `Cartesian_NLFFF` 根目录后可创建独立 Python 环境：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r Preprocess/requirements.txt
```

Windows PowerShell 激活命令为 `.venv\Scripts\Activate.ps1`。Python 预处理建议
使用 Python 3.10 或更高版本。势场和 NLFFF 数值计算仍需接收方另外安装并配置
MPI-AMRVAC、MPI、Fortran 编译器和作业调度环境。

## 输入数据

### 从全日面 HMI 数据生成 CEA 分量

已有 SHARP CEA 三分量时可跳过本节。对于 `hmi.B_720s`，可先运行：

```bash
python Preprocess/make_cea_patch.py \
  --field /path/to/hmi.B_720s.TIME.field.fits \
  --inclination /path/to/hmi.B_720s.TIME.inclination.fits \
  --azimuth /path/to/hmi.B_720s.TIME.azimuth.fits \
  --disambig /path/to/hmi.B_720s.TIME.disambig.fits \
  --output-dir /path/to/cea_fits \
  --center-lon 35.0 --center-lat -20.0 \
  --width-degree 12.0 --height-degree 9.0 \
  --resolution-degree 0.03 \
  --sampling-mode sharp
```

所有时刻用于 TMF 时，CEA 中心、经纬范围、分辨率、采样模式和数组方向必须
保持一致。建议先选择一个时刻，与正式 SHARP CEA 产品比较后再批量处理。

### CEA 三分量要求

输入必须包含同一时刻、相同网格的三个二维 CEA FITS 文件：

| 文件 | SHARP 分量 | 转换后的局地笛卡尔分量 |
| --- | --- | --- |
| `Br.fits` | 径向分量 | `Bz = Br` |
| `Bt.fits` | 纬向分量，正方向向南 | `By = -Bt` |
| `Bp.fits` | 经向分量，正方向向西 | `Bx = Bp` |

推荐把三个文件放在同一目录，并用 `--cea-dir` 指定。也可以分别通过
`--br`、`--bt`、`--bp` 指定文件；单独指定的路径会覆盖 `--cea-dir` 中的
同名默认文件。程序兼容图像位于 primary HDU 或压缩图像扩展 HDU 的 FITS。

输入 FITS 头应包含 `CDELT1`，通常还应包含 `CDELT2`、`CRPIX1`、`CRPIX2`
和 `NAXIS1/2`。程序利用这些字段推导输出像素尺度和裁剪区域中心；必要时可用
命令行参数显式覆盖。

## 基本用法

在复制后的 `Cartesian_NLFFF` 根目录直接运行：

```bash
cd /path/to/Cartesian_NLFFF
python Preprocess/prepare_nlfff_boundary.py \
  --cea-dir /path/to/cea_fits \
  --output-dir /path/to/processed/nlfff_boundary \
  --x0 300 --y0 50 --nx 600 --ny 300 \
  --level 2 --quicklook
```

这个示例选取原始磁图中的半开区间
`x=[300,900)`、`y=[50,350)`，再进行 `2×2` binning，最终物理底图尺寸为
`300×150`。

查看全部参数：

```bash
python Preprocess/prepare_nlfff_boundary.py --help
```

## 裁剪、降采样与鬼点

`--x0`、`--y0` 是从 0 开始的 Python 像素坐标，`--nx`、`--ny` 是
**降采样前的物理区域尺寸**。四个参数都不包含 NLFFF 文件鬼点。

multigrid level 与降采样倍率的关系为：

```text
reducer = 2 ** (level - 1)
level=1  -> 原分辨率
level=2  -> 2×2 binning
level=3  -> 4×4 binning
```

`nx` 和 `ny` 必须能被 `reducer` 整除。程序会在原分辨率裁剪窗口外额外读取
`nghost × reducer` 个像素，然后再做 block average，使输出 NLFFF 文件默认
保留每侧 2 个文件鬼点。如果扩展区域超出 FITS 边界，缺失部分使用
`--blank-value` 填充，默认值为 `0.0 G`。

需要特别区分两类鬼点：

- `nlfff_boundary.dat` 中的横向 padding 是磁图输入文件的一部分；
- AMRVAC 计算过程中各 block 的内部 ghost cells 由程序自动维护。

例如物理底图为 `300×150` 且 `--nghost 2` 时，势场文件头尺寸为
`300×150`，NLFFF 文件头尺寸为 `304×154`。AMRVAC 的物理网格仍应设置为
`domain_nx1=300`、`domain_nx2=150`，不能使用包含文件鬼点的尺寸。

## 可选的磁场预处理

独立 NLFFF 入口默认启用 `--preprocess`，执行 Wiegelmann 风格的力、力矩和
平滑优化，并在终端定期输出迭代次数、损失变化以及各损失项。只有明确需要保留
原始观测矢量磁场时才使用 `--no-preprocess`。

| 参数 | 默认值 | 作用 |
| --- | ---: | --- |
| `--preprocess-mu3` | `0.1` | 控制磁场偏离原始观测的惩罚权重。 |
| `--preprocess-mu4` | `0.1` | 控制空间平滑项权重。 |
| `--preprocess-max-iter` | `5000` | 最大迭代次数。 |
| `--preprocess-report-interval` | `50` | 终端进度输出间隔。 |
| `--preprocess-tol` | `1e-4` | 收敛判据。 |

是否启用预处理仍应由具体科学目标决定。独立 NLFFF 默认处理；数据驱动 TMF
则单独保持 `--prep 0 --nlfff-prep 0`，避免静态优化改变观测时序一致性。
NLFFF 运行后应比较 `parameters.json` 中的优化前后诊断量，而不能只凭快速
检查图判断结果是否更适合作为边界。

## 坐标和 AMRVAC 网格参数

默认情况下，程序根据 FITS 的 `CDELT1/2` 和太阳半径计算降采样后的
`dx_cm`、`dy_cm`，并根据裁剪中心及 `CRPIX1/2` 推导 `xc_cm`、`yc_cm`。
以下参数可覆盖自动计算结果：

| 参数 | 含义 |
| --- | --- |
| `--dx-cm`、`--dy-cm` | 输出物理像素尺度，单位 cm。 |
| `--xc-cm`、`--yc-cm` | 裁剪区域中心的局地平面坐标，单位 cm。 |
| `--solar-radius-cm` | 坐标换算使用的太阳半径。 |
| `--z-min-10mm` | AMRVAC `xprobmin3`，单位为 10 Mm，默认 `0.1`。 |
| `--domain-nx3` | 三维计算域的竖直网格数；不指定时由底图尺寸推导。 |
| `--block-nx1/2/3` | 覆盖建议的 AMRVAC block 尺寸。 |

势场、NLFFF 和后续 TMF 必须使用相同的物理范围、基础网格、像素尺度与中心
坐标。程序生成的 `amrvac_parameters.txt` 是配置模板参数时的主要核对依据。

## 输出文件

每次运行会在 `--output-dir` 中生成：

| 文件 | 内容与用途 |
| --- | --- |
| `potential_boundary.dat` | 去除文件鬼点后的 `Bz`，供势场阶段读取。 |
| `nlfff_boundary.dat` | 依次包含 `Bx`、`By`、`Bz`，保留每侧文件鬼点，供 NLFFF 阶段读取。 |
| `allboundaries.dat` | Guo 前处理格式的三分量边界数据。 |
| `grid.ini` | 与输出物理网格相符的网格配置。 |
| `nd.ini` | 原裁剪和鬼点区域对应的维度配置。 |
| `amrvac_parameters.txt` | 建议的 AMRVAC domain、block 和物理范围参数。 |
| `parameters.json` | 输入、裁剪、padding、multigrid、坐标、诊断量和输出路径的完整记录。 |
| `boundary_quicklook.png` | `Bx/By/Bz` 三分量快速检查图；默认生成。 |

二进制边界文件采用本项目 Guo/AMRVAC 写入约定。文件头为：

```text
int32   nx, ny
float64 xc, yc, dx, dy
```

文件头之后是按约定顺序写入的 `float64` 磁场数组。可用 `--no-quicklook`
关闭检查图，或通过 `--vmax` 调整三分量图的对称色标范围，默认是
`[-500, 500] G`。

## 放入 AMRVAC 事件目录

预处理完成后，将两个主要边界文件分别放入复制出的事件运行目录：

```text
<event>/nlfff/potential/potential_boundary/potential_boundary.dat
<event>/nlfff/nlfff_boundary/nlfff_boundary.dat
```

然后按以下顺序运行：

1. 使用 `potential/` 模板编译并运行势场阶段；
2. 确认生成 `potential/data/lfffNNNN.dat`；
3. 确认上层 `amrvac.par` 的 `restart_from_file` 指向该势场结果；
4. 编译并运行 NLFFF 磁摩擦松弛；
5. 从磁摩擦日志和输出序列中选择真正完成松弛的最终 DAT，供后续 TMF 使用。

## 与 TMF 预处理的关系

本模块导出 `write_nlfff_from_cartesian()`。相邻的 TMF 预处理入口会调用该函数，
用连续 SHARP 序列的首个有效时刻生成势场和 NLFFF 初始边界，从而保证单帧
NLFFF 与 TMF 使用相同的分量、multigrid、像素尺度和文件写入约定。

两个项目的降采样均调用
`Preprocess.nlfff_core.cartesian.multigrid_reduce_components()`，可选磁场预处理均
调用同模块的 `preprocess_cartesian_field()`；不再分别维护不同的插值实现。

## 运行后检查

建议至少完成以下检查：

1. 三个输入 FITS 的形状、观测时刻和 WCS 信息一致；
2. `parameters.json` 中的物理区域尺寸等于 `nx/reducer × ny/reducer`；
3. `potential_boundary.dat` 与 `nlfff_boundary.dat` 的文件头仅相差两侧鬼点；
4. 快速检查图中的 `Bx/By/Bz` 方向和活动区位置符合预期；
5. `dx_cm`、`dy_cm`、`xc_cm`、`yc_cm` 与 AMRVAC 参数文件一致；
6. potential、NLFFF 与 TMF 使用完全相同的物理三维网格。

常见报错包括：`nx/ny` 不能被 multigrid 倍率整除、缺少某个分量 FITS、
FITS 中不存在二维图像，或缺少自动计算像素尺度所需的 `CDELT1`。遇到问题时，
优先检查终端错误、`parameters.json` 和输入 FITS 头信息。

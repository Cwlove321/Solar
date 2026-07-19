# 数据驱动 TMF（磁摩擦）模板

本目录是 MPI-AMRVAC 的三维笛卡尔数据驱动磁摩擦（TMF，time-dependent magnetofrictional）计算模板。它以 NLFFF 外推的三维磁场作为初态，并用连续 SHARP 矢量磁图推导出的底边界磁场与速度场驱动系统演化。

目录仅保存**可复用代码和参数模板**。每个活动区/爆发事件都应复制为独立运行目录；不要把可执行文件、DAT、VTU 或历史日志复制进模板。

## 文件说明

| 文件/目录 | 作用 |
| --- | --- |
| `mod_usr.t` | TMF 用户模块：读入 B/V 时间序列，设置磁摩擦物理、底边界插值和侧/顶边界。 |
| `amrvac.par` | 正式使用的 AMRVAC 参数文件。 |
| `amrvac_note.par` | 与 `amrvac.par` 对应的逐行中文注释版本，仅用于理解和修改参数。 |
| `job` | Guoguo LSF 提交模板，默认 5 个节点 × 64 核 = 320 核。 |
| `convert.par`、`convert_selected.job` | 选择性 DAT→VTU 转换模板。 |
| `Preprocess/prepare_hmi_b720_timeseries.py` | 将连续全日面 HMI `B_720s` 分段批量转换为统一的 SHARP-like CEA `Br/Bt/Bp`。 |
| `Preprocess/prepare_tmf_boundaries.py` | 将连续 CEA 磁图转换为 TMF 磁场边界与 DAVE4VM 速度边界。 |
| `Preprocess/nlfff_core/`、`Preprocess/tmf_core/` | 随项目打包的 NLFFF 边界、SHARP 和 DAVE4VM 数值实现。 |
| `data/` | 放置事件的 `B_boundary.dat`、`V_boundary.dat` 与运行输出。 |
| `nlfff/data/` | 放置用于重启的最终 NLFFF DAT 文件。 |
| `log/` | 标准输出和标准错误日志。 |

## 每个事件的推荐目录结构

```text
case_YYYYMMDD_HHMMUTC/
├── amrvac.par
├── mod_usr.t
├── job
├── data/
│   ├── B_boundary.dat
│   ├── V_boundary.dat
│   ├── driven0000.dat          # 计算后生成
│   └── driven0000.vtu          # autoconvert 开启后生成
├── nlfff/data/
│   └── nlfff0000.dat           # 最终 NLFFF 状态
└── log/
```

## 前处理和输入数据

### Python 预处理程序的简单用法

本目录可以作为独立项目整体复制。接收方先安装随包依赖：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r Preprocess/requirements.txt
```

然后在 `Cartesian_TMF` 根目录运行 `Preprocess/prepare_tmf_boundaries.py`。
输入目录中应包含按
时间排列且相互配对的 SHARP CEA `*.Br.fits`、`*.Bt.fits`、`*.Bp.fits`。
例如裁剪 `x=[300,900)`、`y=[50,350)` 并做 bin2：

```bash
cd /path/to/Cartesian_TMF
python Preprocess/prepare_tmf_boundaries.py \
  --input-dir /path/to/sharp_time_series \
  --output-dir /path/to/processed/event \
  --x1 300 --x2 900 --y1 50 --y2 350 \
  --level 2 --nlfff-level 2 \
  --prep 0 --nlfff-prep 0 \
  --driven-fast 12 --window-size 23 --threshold 100
```

程序一次生成连续驱动所需的 `B_boundary.dat`、`V_boundary.dat`，以及首个驱动
时刻对应的 `potential_boundary.dat`、`nlfff_boundary.dat`。将 B/V 文件复制到
TMF 事件目录的 `data/`，将势场和 NLFFF 边界放入相应阶段目录；先完成势场和
NLFFF，再把最终 `nlfffNNNN.dat` 作为 TMF 的 restart。默认不做 Wiegelmann
力矩预处理；需要检查 DAVE4VM 逐帧结果时再加 `--diagnostic-plots`。完整参数可
运行：

```bash
python Preprocess/prepare_tmf_boundaries.py --help
```

### 1. 生成 TMF 边界

使用 `Preprocess/driven_boundary.py` 中的 `write_driven_boundary()` 对连续 SHARP
`Br`、`Bt`、`Bp` 磁图完成裁剪、binning、预处理和 DAVE4VM 速度反演。模板的坐标约定为：

降采样与可选磁场预处理使用随包提供的 `Preprocess/nlfff_core/cartesian.py`，
与独立 `Cartesian_NLFFF` 包使用同一套分块面积平均和 Wiegelmann 风格实现。

```text
Bx = Bp
By = -Bt
Bz = Br
```

`B_boundary.dat` 和 `V_boundary.dat` 均为 stream binary 文件，内容依次为：

```text
int32   nx, ny, nt
float64 dx, dy, dz                 # km
float64 vector[nt, ny, nx, 3]
float64 time[nt]                   # s，相对 TMF 起始帧
```

磁场单位为 G，速度单位为 km/s，时间单位为 s。用户模块内固定的 `driven_fast_times=12` 会将观测演化加快 12 倍。

注意：NLFFF/势场初始边界使用 TMF 时间序列的第 1 帧；按目前
`Preprocess/driven_boundary.py` 的实现，它对应原始下载序列的第 2 张磁图
（第 1 张仅用于计算第一段 DAVE4VM 速度）。

### 2. 准备 NLFFF 重启文件

先对同一裁剪范围完成势场和 NLFFF 松弛，并将最终输出放入 `nlfff/data/`。在 `amrvac.par` 中设置：

```fortran
restart_from_file='nlfff/data/nlfff0000.dat'
firstprocess=.true.
```

`firstprocess=.true.` 配合 `usr_transform_w`，将 NLFFF 的磁场变量转换为 TMF 使用的变量布局；NLFFF 和 TMF 的网格范围、基础分辨率、AMR 设置必须一致。

## 网格与物理范围

模板示例对应裁剪后 bin2 的 `300 × 150` 底图，三维基础网格为：

```text
domain_nx1 = 300
domain_nx2 = 150
domain_nx3 = 150
block_nx* = 10
refine_max_level = 1
```

若 bin2 后像素尺寸为 \(\Delta x=\Delta y=732.825\) km、磁图中心取 `(xc, yc)=(0,0)`，则：

\[
x_{min/max}=x_c\mp N_x\Delta x/2,\qquad
y_{min/max}=y_c\mp N_y\Delta y/2.
\]

AMRVAC 使用 `unit_length=10^9 cm`（10 Mm），所以该例的横向范围为：

```text
x: -10.992375 到 10.992375
y:  -5.4961875 到  5.4961875
z:   0.1       到 10.992375
```

`nlfff_boundary.dat` 会比物理磁图四周各多两格鬼点，例如 `300×150` 的物理磁图对应 `304×154` 的 NLFFF 输入。鬼点不计入上述 TMF 物理计算盒子。

## 输出时间和结束时间

`dtsave_dat` 的单位是归一化后的 Alfvén 时间，而不是物理秒。对于此模板的单位和加速因子 12：

```text
每物理 1 小时输出一次：dtsave_dat = 3.4934654d0
```

一般公式为：

\[
dtsave\_dat = \frac{\Delta t_{physical}}{unit\_time\times driven\_fast\_times}.
\]

`time_max` 必须小于最后一帧驱动边界的时间。当前底边界插值只在相邻两帧之间有定义；若推进到或越过最后一个时间戳，会报：

```text
Special boundary is not defined for this time
```

应留出至少一个数值时间步的裕量，或在 `mod_usr.t` 中将最后一帧边界钳制为常量。

## DAT、VTU 与编号

默认 `autoconvert=.false.`，只生成 DAT。若希望每个 DAT 同时转换为 VTU，设为：

```fortran
autoconvert=.true.
convert_type='vtuBmpi'
```

从 NLFFF restart 后，AMRVAC 往往继承 restart 文件内部的快照编号；若希望本次 TMF 输出从 `driven0000` 开始，需在 `&filelist` 内显式加入：

```fortran
snapshotnext=0
```

`w_write(i)=.false.` 仅控制第 `i` 个基础状态变量是否写入输出，不影响求解本身、边界条件或 restart。

`aiconvert n` 中的 `n` 是 DAT 文件编号，而非物理第 n 小时。它实际调用：

```bash
./amrvac -i amrvac.par -if data/drivenNNNN.dat -convert
```

本模板的 `mod_usr.t` 已在 `convert` 模式跳过完整 B/V 时序边界的读取和时间驱动边界条件，避免单帧转换重复加载全时序数据。

## 编译、提交和队列规则

每次修改 `mod_usr.t` 后必须重新生成和编译：

```bash
export AMRVAC_DIR=/share/home/guoyang1/codes/amrvac
$AMRVAC_DIR/setup.pl -d=3
make -j 8
```

正常 TMF 计算提交：

```bash
bsub < job
```

Guoguo 队列约束：

```text
mpi    队列只接受 64 的整数倍核心数，例如 64、128、320。
serial 队列用于单核或小规模串行任务。
```

因此，完整 TMF 可使用 `#BSUB -n 320` 与 `#BSUB -R "span[ptile=64]"`；单帧 `aiconvert` 可使用 `#BSUB -n 1`、`#BSUB -q serial`。所有 MPI 运行均应通过 LSF 提交，不在登录节点直接执行 `mpirun`。

## 推荐运行流程

1. 复制模板到按事件时间命名的新目录。
2. 写入/复制 `B_boundary.dat`、`V_boundary.dat` 和最终 NLFFF DAT。
3. 依据边界文件核对 `domain_nx*`、`xprobmin/max*`、AMR、`time_max`、`dtsave_dat`。
4. 在新的输出目录运行编译命令。
5. 先检查 `amrvac.par` 中的 `restart_from_file`、`firstprocess`、磁场尺度、输出设置和快照编号。
6. 提交 `bsub < job>`，用 `bjobs`、`log/*.out` 和 `data/driven*.dat` 监测。
7. 若需手动转换，使用批处理脚本调用 `aiconvert amrvac.par N`；不要在登录节点直接运行。

## 常见问题

| 现象 | 原因与处理 |
| --- | --- |
| 输出从 `driven0004.dat` 开始 | restart 文件保留了快照号；在 `filelist` 加 `snapshotnext=0` 可重置。 |
| 最后一段退出并提示边界未定义 | `time_max` 到达最后边界时刻；提前终止或钳制末帧边界。 |
| `aiconvert` 很慢 | 转换会读完整 3D DAT 并写更大的 VTU；使用已启用 `convert` 保护的用户模块，且通过批处理队列运行。 |
| `aiconvert` 段错误 | 通常是转换模式跳过 B/V 读入后仍调用数据驱动 `specialbound_usr`；需同样在该子程序中对 `convert` 直接返回。 |
| `mpi` 队列拒绝单核任务 | 改投 `serial` 队列，或申请 64 的整数倍核心数。 |

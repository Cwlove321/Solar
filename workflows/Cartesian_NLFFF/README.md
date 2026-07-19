# Cartesian 势场与 NLFFF 模板

本目录是 Guo 版 MPI-AMRVAC 的两阶段笛卡尔磁场外推模板：先由观测 `Bz`
构造势场，再以观测矢量磁场为底边界进行磁摩擦 NLFFF 松弛。每个活动区
必须复制到独立事件目录后运行；模板不保存边界数据、可执行文件、DAT、VTU
或旧日志。

## 两阶段流程

```text
SHARP Br/Bt/Bp 时序
       │
       ├─ 首个 TMF 时刻：Bz ─────────→ potential_boundary.dat
       │                                 │
       │                                 ▼
       │                           potential/lfff0000.dat
       │                                 │ restart
       └─ 首个 TMF 时刻：Bx/By/Bz ─→ nlfff_boundary.dat
                                         │
                                         ▼
                                  NLFFF/nlfff0000.dat
                                         │
                                         ▼
                                 数据驱动 TMF 初态
```

局地笛卡尔分量约定固定为：

```text
Bx = Bp
By = -Bt
Bz = Br
```

## 推荐事件目录结构

```text
case_YYYYMMDD_HHMMUTC/
├── nlfff/
│   ├── amrvac.par
│   ├── amrvac_note.par
│   ├── mod_usr.t
│   ├── job
│   ├── nlfff_boundary/nlfff_boundary.dat
│   ├── potential/
│   │   ├── amrvac.par
│   │   ├── amrvac_note.par
│   │   ├── mod_usr.t
│   │   ├── job
│   │   └── potential_boundary/potential_boundary.dat
│   ├── data/                  # NLFFF 输出
│   └── log/
```

先运行 `potential/`，确认 `potential/data/lfff0000.dat` 已生成后，再运行上层
`nlfff/`。NLFFF 的 `restart_from_file` 必须指向该势场输出。

## Python 预处理程序的简单用法

本目录可作为独立项目整体复制。接收方先在 `Cartesian_NLFFF` 根目录安装依赖：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install -r Preprocess/requirements.txt
```

然后直接运行 `Preprocess/prepare_nlfff_boundary.py`。输入可以是包含
`Br.fits`、`Bt.fits`、`Bp.fits` 的 CEA 目录，也可以分别用 `--br`、`--bt`、
`--bp` 指定三个文件。例如将原图的 `x=[300,900)`、`y=[50,350)` 区域做
bin2：

```bash
cd /path/to/Cartesian_NLFFF
python Preprocess/prepare_nlfff_boundary.py \
  --cea-dir /path/to/cea_fits \
  --output-dir /path/to/processed/nlfff_boundary \
  --x0 300 --y0 50 --nx 600 --ny 300 --level 2 \
  --quicklook
```

其中 `x0/y0` 是裁剪起点，`nx/ny` 是 **binning 前**的裁剪宽度和高度；
`level=2` 表示 `2×2` binning。程序会生成
`potential_boundary.dat`、`nlfff_boundary.dat`、参数摘要和快速检查图。
将两个 DAT 分别复制到事件模板的
`potential/potential_boundary/` 和 `nlfff_boundary/`，再按“势场 → NLFFF”顺序
编译并提交。完整参数可运行：

```bash
python Preprocess/prepare_nlfff_boundary.py --help
```

## 边界文件与鬼点

`Preprocess/prepare_nlfff_boundary.py` 的边界写入流程会生成两个主要二进制文件：

| 文件 | 内容 | 尺寸示例 | 鬼点 |
| --- | --- | --- | --- |
| `potential_boundary.dat` | 仅 Bz，用于势场 | `300×150` | 无 |
| `nlfff_boundary.dat` | Bx、By、Bz，用于 NLFFF | `304×154` | 四周各 2 格 |

文件头顺序是：`int32 nx, ny`、`float64 xc, yc, dx, dy`，随后为磁场数组。
势场文件仅跟随一个 Bz 数组；NLFFF 文件依次写 Bx、By、Bz。

NLFFF 输入的两层鬼点是**文件层面**的横向 padding，用于观测磁图的边缘处理；
它不同于 AMRVAC 在计算时为每个 block 自动分配的内部 ghost cells。物理
模拟盒子仍由 `domain_nx1/2/3` 与 `xprobmin/max*` 决定，不应把 `304×154`
误设为物理网格大小。

## 网格和坐标范围

势场、NLFFF 与后续 TMF 的物理底边界必须一致。对于原始裁剪
`x=[300,900)`、`y=[50,350)` 后再 bin2 的事件：

```text
物理底图            300 × 150
建议三维基础网格      300 × 150 × 150
block                10 × 10 × 10
像素间距             732.825 km
```

当磁图中心设为零、`unit_length=10^9 cm` 时，物理范围为：

```text
xprobmin1 = -10.992375d0      xprobmax1 = 10.992375d0
xprobmin2 =  -5.4961875d0     xprobmax2 =  5.4961875d0
xprobmin3 =   0.1d0           xprobmax3 = 10.992375d0
```

模板中保留的 `300×200×200`、旧范围只是历史示例值；运行新事件前必须同时
修改 potential 与 NLFFF 两个 `amrvac.par`，并与 TMF 参数保持一致。

## 势场阶段（potential）

`potential/mod_usr.t` 读取：

```text
potential_boundary/potential_boundary.dat
```

并写出：

```text
data/lfff0000.dat
data/lfff0000.vtu       # autoconvert=.true. 时
```

该阶段只用观测 `Bz` 产生势场初态。完成后检查输出 DAT 是否存在、网格是否
匹配，以及日志中是否有边界文件或 MPI 错误。

## NLFFF 阶段

`mod_usr.t` 读取：

```text
nlfff_boundary/nlfff_boundary.dat
```

且 `amrvac.par` 从：

```fortran
restart_from_file='potential/data/lfff0000.dat'
```

开始。底部的 Bx、By、Bz 被固定为输入磁图；侧面/顶部用用户边界条件与散度
控制处理。`mf_it_max=60000` 是磁摩擦最大迭代步数，`mf_ditsave=20000` 控制
磁摩擦诊断保存间隔。

典型输出：

```text
data/nlfff0000.dat
data/nlfff0000.vtu
data/nlfff_mflog.csv
```

最终用于 TMF 的 restart 应选择完成松弛后的最后一个 `nlfffNNNN.dat`，而非
仅按文件名假设 `0000` 一定是最终状态。

## 参数文件说明

- `amrvac.par` 与 `potential/amrvac.par` 是实际运行参数。
- `amrvac_note.par` 与 `potential/amrvac_note.par` 是逐行中文注释版，用于
  查阅参数含义；修改后应把需要的值同步回正式 `amrvac.par`。
- 这套历史模板含有少数位于 namelist 结束 `/` 之后的旧式重复项，例如
  `restart_from_file`、`dtsave_dat`、`time_max`。新的事件参数应将实际生效值
  写在对应 namelist 内，避免依赖这些历史兼容写法。

最重要的运行前核对项：

```text
potential/NLFFF 的 domain_nx*、xprobmin/max*、block_nx*、refine_max_level
potential_boundary 与 nlfff_boundary 的文件头、像素尺度和中心坐标
NLFFF restart_from_file 是否为本事件的 lfff 输出
time_max、dtsave_dat、autoconvert、snapshotnext
```

## 编译与提交

每次修改任一 `mod_usr.t` 后，在对应阶段目录重新编译：

```bash
export AMRVAC_DIR=/share/home/guoyang1/codes/amrvac
$AMRVAC_DIR/setup.pl -d=3
make -j 8
```

然后使用 LSF：

```bash
bsub < job
bjobs
```

Guoguo 的 `mpi` 队列只接受 64 的整数倍核心数；势场和 NLFFF 的模板 job
默认使用 128 核。小型单核后处理应提交 `serial` 队列，不在登录节点直接执行
`mpirun`。

## 常见问题

| 现象 | 原因与处理 |
| --- | --- |
| NLFFF 找不到边界文件 | 运行目录中缺少 `nlfff_boundary/nlfff_boundary.dat`，或相对路径不匹配。 |
| 势场/NLFFF 网格不一致 | 两阶段的网格参数没有一起更新；应以 `potential_boundary.dat` 的物理尺寸为准。 |
| 把 `304×154` 当作物理网格 | 这包含每侧两层输入鬼点；物理区域仍为 `300×150`。 |
| NLFFF restart 失败 | `potential/data/lfff0000.dat` 不存在、使用了不同网格，或没有先运行势场阶段。 |
| 输出编号不是从 0000 开始 | restart 继承了快照号；在 `&filelist` 显式设 `snapshotnext=0`。 |
| 自动 VTU 过慢 | `autoconvert=.true.` 会在每个 DAT 后写较大的 VTU；可先只存 DAT，后续批处理转换。 |

## 诊断工具

`Preprocess/diagnostics/` 保存可重用的检查脚本，诊断结果应写入对应案例的
`products/`，不写回工作流目录。

```bash
# 比较 make_cea_patch 与官方 SHARP CEA 的 Br/Bt/Bp
PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Cartesian_NLFFF/Preprocess/diagnostics/compare_cea_patch_to_sharp.py \
  --patch-dir /path/to/generated_cea \
  --official-dir /path/to/official_sharp \
  --output-dir cases/<event>/products/cea_comparison

# 绘制 NLFFF 磁摩擦日志中的无力与散度指标
PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Cartesian_NLFFF/Preprocess/diagnostics/plot_nlfff_forcefree_divergence.py \
  --csv /path/to/nlfff_mflog.csv \
  --out cases/<event>/products/figures/nlfff_forcefree_divergence.png
```

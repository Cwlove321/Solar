# TMF 时间序列预处理

本目录包含数据驱动 Cartesian TMF 所需的完整 Python 预处理流程，可以随
`Cartesian_TMF` 整体复制到另一台机器运行，不需要同时复制 `solar_code/src`、
`Cartesian_NLFFF` 或旧的 `solar_package/tmf_prep`。

## 独立安装

在复制后的 `Cartesian_TMF` 根目录创建 Python 环境：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r Preprocess/requirements.txt
```

Windows PowerShell 激活命令为 `.venv\Scripts\Activate.ps1`。建议使用
Python 3.10 或更高版本。真正编译和运行 TMF 仍需另外安装 MPI-AMRVAC、MPI、
Fortran 编译器和对应的作业调度环境。

## 目录内容

| 文件或目录 | 作用 |
| --- | --- |
| `prepare_hmi_b720_timeseries.py` | 把连续全日面 HMI `B_720s` 分段批量转换为统一网格的 CEA `Br/Bt/Bp`。 |
| `make_cea_patch.py` | 单时刻 HMI `field/inclination/azimuth/disambig` 到 CEA 三分量的转换入口。 |
| `prepare_tmf_boundaries.py` | 连续 CEA 三分量到 TMF `B/V_boundary.dat` 的主入口。 |
| `prepare_nlfff_boundary.py` | 为 TMF 初始时刻生成 potential 和 NLFFF 边界。 |
| `driven_boundary.py` | 时序组织、DAVE4VM 反演和 TMF 二进制文件写入。 |
| `nlfff_core/` | 内置 CEA、block-average、Wiegelmann 预处理和 NLFFF 边界写入代码。 |
| `tmf_core/dave4vm/` | 内置 DAVE4VM 求解器、矩阵和导数代码。 |
| `tmf_core/sharp_boundary.py` | 内置 SHARP CEA 裁剪、分量转换和时间读取。 |
| `requirements.txt` | 独立运行需要安装的第三方 Python 包。 |

## 完整处理流程

```text
HMI B_720s field/inclination/azimuth/disambig 时间序列
        ↓ prepare_hmi_b720_timeseries.py
统一范围与分辨率的 SHARP-like CEA Br/Bt/Bp 时间序列
        ↓ prepare_tmf_boundaries.py
        ├─ B_boundary.dat
        ├─ V_boundary.dat
        ├─ UTC.npy
        ├─ potential_boundary.dat
        └─ nlfff_boundary.dat
```

如果输入已经是连续的 SHARP CEA `Br/Bt/Bp`，可以跳过第一步。

## 第一步：批量生成 CEA 时间序列

输入文件名应符合以下形式，并保证每个时刻四个分段完整：

```text
hmi.b_720s.YYYYMMDD_HHMMSS_TAI.field.fits
hmi.b_720s.YYYYMMDD_HHMMSS_TAI.inclination.fits
hmi.b_720s.YYYYMMDD_HHMMSS_TAI.azimuth.fits
hmi.b_720s.YYYYMMDD_HHMMSS_TAI.disambig.fits
```

查看参数：

```bash
python Preprocess/prepare_hmi_b720_timeseries.py --help
```

示例：

```bash
python Preprocess/prepare_hmi_b720_timeseries.py \
  --input-dir /path/to/hmi_b720_segments \
  --output-dir /path/to/cea_timeseries \
  --harpnum 12345 \
  --center-lon 35.0 --center-lat -20.0 \
  --width-degree 12.0 --height-degree 9.0 \
  --resolution-degree 0.03 --sampling-mode sharp
```

程序会逐时刻调用同目录的 `make_cea_patch.py`，输出按时间命名的
`*.Br.fits`、`*.Bt.fits`、`*.Bp.fits`，并把每个时刻的转换参数保存到
`metadata/`。已完整生成的时刻会自动跳过。

## 第二步：生成 TMF 与初始 NLFFF 边界

```bash
python Preprocess/prepare_tmf_boundaries.py \
  --input-dir /path/to/cea_timeseries \
  --output-dir /path/to/processed/event \
  --x1 300 --x2 900 --y1 50 --y2 350 \
  --level 2 --nlfff-level 2 \
  --prep 0 --nlfff-prep 0 \
  --driven-fast 12 --window-size 23 --threshold 100
```

`x1/x2/y1/y2` 是降采样前、从 0 开始、末端不包含的 Python 像素范围。
`level=2` 表示 `2×2` block average，`level=3` 表示 `4×4` block average。

默认保持：

```text
--prep 0
--nlfff-prep 0
```

也就是连续 TMF 磁场和初始 NLFFF 边界均不做 Wiegelmann 静态预处理，避免改变
观测时序。降采样使用随包提供的 `nlfff_core.cartesian` block-average 实现，
与独立 `Cartesian_NLFFF` 包使用的算法一致。

## 输出目录

主入口会在 `--output-dir` 下生成：

```text
<output>/
├── tmf/bd/
│   ├── B_boundary.dat
│   ├── V_boundary.dat
│   ├── UTC.npy
│   ├── potential_boundary.dat
│   └── nlfff_boundary.dat
├── fig/v_fig/                  # 仅 --diagnostic-plots 时保留逐帧图
└── tmf_preprocess_metadata.json
```

`B_boundary.dat` 和 `V_boundary.dat` 的文件头依次为 `int32 nx,ny,nt` 和
`float64 dx,dy,dz`，随后是三分量时序数组与模拟时间。磁场单位为 G，速度单位
为 km/s，时间单位为 s。

## 独立性说明

本目录运行时只使用相对导入和随包代码，不依赖：

- `solar_code/src/solarutils`；
- 相邻的 `Cartesian_NLFFF`；
- `/Users/.../solar_package/tmf_prep`；
- `PYTHONPATH=src`。

第三方依赖统一记录在 `requirements.txt`。发送前建议保留整个
`Cartesian_TMF` 目录结构，不要只复制单个入口脚本。

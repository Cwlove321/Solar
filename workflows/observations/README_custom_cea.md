# 自定义 SHARP 类磁场与 AIA 共用 CEA 网格

磁场生成入口已接入 `Packages/sharp/cea_pipeline.py` 的官方 JSOC 数值内核。
`Cartesian_NLFFF/Preprocess/make_cea_patch.py` 与
`Cartesian_TMF/Preprocess/make_cea_patch.py` 共用这个实现。
旧的 spline + Gaussian + block average 近似已从这两个入口移除。

## 构建和输入

本机使用 `/opt/anaconda3/bin/python`。首次使用或 Python 架构改变后运行：

```bash
/opt/anaconda3/bin/python /Users/xiaoxingyun/Desktop/solar_code/Packages/sharp/build_core.py
```

当前编译适配为 macOS Clang + Accelerate；尚未验证 Linux/Windows 构建。
复制工作流时，需要保留 `solar_code/Packages/sharp` 与
`solar_code/workflows` 的相对位置，包括 upstream 源码和 Wiener 数据表。
旧的已生成 SHARP 文件作为边界输入时，不需要该编译内核。

输入为同一记录、相同几何的四个全日面 HMI `B_720s` 分段：
`field`、`inclination`、`azimuth`、`disambig`。默认使用消歧 bit 2。
先转换为观测坐标系矢量，再执行官方 Wiener6 插值、3 倍过采样、
sigma=1/半宽=2 的官方 Gaussian 降采样，最后转换为 Br、Bt、Bp。
Bt 正方向为南，Bp 正方向为西。该流程生成自定义 SHARP 类产品，
不包含官方 HARP 自动识别、边框追踪或全套误差/质量产品。

## 示例：2024-10-09 00:00 TAI

下列范围单位是角秒，顺序为 `[Solar-X_min, Solar-X_max, Solar-Y_min, Solar-Y_max]`。
从仓库根目录运行：

```bash
/opt/anaconda3/bin/python workflows/Cartesian_NLFFF/Preprocess/make_cea_patch.py \
  --field /Volumes/Seagate/w3_sphere_rbsl/data/hmi/m_720/hmi.B_720s.20241009_000000_TAI.field.fits \
  --inclination /Volumes/Seagate/w3_sphere_rbsl/data/hmi/m_720/hmi.B_720s.20241009_000000_TAI.inclination.fits \
  --azimuth /Volumes/Seagate/w3_sphere_rbsl/data/hmi/m_720/hmi.B_720s.20241009_000000_TAI.azimuth.fits \
  --disambig /Volumes/Seagate/w3_sphere_rbsl/data/hmi/m_720/hmi.B_720s.20241009_000000_TAI.disambig.fits \
  --bounds 100 800 150 500 --resolution-degree 0.03 \
  --components-only --output-dir products/custom_cea

/opt/anaconda3/bin/python workflows/observations/project_aia_to_sharp.py \
  --sharp products/custom_cea/Br.fits \
  --aia /Volumes/Seagate/w3_sphere_rbsl/data/aia/193/data/aia.lev1_euv_12s.2024-10-08T235918Z.193.image_lev1.fits \
  --out products/custom_cea/aia
```

磁场输出 `Br.fits`、`Bt.fits`、`Bp.fits`、`roi_mask.fits`、`parameters.json`。
默认另有 quicklook；`--no-quicklook` 可关闭。
去掉 `--components-only` 还会输出原有 NPZ、局地边界 ASCII 和最近邻采样的消歧位图。
NPZ 的 lon/lat 和磁场 FITS 现在明确使用 Carrington 坐标。

## 选择中心与范围

- `--bounds`：自动取日面矩形中点作为中心，建立覆盖其投影边界的 CEA 网格，
  矩形之外填 NaN。可同时给 `--center-lon/--center-lat` 指定另一中心。
- 原有 `--center-lon/--center-lat --width-degree/--height-degree` 用法仍可用；
  输入中心默认是 Stonyhurst 坐标，内部转换为 Carrington CEA。
- `--center-frame carrington`：输入中心直接使用 Carrington 经纬度。
- `--target-cea reference.Br.fits`：复用已有 Carrington CEA 的像素网格和有效区域，
  不再传范围或中心。目标与全日面 HMI 必须具有相同 T_OBS。

HPC 矩形的边界在 CEA 中通常是曲线，输出数组为包围它的矩形。
用于 NLFFF/TMF 时，应从有效区域内选择完整矩形，不要把掩膜外的 NaN 当作零磁场。
分辨率、中心和网格在每一组磁场/AIA 产品中完全共用。

## AIA 输出与兼容性

`project_aia_to_sharp.py` 可用 `--aia file1.fits file2.fits` 选择文件，或保留
`--aia-dir` 批量处理。默认使用整个目标网格，支持奇数尺寸；原来的
`--x0/--x1/--y0/--y1` 像素裁剪仍可用，末端不包含。

每幅 AIA 输出全分辨率 FITS、NPY、来源与参考时间元数据及预览图。
默认按曝光时间归一化为 DN/s，已为 DN/s 的数据不会重复除曝光时间；
`--raw-dn` 保留原单位。默认 Howard 差分自转对齐到磁图 T_OBS，
保留 `--rotation-model` 和 `--no-differential-rotation` 参数。
不额外执行 aia_prep、内容配准或日冕高度校正。

FITS 中 `DATE-OBS` / `T_OBS` 表示目标参考时间，`SRC_TIME` 是 AIA 原始观测中点，
`SRC_EXPT` 是原始曝光时间。磁场保留 HMI 仪器与 T_OBS 信息，避免重读后
错误地采用曝光开始时间。AIA 使用磁图的有效像素掩膜；自身无覆盖位置仍为 NaN。

旧版默认生成 `_bin2.npy`，新版需要显式传 `--bin2`，且要求偶数尺寸，
不会为了降采样丢掉奇数边缘。旧版 NPY 为 DN，新版默认 DN/s，下游强度阈值应匹配单位。
已有 Python 调用的 `reproject_to_sharp()` 保留原来的低层接口；新的
`project_aia_product()` 负责单位、掩膜和可重读 FITS 的完整处理。

两个时间序列入口也使用官方内核，恢复计算检查包含后端版本。
原有近似结果不会被当作新版本结果跳过。NLFFF 入口继续采用中心差分自转追踪；
TMF 入口继续使用其原有固定 Stonyhurst 中心策略。

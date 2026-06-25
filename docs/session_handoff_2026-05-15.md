# Session Handoff: Observation Workflow and AIA Preprocessing

日期：2026-05-15  
项目路径：`/Users/xiaoxingyun/Desktop/solar_code`

这份文档用于下一个对话快速恢复上下文。下一次可以直接说：

```text
请读取 docs/session_handoff_2026-05-15.md，然后继续上次 solar_code 项目。
```

## 已完成工作

### 1. 观测数据下载 workflow

已经搭建并测试了观测数据下载 workflow。

测试数据时间范围：

```text
2024-10-08T23:00:00 到 2024-10-09T00:00:00
```

测试数据路径：

```text
data/test/observations/20241008_2300_2400/
```

内容包括：

- SDO/AIA：94, 131, 171, 193, 211, 304, 335 Å，共 49 个 FITS。
- HMI 45s：7 个 FITS。
- HMI 720s：6 个 FITS。
- SHARP：`HARPNUM=12001`，共 30 个 FITS。
- STEREO-A EUVI 195 Å：6 个 FITS。
- Solar Orbiter EUI FSI 174 Å L2：6 个 FITS。

中文下载说明：

```text
workflows/observations/download_zh.md
```

英文下载说明入口：

```text
workflows/observations/README.md
```

### 2. AIA 子母图绘图

已经实现 AIA 子母图绘图：

- 左边母图为 full-disk。
- 右边子图为 submap。
- 母图上用 cyan 虚线框出子图区域。
- 用 cyan 虚线连接母图框和子图。
- AIA 标题格式为：

```text
AIA 171 Å 23:00:09 UT
```

- 坐标标签统一为：

```text
Solar-X (arcsec)
Solar-Y (arcsec)
```

相关基础代码：

```text
src/solarutils/plotting/maps.py
src/solarutils/plotting/__init__.py
```

视觉测试脚本：

```text
tests/visual/plot_aia_context_examples.py
```

测试图归类路径：

```text
data/test/figures/observations/sdo/aia/171/context_submap.png
data/test/figures/observations/sdo/aia/193/context_rectangle.png
```

### 3. AIA L1 到 L1.5 预处理

已经在基础代码中实现 AIA 预处理：

```python
from solarutils.obs.maps import aia_prep
```

处理步骤：

1. `update_pointing`
2. 可选 `respike`
3. 可选 PSF deconvolution
4. `register`

重要决定：

- PSF 去卷积默认关闭。
- 只有显式传入 `deconvolve_psf=True` 才会做 PSF 去卷积。

相关代码：

```text
src/solarutils/obs/maps.py
src/solarutils/obs/__init__.py
```

### 4. AIA PSF 缓存

已经生成 7 个 AIA EUV 波段的 PSF 缓存：

```text
data/test/calibration/sdo/aia/psf/
├── aia_0094_psf.npz
├── aia_0131_psf.npz
├── aia_0171_psf.npz
├── aia_0193_psf.npz
├── aia_0211_psf.npz
├── aia_0304_psf.npz
└── aia_0335_psf.npz
```

每个 PSF：

- shape: `(4096, 4096)`
- dtype: `float32`
- 文件大小约 58 MB
- 归一化总和约为 1

相关函数：

```python
from solarutils.obs.maps import (
    aia_psf_path,
    cache_aia_psfs,
    load_aia_psf,
)
```

缓存脚本：

```text
tests/visual/cache_aia_psf_examples.py
```

### 5. 与旧 notebook 的对比

已经和旧 notebook 对比过：

```text
/Users/xiaoxingyun/Desktop/20241008/cme/193.ipynb
```

旧流程核心为：

```text
update_pointing -> aiapy.psf.deconvolve -> register
```

使用同一张 AIA 193 FITS、同一个 JSOC pointing table、同一个 193 PSF、同样 5 次迭代时：

```text
max_abs_difference = 0
```

说明 `solarutils.obs.maps.aia_prep` 与旧 notebook 预处理结果完全一致。

对比脚本：

```text
tests/visual/compare_aia193_notebook_flow.py
```

对比图：

```text
data/test/figures/observations/sdo/aia/193/notebook_flow_compare.png
```

### 6. 局部 padded PSF 去卷积

已经实现局部 padded PSF 去卷积：

```python
from solarutils.obs.maps import deconvolve_submap_with_padding
```

默认参数：

```python
padding = 200.0  # arcsec
psf = None       # 可选；不传则按 aia_map.wavelength 自动计算
```

示例：

```python
from solarutils.obs.maps import deconvolve_submap_with_padding, aia_psf_path

local = deconvolve_submap_with_padding(
    aia_map,
    extent=(100, 800, 150, 500),
    psf=aia_psf_path(193),
    padding=200,
    iterations=25,
)
```

相关辅助函数：

```python
pad_extent
ensure_even_shape
crop_psf_to_shape
```

AIA 193 测试区域：

```text
x = 100 到 800 arcsec
y = 150 到 500 arcsec
```

和全日面去卷积后裁剪的对比：

```text
padding=0      mean |diff|=3.703   max=127.384
padding=100    mean |diff|=0.757   max=10.851
padding=200    mean |diff|=0.436   max=6.569
padding=400    mean |diff|=0.193   max=1.722
```

结论：

- `padding=200 arcsec` 是速度和误差的较好折中。
- 快速预览可用局部去卷积。
- 正式定量分析仍建议 full-frame PSF deconvolution 后再裁剪，或用更大的 padding。

测试脚本：

```text
tests/visual/compare_aia193_local_psf_padding.py
```

测试图：

```text
data/test/figures/observations/sdo/aia/193/local_psf_padding_compare.png
```

### 7. AIA 锐化验证

已经把旧 notebook 中的锐化程序整理为基础函数，并增加了多尺度版本：

```python
from solarutils.obs.maps import multiscale_sharpen_map, sharpen_map
```

公式为：

```text
sharpened = data + (data - gaussian_filter(data, sigma)) * weight
```

默认参数对应旧 notebook：

```python
sigma = 1.0
weight = 1.5
```

注意：这是 unsharp mask，主要用于视觉增强，不能替代 PSF 去卷积。

多尺度默认参数：

```python
sigmas = (1.0, 2.0, 4.0)
weights = (0.8, 0.35, 0.15)
```

验证脚本：

```text
tests/visual/verify_aia_sharpen.py
```

验证图：

```text
data/test/figures/observations/sdo/aia/193/sharpen_compare.png
data/test/figures/observations/sdo/aia/193/multiscale_sharpen_compare.png
```

本次 AIA 193 子图验证结果：

```text
original range: 76.000 to 2570.000
sharpened range: 49.621 to 3592.817
difference mean/std: -0.000 / 20.189
```

多尺度默认参数验证结果：

```text
original range: 76.000 to 2570.000
sharpened range: 43.323 to 3700.406
difference mean/std: -0.000 / 20.931
```

### 8. AIA derotate FITS metadata 规则

已经调整 `derotate_to_reference`：

- 输出数据尺寸裁剪回参考图尺寸，不再保留旧 notebook 中 `dy+1, dx+1` 的最终尺寸。
- 输出 map 的 WCS/几何 meta 来自参考图。
- 输出 map 的时间关键字来自当前输入图，例如 `DATE-OBS` 和 `T_OBS`。
- derotate 时数据除以当前输入图曝光时间，再乘以参考图曝光时间。
- 输出数据统一到参考曝光时间后的 `DN` 标度，`BUNIT='DN'`。
- 输出 `EXPTIME` 保持和参考图一致。
- 太阳自转计算改为显式使用当前图 `observer_coordinate`，不再触发 SunPy 的 Earth-observer warning。
- 对 float 输出清理参考 header 中不适用的 `BLANK` keyword。

验证脚本：

```text
tests/visual/verify_aia_derotate.py
```

验证输出：

```text
data/test/processed/observations/sdo/aia/193/derotate_20241008T231006_to_230006.fits
```

落盘 FITS header 验证：

```text
shape (584, 1167)
DATE-OBS 2024-10-08T23:10:04.843
T_OBS 2024-10-08T23:10:05.843
EXPTIME 2.000075
BUNIT DN
BLANK None
```

批量处理脚本：

```text
workflows/observations/process_aia_prep_derotate_submaps.py
```

运行命令：

```bash
PYTHONPATH=src python workflows/observations/process_aia_prep_derotate_submaps.py --overwrite
```

已生成测试数据库：

```text
data/test/processed/observations/20241008_2300_2400/sdo/aia/
```

内容：

- 7 个波段：0094, 0131, 0171, 0193, 0211, 0304, 0335。
- 每个波段 7 个 FITS。
- 共 49 个 FITS。
- 每张图都先 `aia_prep(deconvolve_psf=True)`，再裁剪/derotate 到该波段第一帧参考 submap。
- 默认 submap：`x=100-800 arcsec, y=150-500 arcsec`。
- 输出 shape：`(584, 1168)`。
- `BLANK` keyword 已清理，`BUNIT='DN'`。
- 总大小约 238M。

绘图脚本：

```text
workflows/observations/plot_processed_aia_submaps.py
```

运行命令：

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py --overwrite
```

已生成图片：

```text
data/test/figures/observations/20241008_2300_2400/sdo/aia/processed/
```

- 单帧 quicklook PNG：49 张。
- 每个波段一张时间序列总览 PNG：7 张。
- 总大小约 44M。
- 例图：`0193/0193_prep_psf_derotate_submap_sequence.png`。
- 绘图默认使用每个波段第一张图的 `1` 和 `99.5` 百分位作为该波段所有帧共同的 `vmin/vmax`。
- 可通过 `--norm-source first|all` 和 `--percentiles 1 99.5` 修改。
- 绘图支持 `--sharpen none|single|multiscale`。
- 如果启用锐化，`vmin/vmax` 会基于锐化后的图像计算。
- 锐化版本保存到 `processed/single_sharpen/` 或 `processed/multiscale_sharpen/`，不覆盖原始 quicklook。
- 已生成多尺度锐化版：49 张单帧 PNG + 7 张时间序列 PNG。

电影基础函数：

```python
from solarutils.plotting import make_movie_from_images
```

相关文件：

```text
src/solarutils/plotting/movie.py
workflows/observations/make_processed_aia_movies.py
tests/test_plotting_movie.py
```

运行命令：

```bash
PYTHONPATH=src python workflows/observations/make_processed_aia_movies.py \
  --variant multiscale_sharpen \
  --backend imageio \
  --fps 3 \
  --overwrite
```

已生成 7 个多尺度锐化 MP4：

```text
data/test/movies/observations/20241008_2300_2400/sdo/aia/
```

当前机器没有系统 `ffmpeg`，使用 `imageio` 后端。`imageio` 后端会自动把奇数宽/高 padding 到偶数，避免 H.264 编码错误。

### 9. HMI 45s LOS submap derotate workflow

已新增 HMI 45s LOS 磁图 workflow：

```text
workflows/observations/process_hmi45_derotate_submaps.py
workflows/observations/plot_processed_hmi45_submaps.py
```

处理规则：

- 第一张 HMI 45s 磁图 `rotate()` 后作为参考。
- 从参考图裁剪 submap，默认 `x=100-800 arcsec, y=150-500 arcsec`。
- 后续每张 HMI 45s 磁图 `rotate()` 后 derotate 到参考 submap。
- 支持 `--method center|pixel`。`center` 只旋转参考中心点；`pixel` 会将参考 submap 每个像素旋转到目标时刻，再在目标图像上插值取值。
- 默认推荐 `pixel`，更适合大范围 submap。
- 输出 FITS 使用参考 WCS/几何，时间关键字使用当前帧自己的时间。
- HMI 不做曝光时间缩放，不写 `EXPTIME`。
- 输出单位 `BUNIT='Gauss'`。

命令：

```bash
PYTHONPATH=src python workflows/observations/process_hmi45_derotate_submaps.py --method pixel --overwrite
PYTHONPATH=src python workflows/observations/plot_processed_hmi45_submaps.py --method pixel --overwrite
```

输出：

```text
data/test/processed/observations/20241008_2300_2400/sdo/hmi_45s/los_derotate_submap_pixel/
data/test/figures/observations/20241008_2300_2400/sdo/hmi_45s/los_derotate_submap_pixel/
```

验证结果：

```text
count: 7 FITS
shape: (695, 1390)
blank_count: 0
exptime_count: 0
bunit: Gauss
processed FITS size: 52M
single_png: 7
sequence_png: 1
```

center-vs-pixel 对比：

```text
tests/visual/compare_hmi45_derotation_methods.py
data/test/figures/observations/20241008_2300_2400/sdo/hmi_45s/compare_derotation_center_pixel.png
mean_abs_diff=12.9514 G
p95_abs_diff=35.948 G
p99_abs_diff=128.035 G
max_abs_diff=2127.91 G
```

### 10. 统一 derotate 方法入口和速度测试

`derotate_to_reference` 现在统一支持两种方法：

```python
from solarutils.obs.maps import derotate_to_reference

aligned = derotate_to_reference(path, reference_map)
aligned = derotate_to_reference(path, reference_map, method="center")
```

含义：

- `center`：只对参考 submap 中心点做太阳自转，然后裁剪和平移，速度快。
- `pixel`：默认方法，把参考 submap 的每个像素都旋转到目标时刻，再在目标图像上插值，更适合大范围区域。

AIA 批处理现在默认使用 pixel：

```bash
PYTHONPATH=src python workflows/observations/process_aia_prep_derotate_submaps.py \
  --overwrite
```

HMI 批处理已经改为通过统一入口调用：

```bash
PYTHONPATH=src python workflows/observations/process_hmi45_derotate_submaps.py \
  --method pixel \
  --overwrite
```

AIA 绘图工作流也默认读取 pixel 结果：

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --overwrite
```

速度测试脚本：

```text
tests/visual/benchmark_derotation_methods.py
```

测试结果保存到：

```text
data/test/processed/benchmarks/derotation_methods.txt
```

当前测试区域 `x=100-800 arcsec, y=150-500 arcsec`，`repeats=2`：

```text
hmi45 center best=4.645640 s  mean=4.773800 s
hmi45 pixel  best=4.857958 s  mean=5.055192 s
aia193 center best=0.346969 s  mean=0.364392 s
aia193 pixel  best=0.635687 s  mean=0.673612 s
```

## 已读取或参考的外部旧项目文件

主要参考了：

```text
/Users/xiaoxingyun/Desktop/20241008/cme/193.ipynb
```

里面有旧的 AIA 193 处理流程：

- `get_pointing_table("JSOC", ...)`
- `update_pointing`
- `aiapy.psf.calculate_psf`
- `aiapy.psf.deconvolve`
- `register`
- `de_rotate`
- `sharpen`

旧 notebook 中常用 submap 范围：

```text
[100, 800, 150, 500]
```

## 主要修改/新增文件

### 基础代码

```text
src/solarutils/obs/maps.py
src/solarutils/obs/__init__.py
src/solarutils/plotting/maps.py
src/solarutils/plotting/__init__.py
```

### 测试

```text
tests/test_obs_maps.py
tests/test_download.py
tests/visual/plot_aia_context_examples.py
tests/visual/preprocess_aia_examples.py
tests/visual/cache_aia_psf_examples.py
tests/visual/compare_aia193_notebook_flow.py
tests/visual/compare_aia193_local_psf_padding.py
tests/visual/verify_aia_sharpen.py
tests/visual/verify_aia_derotate.py
workflows/observations/process_aia_prep_derotate_submaps.py
workflows/observations/plot_processed_aia_submaps.py
workflows/observations/make_processed_aia_movies.py
tests/test_plotting_movie.py
```

### 文档

```text
workflows/observations/README.md
workflows/observations/download_zh.md
workflows/observations/aia_preprocessing.md
workflows/observations/aia_preprocessing_zh.md
```

注意：`aia_preprocessing.md` 是英文说明，`aia_preprocessing_zh.md` 是中文说明。`workflows/observations/README.md` 已提供中英文入口。

### 已删除

```text
workflows/observations/plot_aia171_context.py
```

原因：它属于 visual test，不属于正式 workflow。

## 当前约束和组织规则

### 目录规范

测试原始数据：

```text
data/test/observations/
```

测试图：

```text
data/test/figures/observations/<source>/<instrument>/<wavelength_or_series>/
```

校准测试数据：

```text
data/test/calibration/
```

视觉测试脚本：

```text
tests/visual/
```

正式 workflow：

```text
workflows/
```

长期维护基础函数：

```text
src/solarutils/
```

### 绘图规范

- 坐标标签统一：

```text
Solar-X (arcsec)
Solar-Y (arcsec)
```

- AIA 标题格式：

```text
AIA 193 Å 23:00:04 UT
```

- 子母图：
  - cyan 虚线框
  - cyan 虚线连接
  - 不要写成“字母图”

### 数据和 Git

`.gitignore` 当前会忽略：

- FITS/FTS
- NPZ/NPY
- PNG/PDF/EPS
- CDF
- MP4 等大文件

所以 PSF 缓存、测试 FITS、测试图默认不提交 Git。

## 已验证状态

最近一次测试结果：

```text
21 passed
```

## 尚未完成 / 下一步建议

1. 后续可以继续 AIA 时间序列处理：

- AIA register / optional PSF
- derotate 到参考 submap
- 可选 sharpen_map
- 生成序列图或 movie

2. 后续可以开始 HMI 数据处理 workflow。

3. 注意：当前 open tabs 中仍显示 `workflows/observations/plot_aia171_context.py`，但文件已经删除；这是 IDE 未刷新导致。

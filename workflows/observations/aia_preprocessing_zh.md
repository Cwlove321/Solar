# AIA 预处理说明

English version: [aia_preprocessing.md](aia_preprocessing.md).

这份文档说明本代码库中维护的 SDO/AIA level-1 到 level-1.5 预处理流程。核心目标是把散落在 notebook 里的常用流程整理成可复用、可测试、可长期维护的工具。

## 1. 目的

AIA level-1 FITS 可以用于快速查看图像，但正式分析通常需要：

- 修正 pointing metadata。
- 统一到 AIA level-1.5 的标准 plate scale 和方向。
- 根据需要做 PSF 去卷积。
- 在时间序列中把不同时间的图像裁剪、对齐到同一个参考区域。

当前维护的核心函数是：

```python
from solarutils.obs.maps import aia_prep
```

处理顺序为：

1. `update_pointing`：修正 AIA 指向信息。
2. 可选 `respike`：如果提供 spike 位置和值，则恢复 despiked pixels。
3. 可选 `aiapy.psf.deconvolve`：做 AIA PSF 去卷积。
4. `register`：旋转并重采样到标准 AIA level-1.5 网格。

重要约定：PSF 去卷积默认关闭。只有显式设置 `deconvolve_psf=True` 时才会执行。

## 2. 环境安装

在代码库根目录运行：

```bash
pip install -e ".[aia,test]"
```

如果还需要下载 Solar Orbiter 数据，可以安装完整观测依赖：

```bash
pip install -e ".[aia,solar-orbiter,test]"
```

## 3. PSF 缓存

AIA PSF 是全日面的 `4096 x 4096` 数组，重复计算比较慢，所以测试流程已经把常用 EUV 波段缓存到本地：

```text
data/calibration/sdo/aia/psf/
├── aia_0094_psf.npz
├── aia_0131_psf.npz
├── aia_0171_psf.npz
├── aia_0193_psf.npz
├── aia_0211_psf.npz
├── aia_0304_psf.npz
└── aia_0335_psf.npz
```

重新生成或补充缓存：

```bash
PYTHONPATH=src python tests/visual/cache_aia_psf_examples.py
```

在代码中使用缓存：

```python
from solarutils.obs.maps import aia_psf_path, cache_aia_psfs, load_aia_psf

cache_aia_psfs([171, 193])
psf_path = aia_psf_path(193)
psf = load_aia_psf(193)
```

这些 PSF 文件属于本地校准/测试数据，不建议提交到 Git。

## 4. 基本用法

只做 pointing 修正和 registration，不做 PSF 去卷积：

```python
import astropy.units as u
from sunpy.map import Map
from aiapy.calibrate.utils import get_pointing_table

from solarutils.obs.maps import aia_prep

aia_map = Map("aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits")
pointing_table = get_pointing_table(
    "JSOC",
    time_range=(aia_map.date - 12 * u.h, aia_map.date + 12 * u.h),
)

prepared = aia_prep(aia_map, pointing_table=pointing_table)
```

显式打开 PSF 去卷积：

```python
from solarutils.obs.maps import aia_prep, aia_psf_path

prepared = aia_prep(
    aia_map,
    pointing_table=pointing_table,
    deconvolve_psf=True,
    psf=aia_psf_path(193),
    deconvolve_iterations=25,
)
```

做快速 smoke test 时，可以跳过 pointing 和 PSF，只测试 registration 是否可运行：

```python
prepared = aia_prep(
    aia_map,
    correct_pointing=False,
    deconvolve_psf=False,
    register_map=True,
)
```

## 5. 批量处理建议

处理同一时间段的大量 AIA 图像时：

- `pointing_table` 只查询一次，然后复用。
- PSF 使用本地缓存路径，不要每张图重新计算。
- 正式处理建议 `deconvolve_iterations=25`。
- 快速看效果可以临时用 `deconvolve_iterations=5`。

示例：

```python
from pathlib import Path

import astropy.units as u
from sunpy.map import Map
from aiapy.calibrate.utils import get_pointing_table

from solarutils.obs.maps import aia_prep, aia_psf_path

files = sorted(Path("cases/20241008/data/2300_2400/sdo/aia").glob("*193*.fits"))
first_map = Map(files[0])
pointing_table = get_pointing_table(
    "JSOC",
    time_range=(first_map.date - 12 * u.h, first_map.date + 12 * u.h),
)
psf = aia_psf_path(193)

prepared_maps = [
    aia_prep(
        Map(path),
        pointing_table=pointing_table,
        deconvolve_psf=True,
        psf=psf,
    )
    for path in files
]
```

## 6. 局部 PSF 去卷积

局部 PSF 去卷积比全日面去卷积快很多，但它是近似方法。差别主要来自边界：全日面去卷积会考虑目标区域外的大尺度散射贡献，而直接对子图去卷积会截断区域外的信息，所以靠近子图边缘时误差更明显。

当前推荐函数是：

```python
from solarutils.obs.maps import deconvolve_submap_with_padding
```

它会先裁一个比目标区域更大的 padded submap，在这个大区域上做局部去卷积，然后只保留中间目标区域。默认设置：

```python
padding = 200.0  # arcsec
psf = None       # 可选；不传则按 aia_map.wavelength 自动计算
```

推荐用法：

```python
from sunpy.map import Map

from solarutils.obs.maps import deconvolve_submap_with_padding, aia_psf_path

aia_map = Map("aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits")
local = deconvolve_submap_with_padding(
    aia_map,
    extent=(100, 800, 150, 500),
    psf=aia_psf_path(193),
    padding=200,
    iterations=25,
)
```

如果不传 `psf`，函数会根据 `aia_map.wavelength` 自动计算 PSF：

```python
local = deconvolve_submap_with_padding(
    aia_map,
    extent=(100, 800, 150, 500),
    padding=200,
    iterations=25,
)
```

注意：

- `extent` 和 `padding` 的单位都是 arcsec。
- 子图尺寸需要是偶数，因为 `aiapy.psf.deconvolve` 内部使用 `rfft2/irfft2`。
- 函数内部会用 `ensure_even_shape` 自动裁成偶数尺寸。
- 直接拿 `4096 x 4096` 全日面 PSF 去卷积小子图会因为 FFT shape 不匹配而失败，所以局部处理会用 `crop_psf_to_shape` 把 PSF 中心裁剪并重新归一化。

## 7. 局部去卷积和全日面去卷积的差别

对 AIA 193 的目标区域：

```text
x = 100 到 800 arcsec
y = 150 到 500 arcsec
```

比较“局部 padded 去卷积”与“全日面去卷积后再裁剪”的结果：

```text
padding=0      mean |diff|=3.703   max=127.384
padding=100    mean |diff|=0.757   max=10.851
padding=200    mean |diff|=0.436   max=6.569
padding=400    mean |diff|=0.193   max=1.722
```

结论：

- `padding=200 arcsec` 是当前速度和误差的较好折中，因此作为默认值。
- 快速预览、交互式探索、调试裁剪区域时，可以使用局部 padded 去卷积。
- 正式定量分析更建议全日面 PSF 去卷积后再裁剪，或者使用更大的 padding 并检查误差。

测试脚本：

```text
tests/visual/compare_aia193_local_psf_padding.py
```

测试图：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/local_psf_padding_compare.png
```

## 8. 时间序列 derotation

旧 case notebook 中常见的 `de_rotate` 思路，是把每个时间的图像裁剪到参考图像的太阳自转后位置。当前维护函数为：

```python
from solarutils.obs.maps import derotate_to_reference, aia_psf_path

aligned = derotate_to_reference(
    path,
    reference_map,
    aia_prep=True,
    deconvolve_psf=True,
    psf=aia_psf_path(193),
    deconvolve_iterations=25,
    method="pixel",
)
```

如果只想做 AIA registration，不做 PSF 去卷积：

```python
aligned = derotate_to_reference(
    path,
    reference_map,
    aia_prep=True,
    deconvolve_psf=False,
)
```

`method="pixel"` 现在是默认方法。它会把参考 submap 的每个像素都旋转到当前输入图的时刻，再从当前图插值取值，更适合大范围 submap 和较长时间间隔；当前速度测试显示额外耗时不大。`method="center"` 会根据参考图像中心的太阳自转位置裁剪，并做亚像素平移，仍可用于快速预览：

```python
aligned = derotate_to_reference(
    path,
    reference_map,
    aia_prep=True,
    deconvolve_psf=True,
    psf=aia_psf_path(193),
    method="center",
)
```

这个函数会先按需执行 AIA 预处理，然后 derotate 到参考图，使输出图像尺寸与参考图一致。

输出 FITS metadata 规则：

- WCS 和几何信息使用参考图的 meta。
- `DATE-OBS`、`T_OBS` 等时间关键字改成当前输入图自己的时间。
- 数据会先除以当前输入图自己的曝光时间，再乘以参考图曝光时间。
- 因此输出数据是统一到参考曝光时间后的 `DN`，不是 `DN / s`。
- `EXPTIME` 保持和参考图一致，使同一组 derotate 后的 FITS 拥有统一参考曝光关键字。
- `BUNIT` 写为 `DN`。
- 太阳自转计算显式使用当前图的 `observer_coordinate`，因此不会再触发 SunPy 的 Earth-observer warning。

FITS 验证脚本：

```bash
PYTHONPATH=src python tests/visual/verify_aia_derotate.py --overwrite
```

验证输出：

```text
cases/20241008/products/processed/observations/2300_2400/sdo/aia/193/derotate_20241008T231006_to_230006.fits
```

本次 case 验证结果：

```text
reference shape: (584, 1167)
derotated shape: (584, 1167)
reference date-obs: 2024-10-08T23:00:04.843
source date-obs: 2024-10-08T23:10:04.843
output date-obs: 2024-10-08T23:10:04.843
output t_obs: 2024-10-08T23:10:05.843
output exptime: 2.000075
output bunit: DN
```

批量处理测试 AIA 数据：

```bash
PYTHONPATH=src python workflows/observations/process_aia_prep_derotate_submaps.py --overwrite
```

可以通过 `--method center|pixel` 选择 derotation 方法。现在默认使用 `pixel`，输出到独立目录；如果需要旧的快速近似结果，可以显式使用 `center`：

```bash
PYTHONPATH=src python workflows/observations/process_aia_prep_derotate_submaps.py \
  --method center \
  --overwrite
```

默认处理：

- 输入：`cases/20241008/data/2300_2400/sdo/aia/`
- 输出：`cases/20241008/products/processed/observations/2300_2400/sdo/aia/`
- 波段：94, 131, 171, 193, 211, 304, 335 Å。
- 每个波段第一张作为参考。
- 每张图先执行 `aia_prep(deconvolve_psf=True)`，再裁剪/derotate 到参考 submap。
- 默认 submap 范围：`x=100-800 arcsec, y=150-500 arcsec`。
- 默认 PSF deconvolution 迭代次数：`25`。

当前测试数据库已经生成：

```text
cases/20241008/products/processed/observations/2300_2400/sdo/aia/
├── 0094/prep_psf_derotate_submap/
├── 0131/prep_psf_derotate_submap/
├── 0171/prep_psf_derotate_submap/
├── 0193/prep_psf_derotate_submap/
├── 0211/prep_psf_derotate_submap/
├── 0304/prep_psf_derotate_submap/
└── 0335/prep_psf_derotate_submap/
```

如果使用 `--method pixel`，输出目录为：

```text
cases/20241008/products/processed/observations/2300_2400/sdo/aia/<wavelength>/prep_psf_derotate_submap_pixel/
```

验证结果：

```text
count: 49 FITS
shape: (584, 1168)
blank_count: 0
missing_bunit_count: 0
size: 238M
```

绘制并保存 processed AIA submap 图片：

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py --overwrite
```

绘图工作流也默认读取 pixel derotation 的结果。如果要绘制 center 结果：

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --derotate-method center \
  --overwrite
```

默认显示范围使用每个波段第一张图的 `1` 和 `99.5` 百分位，并把同一波段所有帧使用同一个 `vmin/vmax`。可以修改为：

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --norm-source first \
  --percentiles 1 99.5 \
  --overwrite
```

如果想用同一波段全部时间帧共同计算百分位：

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --norm-source all \
  --percentiles 1 99.5 \
  --overwrite
```

绘图时可以选择是否锐化，且 `vmin/vmax` 会使用锐化后的图像计算：

```bash
# 不锐化，默认
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --sharpen none \
  --percentiles 1 99.5 \
  --overwrite

# 单尺度 unsharp mask
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --sharpen single \
  --sigma 1 \
  --weight 1.5 \
  --percentiles 1 99.5 \
  --overwrite

# 多尺度 unsharp mask
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --sharpen multiscale \
  --sigmas 1 2 4 \
  --weights 0.8 0.35 0.15 \
  --percentiles 1 99.5 \
  --overwrite
```

锐化版本会保存到独立目录，避免覆盖未锐化版本：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/processed/single_sharpen/
cases/20241008/products/figures/observations/2300_2400/sdo/aia/processed/multiscale_sharpen/
```

输出：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/processed/
```

当前已生成：

```text
single_png: 49
sequence_png: 7
size: 44M
```

本次重画使用的默认设置：

```text
norm_source: first
percentiles: 1, 99.5
```

当前已额外生成多尺度锐化 quicklook：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/processed/multiscale_sharpen/
single_png: 49
sequence_png: 7
```

制作 AIA quicklook 电影：

```bash
PYTHONPATH=src python workflows/observations/make_processed_aia_movies.py \
  --variant multiscale_sharpen \
  --backend imageio \
  --fps 3 \
  --quality simple \
  --overwrite
```

如果需要更清晰的视频，用 `--quality high`；默认的 `simple` 适合快速预览和小文件。

基础函数：

```python
from solarutils.plotting import make_movie_from_images
```

说明：

- `backend=auto` 会优先使用系统 `ffmpeg`，找不到时使用 `imageio`。
- 当前机器没有系统 `ffmpeg`，测试使用 `imageio` 后端。
- MP4 写入前会自动把奇数宽/高 padding 到偶数，避免 H.264 编码器报错。
- `simple` 和 `high` 两个质量预设在 `ffmpeg` 和 `imageio` 后端都可用。

当前已生成 7 个多尺度锐化 MP4：

```text
cases/20241008/products/movies/observations/2300_2400/sdo/aia/
```

示例：

```text
cases/20241008/products/movies/observations/2300_2400/sdo/aia/0193/0193_prep_psf_derotate_submap_multiscale_sharpen.mp4
```

其中每个波段都有一张时间序列总览图，例如：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/processed/0193/0193_prep_psf_derotate_submap_sequence.png
```

## 9. 锐化

旧 case notebook 中还常用一个轻量的锐化步骤，一般放在 derotation 之后，用于增强图像细节。当前维护了两个函数：

```python
from solarutils.obs.maps import multiscale_sharpen_map, sharpen_map

sharpened = sharpen_map(aligned, sigma=1, weight=1.5)
multiscale = multiscale_sharpen_map(
    aligned,
    sigmas=(1, 2, 4),
    weights=(0.8, 0.35, 0.15),
)
```

单尺度计算公式为：

```text
sharpened = data + (data - gaussian_filter(data, sigma)) * weight
```

`multiscale_sharpen_map` 会在多个高斯尺度上提取细节，再按权重加回原图。当前默认参数比较温和：

```python
sigmas = (1.0, 2.0, 4.0)
weights = (0.8, 0.35, 0.15)
```

说明：

- 这两种方法都是 unsharp mask 类型的视觉增强。
- 它不能替代 PSF 去卷积。
- `sigma` 控制高斯模糊尺度。
- `weight` 控制锐化强度。
- 当前旧 notebook 默认参数是 `sigma=1, weight=1.5`。

验证脚本：

```bash
PYTHONPATH=src python tests/visual/verify_aia_sharpen.py
PYTHONPATH=src python tests/visual/verify_aia_sharpen.py --multiscale \
  --output cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/multiscale_sharpen_compare.png
```

验证图：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/sharpen_compare.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/multiscale_sharpen_compare.png
```

单尺度运行结果示例：

```text
original range: 76.000 to 2570.000
sharpened range: 49.621 to 3592.817
difference mean/std: -0.000 / 20.189
```

多尺度默认参数运行结果示例：

```text
sigmas=(1.0, 2.0, 4.0), weights=(0.8, 0.35, 0.15)
original range: 76.000 to 2570.000
sharpened range: 43.323 to 3700.406
difference mean/std: -0.000 / 20.931
```

## 10. HMI 45s LOS 磁图 submap derotate

HMI 45s LOS magnetogram 不需要 AIA 的 pointing、PSF、register 流程。当前 workflow 是：

1. 读取第一张 HMI 45s 磁图并执行 `rotate()`。
2. 从第一张图裁剪参考 submap。
3. 后续每张 HMI 45s 磁图执行 `rotate()`，再 derotate 到参考 submap。
4. 输出 FITS 使用参考图 WCS/几何信息，时间关键字使用当前帧自己的时间。
5. 保留磁场单位 `BUNIT='Gauss'`，不写 `EXPTIME`。

当前支持两种 derotate 方法：

- `center`：只旋转参考 submap 的中心点，然后做裁剪和平移，速度快，适合小区域快速预览。
- `pixel`：把参考 submap 的每个像素都旋转到目标时刻，再从目标图像插值取值，更适合大范围 submap 和较长时间间隔。

批处理命令：

```bash
PYTHONPATH=src python workflows/observations/process_hmi45_derotate_submaps.py --method pixel --overwrite
```

默认输入：

```text
cases/20241008/data/2300_2400/sdo/hmi_45s/
```

默认输出：

```text
cases/20241008/products/processed/observations/2300_2400/sdo/hmi_45s/los_derotate_submap_pixel/
```

默认 submap 范围：

```text
x = 100 到 800 arcsec
y = 150 到 500 arcsec
```

当前测试结果：

```text
count: 7 FITS
shape: (695, 1390)
blank_count: 0
exptime_count: 0
bunit: Gauss
size: 52M
```

AIA 和 HMI 现在共用同一个基础入口：

```python
from solarutils.obs.maps import derotate_to_reference

aligned = derotate_to_reference(path, reference_map)
aligned = derotate_to_reference(path, reference_map, method="center")
```

在当前测试区域 `x=100-800 arcsec, y=150-500 arcsec` 上，两种方法速度对比为：

```text
hmi45 center best=4.645640 s  mean=4.773800 s
hmi45 pixel  best=4.857958 s  mean=5.055192 s
aia193 center best=0.346969 s  mean=0.364392 s
aia193 pixel  best=0.635687 s  mean=0.673612 s
```

速度测试脚本和结果：

```text
tests/visual/benchmark_derotation_methods.py
data/benchmarks/derotation_methods.txt
```

绘制 HMI 45s quicklook：

```bash
PYTHONPATH=src python workflows/observations/plot_processed_hmi45_submaps.py --method pixel --overwrite
```

默认使用对称显示范围 `vmax=800 G`。pixel 方法输出图：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/hmi_45s/los_derotate_submap_pixel/
```

当前已生成：

```text
single_png: 7
sequence_png: 1
```

center 和 pixel 方法对比图：

```bash
PYTHONPATH=src python tests/visual/compare_hmi45_derotation_methods.py
```

最后一帧测试结果：

```text
mean_abs_diff=12.9514 G
p95_abs_diff=35.948 G
p99_abs_diff=128.035 G
max_abs_diff=2127.91 G
```

对比图：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/hmi_45s/compare_derotation_center_pixel.png
```

## 11. 与旧 notebook 的验证

已经和旧 notebook 对比过：

```text
/path/to/193.ipynb
```

旧流程核心为：

```text
update_pointing -> aiapy.psf.deconvolve -> register
```

使用同一张 AIA 193 FITS、同一个 JSOC pointing table、同一个 193 PSF、同样 5 次迭代时：

```text
max_abs_difference = 0
```

说明当前 `solarutils.obs.maps.aia_prep` 与旧 notebook 的核心预处理结果完全一致。

对比脚本：

```text
tests/visual/compare_aia193_notebook_flow.py
```

对比图：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/notebook_flow_compare.png
```

## 12. 相关文件

基础代码：

```text
src/solarutils/obs/maps.py
src/solarutils/obs/__init__.py
```

视觉测试：

```text
tests/visual/preprocess_aia_examples.py
tests/visual/cache_aia_psf_examples.py
tests/visual/compare_aia193_notebook_flow.py
tests/visual/compare_aia193_local_psf_padding.py
tests/visual/verify_aia_sharpen.py
tests/visual/verify_aia_derotate.py
```

HMI workflow：

```text
workflows/observations/process_hmi45_derotate_submaps.py
workflows/observations/plot_processed_hmi45_submaps.py
```

输出图片：

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/171/preprocess_register_smoke.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/171/preprocess_psf_register_submap.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/notebook_flow_compare.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/local_psf_deconvolution_compare.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/local_psf_padding_compare.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/sharpen_compare.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/multiscale_sharpen_compare.png
```

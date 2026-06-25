# prepare_rbsl_path.py 中文说明

`prepare_rbsl_path.py` 用于为球坐标 RBSL/flux-rope 模拟准备磁绳轴线文件 `path.dat`。脚本可以在 AIA 子图上交互式拾取路径点，也可以复用已有的 `pts.npy` 或直接传入 Stonyhurst 经纬度点，然后根据给定的顶点高度 `Hs` 生成三维轴线，并输出 AIA 投影图、经纬度半径文本和 AMRVAC 可读的二进制 `path.dat`。

## 基本功能

- 在 AIA 子图上交互取点，得到二维路径点 `pts.npy`。
- 将二维 Stonyhurst 路径重采样为固定点数的三维磁绳轴线，默认 `npoints=6000`。
- 使用 `--height-mm` 指定默认圆弧高度曲线的顶点高度 `Hs`，单位为 Mm。
- 可选叠加参考路径点，便于重新取点时对照旧路径。
- 可选生成半径为 `--rope-radius-mm`、扭转圈数为 `--rope-twist` 的边界螺旋线。
- 可选在 AIA 经纬度底图、HMI Br 底图和 3D Plotly 曲面上绘制路径。

## 常用运行环境

建议从仓库根目录运行：

```bash
cd /Users/xiaoxingyun/Desktop/solar_code

SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mplconfig PYTHONPATH=src \
python workflows/simulations/spherical/prepare_rbsl_path.py ...
```

其中：

- `PYTHONPATH=src` 用于找到本仓库的 `solarutils` 包。
- `SUNPY_CONFIGDIR` 和 `MPLCONFIGDIR` 可避免 SunPy/Matplotlib 配置写到系统默认目录。

## 交互取点

如果不提供 `--points` 或 `--path-points-stonyhurst`，脚本会打开 AIA 子图窗口用于取点。

操作方式：

- 鼠标左键按路径顺序点击控制点。
- `delete` 或 `backspace` 删除上一个点。
- `enter` 或 `return` 完成取点。
- 至少需要两个点，首尾点会作为磁绳两个足点。

如果设置：

```bash
--smooth-path-points 120 --smooth-path-smoothing 0
```

交互点击的原始控制点会在像素坐标中平滑/重采样为 120 个路径点，再转换到 Stonyhurst 经纬度。这个平滑只作用于交互点击得到的点；如果使用 `--points` 读取已有 `pts.npy`，脚本会直接使用该文件中的点。

## 关键参数

必需参数：

- `--aia`：用于交互取点和投影的 AIA FITS 文件。
- `--output-dir`：输出目录。
- `--submap-arcsec X0 X1 Y0 Y1`：AIA 取点/投影子图范围，单位 arcsec。
- `--height-mm Hs`：默认圆弧高度曲线的顶点高度，单位 Mm。

常用可选参数：

- `--base-map`：经纬度底图采样使用的 FITS；只有绘制经纬度/HMI/3D 图时需要。
- `--points`：读取已有 `pts.npy`，跳过交互取点。
- `--path-points-stonyhurst`：直接在命令行传入 `lon lat lon lat ...`。
- `--reference-points`：交互窗口中叠加旧的 `pts.npy` 作为二维参考线。
- `--npoints`：输出三维轴线点数，默认 6000。
- `--cmap`：底图色表，例如 `sdoaia193`。
- `--image-scale`：经纬度底图显示缩放，默认 `log`。
- `--image-clip-percent LOW HIGH`：经纬度底图显示百分位范围，默认 `1 99.5`。
- `--image-vmin-percent`、`--image-vmax-percent`：分别指定显示下限/上限百分位，会覆盖 `--image-clip-percent` 对应的值。
- `--rope-radius-mm`：可选边界螺旋线半径，单位 Mm。
- `--rope-twist`：可选边界螺旋线扭转圈数；需要和 `--rope-radius-mm` 一起使用。
- `--write-3d`：写出 Plotly 3D HTML；如果未提供经纬度范围，脚本会根据路径点自动估算。

经纬度图相关可选参数：

- `--lon-range LON0 LON1`：经纬度底图的经度范围，单位 degree。
- `--lat-range LAT0 LAT1`：经纬度底图的纬度范围，单位 degree。
- `--resolution-degree`：经纬度底图采样分辨率，默认 0.1 度。
- `--nghost`：经纬度网格 ghost cell 数，默认 2。
- `--image-power`、`--image-vmin`、`--image-vmax`：旧版幂律显示参数；一般不需要再给。只要显式给其中任意一个，脚本会使用旧的 power 模式。

如果不需要 `base_footpoints.png`、`picked_points_lonlat.png`、HMI Br 经纬度图或 `rope_3d.html`，可以省略 `--lon-range`、`--lat-range`、`--resolution-degree` 和 `--nghost`。如果只要求 `--write-3d` 且没有给 `--lon-range/--lat-range`，脚本会从路径点自动估算经纬度范围，默认在路径最小/最大经纬度外加 padding。

经纬度 AIA 底图默认使用 `log10(abs(data))`，并按 log 图像的 `[1, 99.5]` 百分位设置 `vmin/vmax`。因此通常不需要手动给 `--image-power 0.25 --image-vmin 2 --image-vmax 6`。

如果要调整百分位，可以任选一种写法：

```bash
--image-clip-percent 2 99
```

或：

```bash
--image-vmin-percent 2 --image-vmax-percent 99
```

HMI 相关可选参数：

- `--hmi-field`
- `--hmi-inclination`
- `--hmi-azimuth`
- `--hmi-disambig`

这些参数成组使用时，脚本会计算并绘制 HMI Br 经纬度底图。

HMI Br 图建议显式提供 `--lon-range` 和 `--lat-range`，因为磁图范围通常需要人工控制。若只给了其中一个经纬度范围，脚本会报错并要求二者同时提供。

## 主要输出文件

输出目录中常见文件如下：

- `pts.npy`：路径点，形状为 `(N, 2)`，列为 `(lon_deg, lat_deg)`。
- `clicked_points_raw.npy`：交互点击的原始像素点。
- `path.dat`：AMRVAC/RBSL 使用的二进制路径文件，写出时经度会加 180 度。
- `path_lonlat_radius.txt`：三维轴线文本文件，列为 `(lon_deg, lat_deg, radius_rsun)`。
- `boundary_lonlat_radius.txt`：可选边界螺旋线，只有同时设置 `--rope-radius-mm` 和 `--rope-twist` 时生成。
- `path_summary.json`：本次运行参数、足点、几何诊断和输出路径摘要。
- `picked_points.png`：AIA 子图上的拾取点或路径点图。
- `picked_points_lonlat.png`：可选，经纬度底图上的路径点图。
- `selected_footpoints.png`：AIA 子图上的两个足点。
- `base_footpoints.png`：可选，经纬度底图上的两个足点。
- `projected_path_aia.png`：三维路径投影回 AIA 子图后的图。
- `rope_3d.html`：可选 3D 可视化图，仅 `--write-3d` 时生成。

## 示例：交互重新取点，Hs = 50

下面示例使用 `193_2236` AIA 193 底图，参考旧路径 `hs40_smooth_repick5/pts.npy`，交互重新取点并生成 `Hs=50 Mm` 的默认高度路径。

```bash
SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mplconfig PYTHONPATH=src \
python workflows/simulations/spherical/prepare_rbsl_path.py \
  --aia cases/20241008/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --base-map cases/20241008/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --output-dir cases/20241008/products/rbsl_path_193_2236_pathonly_hs50_default_repick16 \
  --submap-arcsec 400 800 100 400 \
  --height-mm 50 \
  --reference-points cases/20241008/products/rbsl_path_193_2236_pathonly_hs40_smooth_repick5/pts.npy \
  --smooth-path-points 120 \
  --smooth-path-smoothing 0 \
  --cmap sdoaia193
```

## 示例：复用已有 pts.npy，跳过交互

如果已经有满意的 `pts.npy`，可以跳过取点，直接重新生成不同 `Hs` 的路径。

```bash
SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mplconfig PYTHONPATH=src \
python workflows/simulations/spherical/prepare_rbsl_path.py \
  --aia cases/20241008/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --base-map cases/20241008/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --output-dir cases/20241008/products/rbsl_path_193_2236_pathonly_hs50_from_existing_points \
  --submap-arcsec 400 800 100 400 \
  --height-mm 50 \
  --points cases/20241008/products/rbsl_path_193_2236_pathonly_hs50_default_repick16/pts.npy \
  --cmap sdoaia193
```

## 示例：生成半径边界线

如果需要在 `projected_path_aia.png` 中同时显示磁绳轴线和一条边界螺旋线，可以设置半径和 twist。例如半径 `25 Mm`、扭转 `5` 圈：

```bash
SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mplconfig PYTHONPATH=src \
python workflows/simulations/spherical/prepare_rbsl_path.py \
  --aia cases/20241008/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --base-map cases/20241008/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --output-dir cases/20241008/products/rbsl_path_193_2236_pathonly_hs50_r25_twist5 \
  --submap-arcsec 400 800 100 400 \
  --height-mm 50 \
  --points cases/20241008/products/rbsl_path_193_2236_pathonly_hs50_default_repick16/pts.npy \
  --rope-radius-mm 25 \
  --rope-twist 5 \
  --cmap sdoaia193
```

这会额外输出：

- `boundary_lonlat_radius.txt`
- 带边界线的 `projected_path_aia.png`

## 注意事项

- `--height-mm` 生成的是脚本内置的默认圆弧高度曲线，顶点高度为指定的 `Hs`。
- `path.dat` 中经度会自动加 180 度；读取时通常需要用对应工具减回 180 度。
- `--reference-points` 需要的是二维路径点 `pts.npy`，不是 `path_lonlat_radius.txt`。
- 如果交互取点后报错 `path points must have shape (N, 2) with at least two points`，通常是没有成功点到至少两个点，或直接关闭了窗口。重新运行并在图轴区域内点击后按 Enter 即可。
- `--rope-radius-mm` 只有在 `--rope-twist` 非零时才会生成边界线。
- 多视角 FSI/EUVI 投影不在本脚本中完成，通常使用同目录下的 `plot_multiview_path.py`。

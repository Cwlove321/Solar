# prepare_rbsl_path.py 中文说明

`prepare_rbsl_path.py` 用于为球坐标 RBSL/flux-rope 模拟准备磁绳轴线文件 `path.dat`。脚本现在只保留两类图像输入：`--aia` 和 `--hmi`。它可以在 AIA 或 HMI 子图上交互式拾取路径点，也可以复用已有的 `pts.npy` 或直接传入 Stonyhurst 经纬度点，然后根据给定的顶点高度 `Hs` 生成三维轴线，并输出投影图、经纬度半径文本和 AMRVAC 可读的二进制 `path.dat`。

## 基本功能

- 在 AIA 子图上交互取点，得到二维路径点 `pts.npy`。
- 将二维 Stonyhurst 路径重采样为固定点数的三维磁绳轴线，默认 `npoints=6000`。
- 使用 `--height-mm` 指定默认圆弧高度曲线的顶点高度 `Hs`，单位为 Mm。
- 可选叠加参考路径点，便于重新取点时对照旧路径。
- 可选生成半径为 `--rope-radius-mm`、扭转圈数为 `--rope-twist` 的边界螺旋线。
- 可选在 AIA 经纬度底图、HMI 经纬度底图和 3D Plotly 曲面上绘制路径。
- 如果提供单张 `--hmi` FITS，默认会计算两个足点在球面圆区域内的磁通量，并写入 `path_summary.json`。

## 常用运行环境

建议从复制后的 `Spherical_RBSL` 根目录运行：

```bash
cd /path/to/Spherical_RBSL

python Preprocess/prepare_rbsl_path.py ...
```

其中：

- 先使用 `python -m pip install -r Preprocess/requirements.txt` 安装依赖；脚本使用内置 `rbsl_core`，不需要 `PYTHONPATH`。
- `SUNPY_CONFIGDIR` 和 `MPLCONFIGDIR` 可避免 SunPy/Matplotlib 配置写到系统默认目录。

## 交互取点

如果不提供 `--points` 或 `--path-points-stonyhurst`，脚本会打开取点窗口。

- 提供 `--aia` 时，默认在 AIA 子图上取点。
- 如果不提供 `--aia` 但提供了 `--hmi`，则直接在 HMI 子图上取点。
- `--hmi` 可以是 HMI `Br`，也可以是 HMI45。只要是 HMI 图，脚本都会先自动 `rotate`，保证上方尽量对齐太阳北极。

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

- `--output-dir`：输出目录。
- `--submap-arcsec X0 X1 Y0 Y1`：AIA 取点/投影子图范围，单位 arcsec。
- `--height-mm Hs`：默认圆弧高度曲线的顶点高度，单位 Mm。

背景图相关说明：

- `--aia`：用于交互取点和投影的主背景 FITS，通常是 AIA。
- `--hmi`：可选单张 HMI FITS，用于 HMI 经纬度图、HMI 3D 曲面和足点通量计算；如果未提供 `--aia`，也可直接作为取点背景。
- `--hmi` 可以是 HMI `Br`，也可以是 HMI45 LOS 磁图；脚本会自动 `rotate` 后再显示或采样。

常用可选参数：

- `--base-map`：经纬度底图采样使用的 FITS；只有绘制经纬度/3D 图时需要。
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

如果不需要 HMI 经纬度图或 `rope_3d.html`，可以省略 `--lon-range`、`--lat-range`、`--resolution-degree` 和 `--nghost`。如果只要求 `--write-3d` 且没有给 `--lon-range/--lat-range`，脚本会从路径点自动估算经纬度范围，默认在路径最小/最大经纬度外加 padding。

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

- `--hmi`
- `--hmi-vmin`
- `--hmi-vmax`
- `--hmi-contour-levels -30 30`：同时提供 `--aia` 和 `--hmi` 时，在交互窗口、选点结果图和轴线投影图上叠加 HMI 的 ±30 G 等高线。正值用红色实线，负值用蓝色虚线。原始 HMI 通过 Howard 较差自转对齐到 AIA 的坐标参考时刻，再一次双线性采样到 AIA 网格；不预先旋转磁图，也不额外平滑。输出 `hmi_on_aia_grid.fits`，时间、插值方式写入 `path_summary.json` 的 `hmi_contour_alignment`。该对齐不重建磁场随时间的演化，也不会把 HMI45 LOS 转换为 Br。单独使用 HMI 作为底图时仍自动旋转。

提供 `--hmi` 时，脚本会在给定经纬度范围内绘制 HMI 经纬度底图；如果同时开启 `--write-3d`，也会输出 HMI 曲面上的 3D 磁绳图。

足点通量相关可选参数：

- `--footpoint-flux-radius-mm`：足点通量积分半径，单位 Mm。如果提供了 `--hmi` 但没有显式给这个参数，脚本默认使用 `20 Mm`。
- `--footpoint-flux-resolution-degree`：足点通量积分采样分辨率，默认 `0.02` 度。

只要给了 `--hmi`，脚本就会自动计算足点通量；默认半径为 `20 Mm`，也可以用 `--footpoint-flux-radius-mm` 覆盖。结果会写入 `path_summary.json` 的 `footpoint_fluxes` 字段，并在终端输出两个足点的 signed / unsigned flux。

如果 `--hmi` 是 HMI45 这类 LOS 磁图，而不是严格的 HMI `Br`，脚本会额外输出明确说明：该通量是基于 LOS 磁图值旋转/采样后的近似结果，不应当当作严格的径向 `Br` 通量解释。

HMI 图建议显式提供 `--lon-range` 和 `--lat-range`，因为磁图范围通常需要人工控制。若只给了其中一个经纬度范围，脚本会报错并要求二者同时提供。

## 主要输出文件

输出目录中常见文件如下：

- `pts.npy`：路径点，形状为 `(N, 2)`，列为 `(lon_deg, lat_deg)`。
- `path.dat`：AMRVAC/RBSL 使用的二进制路径文件，写出时经度会加 180 度。
- `path_lonlat_radius.txt`：三维轴线文本文件，列为 `(lon_deg, lat_deg, radius_rsun)`。
- `boundary_lonlat_radius.txt`：可选边界螺旋线，只有同时设置 `--rope-radius-mm` 和 `--rope-twist` 时生成。
- `path_summary.json`：本次运行参数、足点、几何诊断和输出路径摘要。
- `path_summary.json` 中的 `footpoint_fluxes`：如果提供了 `--hmi`，这里会记录两个足点的通量结果，包括 `signed_flux_Mx`、`positive_flux_Mx`、`negative_flux_Mx`、`unsigned_flux_Mx`、`mean_br_G` 等。
- `path_summary.json` 中的 `footpoint_flux_note`：如果输入的是 HMI45 / LOS 磁图，这里会记录近似说明。
- `picked_points.png`：AIA 子图上的拾取点或路径点图。
- `projected_path_aia.png`：三维路径投影回 AIA 子图后的图。
- `rope_3d.html`：可选 3D 可视化图，仅 `--write-3d` 时生成。
- `hmi_footpoints.png`、`hmi_picked_points.png`：如果提供了 `--hmi` 且给定经纬度范围，会输出 HMI 经纬度图。
- `rope_3d_hmi.html`：如果提供了 `--hmi` 且开启 `--write-3d`，会输出 HMI 曲面上的 3D 图。

## 示例：交互重新取点，Hs = 50

下面示例使用 `193_2236` AIA 193 底图，参考旧路径 `hs40_smooth_repick5/pts.npy`，交互重新取点并生成 `Hs=50 Mm` 的默认高度路径。

```bash
python Preprocess/prepare_rbsl_path.py \
  --aia cases/20241008/data/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --base-map cases/20241008/data/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
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
python Preprocess/prepare_rbsl_path.py \
  --aia cases/20241008/data/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --base-map cases/20241008/data/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --output-dir cases/20241008/products/rbsl_path_193_2236_pathonly_hs50_from_existing_points \
  --submap-arcsec 400 800 100 400 \
  --height-mm 50 \
  --points cases/20241008/products/rbsl_path_193_2236_pathonly_hs50_default_repick16/pts.npy \
  --cmap sdoaia193
```

## 示例：生成半径边界线

如果需要在 `projected_path_aia.png` 中同时显示磁绳轴线和一条边界螺旋线，可以设置半径和 twist。例如半径 `25 Mm`、扭转 `5` 圈：

```bash
python Preprocess/prepare_rbsl_path.py \
  --aia cases/20241008/data/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
  --base-map cases/20241008/data/downloads/aia_193_2236/aia.lev1_euv_12s.2024-10-08T223606Z.193.image_lev1.fits \
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

## 示例：计算两个足点在 r = 20 Mm 内的通量

下面示例复用已有 `pts.npy`，不重新交互取点，同时计算两个足点在 `r = 20 Mm` 球面圆区域内的 HMI 通量：

```bash
cd /path/to/Spherical_RBSL

python Preprocess/prepare_rbsl_path.py \
  --output-dir /data/cases/20241007/products/rbsl_path_193_2200_repick7 \
  --submap-arcsec 100 700 50 350 \
  --lon-range -17 47 \
  --lat-range 0 32 \
  --resolution-degree 0.1 \
  --nghost 2 \
  --height-mm 30 \
  --aia /cases/20241007/data/2200/sdo/aia/aia.lev1_euv_12s.2024-10-07T220006Z.193.image_lev1.fits \
  --hmi /data/cases/20241007/products/archive/diagnostics/rbsl_2200_lon-17_47_lat0_32_res02_test/hmi_br_submap.fits \
  --points /data/cases/20241007/products/rbsl_path_193_2200_repick7/pts.npy \
  --reference-points /data/cases/20241007/products/rbsl_path_193_2200_repick6/pts.npy \
  --cmap gray \
  --footpoint-flux-radius-mm 20
```

运行完成后，可在下面文件中查看结果：

- `/data/cases/20241007/products/rbsl_path_193_2200_repick7/path_summary.json`

其中 `footpoint_fluxes.results[0]` 和 `footpoint_fluxes.results[1]` 分别对应两个足点。

如果这里把 `--hmi` 换成 HMI45，例如：

```bash
--hmi /cases/20241007/data/2200/sdo/hmi_45s/hmi.m_45s.20241007_220045_TAI.2.magnetogram.fits
```

脚本仍然会给出通量，但会在终端和 `path_summary.json` 里额外写明：这是 LOS 近似通量。

## 注意事项

- `--height-mm` 生成的是脚本内置的默认圆弧高度曲线，顶点高度为指定的 `Hs`。
- `path.dat` 中经度会自动加 180 度；读取时通常需要用对应工具减回 180 度。
- `--reference-points` 需要的是二维路径点 `pts.npy`，不是 `path_lonlat_radius.txt`。
- 如果交互取点后报错 `path points must have shape (N, 2) with at least two points`，通常是没有成功点到至少两个点，或直接关闭了窗口。重新运行并在图轴区域内点击后按 Enter 即可。
- `--rope-radius-mm` 只有在 `--rope-twist` 非零时才会生成边界线。
- 多视角 FSI/EUVI 投影不在本脚本中完成，通常使用同目录下的 `plot_multiview_path.py`。

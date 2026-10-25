# Spherical RBSL 预处理

本目录包含球面 RBSL 工作流所需的全部 Python 入口和数值辅助代码。它随
`Spherical_RBSL_FullMHD` 整体分发，不依赖外部 `solarutils` 源码。

完整交付入口见根目录[QUICKSTART.md](../QUICKSTART.md)。完成边界预处理后，
可使用 `make_pdfi_brll.py` 将 `cutout.dat` 的无鬼点径向场
转换为 PDFI_SS 所需的 `Brll0000.fits`，并严格验证它与
`pfss_boundary.dat` 中对应 patch 一致：

```bash
python Preprocess/make_pdfi_brll.py /path/to/rbsl_boundary_product
```

多视角绘图默认显示整幅FSI/EUVI图，不再沿用旧事件的固定裁剪窗；需要局部视图时
显式传入`--fsi-submap-arcsec`或`--euvi-submap-arcsec`。

## 文件说明

| 文件或目录 | 作用 |
| --- | --- |
| `prepare_rbsl_boundary.py` | 从 HMI 矢量分段和 synoptic 图生成 PFSS 全球边界、局部三分量 `cutout.dat` 及参数摘要 |
| `prepare_rbsl_path.py` | 交互或复用取点结果，生成 6000 点三维 RBSL 轴线 `path.dat` |
| `plot_multiview_path.py` | 将轴线投影到 Solar Orbiter/FSI 与 STEREO/EUVI 图像 |
| `aia_to_spherical_vtk.py` | 将 AIA 图像采样到球面经纬网格并输出 VTS |
| `make_pdfi_brll.py` | 从边界产品生成无鬼点 PDFI_SS `Brll0000.fits` 并验证方向和数值 |
| `rbsl_core/` | 内置、经过本工作流使用的球面网格、磁场处理、SunPy map 和矢量绘图实现 |
| `requirements.txt` | 独立 Python 环境依赖 |

## 独立安装

进入复制后的 `Spherical_RBSL_FullMHD` 根目录：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r Preprocess/requirements.txt
```

依赖分工如下：

- `numpy`、`scipy`：数组、插值、稀疏矩阵和球面磁场预处理；
- `matplotlib`、`plotly`：二维快速图和可选三维交互图；
- `astropy`、`sunpy[map]`：FITS、坐标变换、WCS 和太阳图像；
- `opencv-python-headless`：将局部patch与全球synoptic图做模板匹配；图形取点由Matplotlib负责；
- `pyevtk`：输出 ParaView 可读的球面 VTS。

脚本支持直接运行，因此不需要 `PYTHONPATH`：

```bash
python Preprocess/prepare_rbsl_boundary.py --help
python Preprocess/prepare_rbsl_path.py --help
python Preprocess/plot_multiview_path.py --help
python Preprocess/aia_to_spherical_vtk.py --help
```

## 边界预处理输入

`prepare_rbsl_boundary.py` 需要：

1. 同一时刻的 HMI `field`、`inclination`、`azimuth`、`disambig`；
2. HMI `mrdailysynframe_polfil` 的 `Mr_polfil` 全球径向磁图；
3. 日面裁剪范围 `--extent-arcsec`；
4. 球面模拟 patch 的 `--lon-range`、`--lat-range` 和角分辨率。

程序执行以下步骤：

1. 应用 HMI 消歧 bit 并构造图像平面磁场分量；
2. 将图像平面分量转换为球面 `Br/Btheta/Bphi`；
3. 在指定 Stonyhurst 经纬网格上采样局部活动区；
4. 默认执行球面力、力矩、观测偏离和平滑联合优化；
5. 对 synoptic 径向场做保守重采样和局部 patch 替换；
6. 写出 AMRVAC 边界文件、参数摘要和快速检查图。

通用示例：

```bash
python Preprocess/prepare_rbsl_boundary.py \
  --field /data/hmi.field.fits \
  --inclination /data/hmi.inclination.fits \
  --azimuth /data/hmi.azimuth.fits \
  --disambig /data/hmi.disambig.fits \
  --synoptic /data/hmi.Mr_polfil.fits \
  --output-dir /data/processed_boundary \
  --extent-arcsec 400 800 100 400 \
  --lon-range 13 61 --lat-range 8 32 \
  --resolution-degree 0.1 \
  --nghost 2 --refine-level 1 \
  --nr 400 --length-r 0.8
```

### 磁场预处理默认值

| 参数 | 默认值 | 含义 |
| --- | ---: | --- |
| `--preprocess-mu3` | `0.1` | 保持接近观测磁场的惩罚权重 |
| `--preprocess-mu4` | `1.0` | 空间平滑权重 |
| `--preprocess-max-iter` | `5000` | 最大迭代数 |
| `--preprocess-report-interval` | `200` | 终端诊断输出间隔 |
| `--preprocess-tol` | `1e-4` | 收敛判据 |

默认开启预处理；仅在科学目标明确要求原始磁图时使用 `--no-preprocess`。应比较
`parameters.json` 中优化前后的力、力矩和磁通诊断，而不是只看快速图。

### 鬼点与 refinement

有效鬼点宽度为：

```text
effective_nghost = nghost * 2**(refine_level - 1)
```

例如 `nghost=2, refine_level=3` 时有效宽度为 8。局地 `cutout.dat` 的形状包含
这部分边缘；AMRVAC 的物理域和 block 内部 ghost cells 是另一层概念，不能混用。

### 边界输出

| 输出 | 内容 |
| --- | --- |
| `pfss_boundary.dat` | 全球经纬网格径向磁场，供 `pfss/` 读取 |
| `cutout.dat` | 局地 `Br/Btheta/Bphi` 三分量二进制边界，供 `mf/` 和 MHD 阶段读取 |
| `hmi_br_submap.fits` | 转换后的局部径向场 FITS |
| `parameters.json` | 输入文件、网格、索引、预处理和输出元数据 |
| `parameters.txt` | 可复制到 `mod_usr.t`、`amrvac.par` 的建议参数 |
| `cut.png`、`projection.png`、`synoptic.png` | 区域、投影和全球边界快速检查图 |

`cutout.dat` 的布局为 `int32 nx, ny`，随后依次写入 `float64 Br`、`Btheta`、
`Bphi`。`path.dat` 则没有文件头，二者不可混用。

## RBSL 路径生成

基本交互示例：

```bash
python Preprocess/prepare_rbsl_path.py \
  --aia /data/aia_193.fits \
  --hmi /data/hmi_br.fits \
  --output-dir /data/rbsl_path \
  --submap-arcsec 400 800 100 400 \
  --lon-range 13 61 --lat-range 8 32 \
  --height-mm 50 \
  --smooth-path-points 120
```

交互窗口中依次点击路径点，按 Enter 完成。若已有满意的 `pts.npy`：

```bash
python Preprocess/prepare_rbsl_path.py \
  --aia /data/aia_193.fits \
  --output-dir /data/rbsl_path_hs60 \
  --submap-arcsec 400 800 100 400 \
  --height-mm 60 \
  --points /data/old_path/pts.npy
```

主要输出：

- `pts.npy`：二维 `(lon_deg, lat_deg)` 路径点；
- `path.dat`：AMRVAC 使用的 6000 点轴线；
- `path_lonlat_radius.txt`：便于检查的三维轴线文本；
- `path_summary.json`：参数、足点、长度、高度和可选通量诊断；
- `projected_path_aia.png`：轴线投影快速图；
- `rope_3d.html`：设置 `--write-3d` 时生成；
- `boundary_lonlat_radius.txt`：同时设置磁绳半径和 twist 时生成。

提供 `--hmi` 后，程序默认计算两个足点半径 20 Mm 范围内的 signed、positive、
negative 和 unsigned flux。若输入是 LOS 磁图，输出只应解释为近似通量。

更多交互取点、复用路径、磁绳半径/twist 和足点通量示例见
`prepare_rbsl_path_README.md`。

## 生成Full-MHD案例

在本工作流根目录确认`case.delivery.example.json`中的科学参数，随后运行：

```bash
python workflow.py new \
  --boundary /data/processed_boundary \
  --axis /data/rbsl_path/path.dat \
  --settings case.delivery.example.json \
  --output /data/runs/my_event
python workflow.py check /data/runs/my_event --config-only
```

生成器在案例副本中准备`Brll0000.fits`，核对其与全球磁图patch一致，自动推导角域/
角向根网格/点数，并向局部和全球PFSS及后续阶段同步配置。无需手工复制到旧的
`pfss/`或`mf/`目录。预处理产生的旧格式`parameters.txt`保留作诊断，新的Full-MHD
入口读取真实数组与`parameters.json`，不会照抄其中旧径向网格建议。

本工作流的默认顺序为`01_pfss（或01_pfss_global_ct）→02_relax→03_insert→04_fullmhd`。
选择无绳方案时由配置明确跳过03，并关闭PFSS扣除磁绳Br。运行环境与后续操作见
[QUICKSTART.md](../QUICKSTART.md)。

## 常见问题

| 现象 | 检查项 |
| --- | --- |
| `ModuleNotFoundError` | 是否在当前Python环境安装了根目录`requirements.txt` |
| 找不到 `rbsl_core` | 应直接运行 `Preprocess/*.py`，不要单独复制一个脚本 |
| OpenCV 导入失败 | 安装`requirements-preprocess.txt`中的opencv-python-headless；图形取点另需Matplotlib图形后端 |
| SunPy 坐标报错 | 检查 FITS WCS、观测时间、observer 和输入分段是否一致 |
| patch 超出全球网格 | 检查经纬范围、角分辨率及数组纬度方向 |
| 预处理不收敛 | 查看迭代诊断，检查 NaN、边界、磁通不平衡并谨慎调整 mu3/mu4 |
| `path.dat` 点数不对 | 新案例生成器从字节数推导点数，并与path_summary及显式设置交叉检查 |
| PFSS/MHD restart 失败 | 查看pipeline.state.json和阶段日志，检查目标时间、网格及DAT验收结果 |

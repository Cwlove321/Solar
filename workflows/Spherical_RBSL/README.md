# 球面 RBSL：PFSS、磁摩擦与 MHD 工作流

本目录是一套可以整体复制和分发的球面 RBSL MPI-AMRVAC 项目。它将 Python
观测数据预处理、RBSL 磁绳轴线构造、PFSS 背景场、磁摩擦松弛和后续 MHD
演化放在同一个目录中，不依赖原始 `solar_code/src`。

## 目录结构

```text
Spherical_RBSL/
├── Preprocess/             # 独立 Python 预处理程序与内置 rbsl_core
├── pfss/                   # PFSS 背景场和 RBSL 磁绳插入
├── mf/                     # 磁摩擦松弛，生成 NLFFF
├── mhd_0beta/              # 零 beta MHD 演化
├── mhd_internal_e/         # 含内能、重力和冷致密物质的 MHD 演化
└── README.md
```

四个数值阶段按以下依赖关系运行：

```text
HMI B_720s + HMI synoptic + AIA
                  │
                  ▼
             Preprocess/
      ┌───────────┼────────────┐
      ▼           ▼            ▼
pfss_boundary  cutout.dat    path.dat
      │           │            │
      └──────┬────┴──────┬─────┘
             ▼           │
           pfss/         │
             │ restart   │
             ▼           │
            mf/ ◀────────┘
             │ restart
       ┌─────┴──────────┐
       ▼                ▼
 mhd_0beta/      mhd_internal_e/
```

## 一、安装 Python 预处理环境

整个 `Spherical_RBSL` 复制到新机器后，在项目根目录执行：

```bash
python -m venv .venv
source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r Preprocess/requirements.txt
```

推荐 Python 3.10 或更高版本。Windows PowerShell 的激活命令是：

```powershell
.venv\Scripts\Activate.ps1
```

预处理所需的球面网格、磁场变换、磁场预处理、路径写入和绘图辅助函数已经放在
`Preprocess/rbsl_core/`，不需要安装 `solarutils`，也不需要设置 `PYTHONPATH`。
MPI-AMRVAC、MPI、Fortran 编译器和集群调度器只用于后面的数值模拟阶段，需要
在服务器上另外配置。

## 二、生成球面磁场边界

`Preprocess/prepare_rbsl_boundary.py` 读取同一时刻的 HMI `field`、
`inclination`、`azimuth`、`disambig` 以及一张 HMI synoptic 径向磁图，生成
局地区域矢量边界和替换活动区后的全球 PFSS 边界。

```bash
python Preprocess/prepare_rbsl_boundary.py \
  --field /data/hmi.B_720s.TIME.field.fits \
  --inclination /data/hmi.B_720s.TIME.inclination.fits \
  --azimuth /data/hmi.B_720s.TIME.azimuth.fits \
  --disambig /data/hmi.B_720s.TIME.disambig.fits \
  --synoptic /data/hmi.mrdailysynframe.TIME.Mr_polfil.fits \
  --output-dir /data/rbsl_boundary \
  --extent-arcsec 400 800 100 400 \
  --lon-range 13 61 --lat-range 8 32 \
  --resolution-degree 0.1 \
  --refine-level 1
```

该入口默认启用球面 Wiegelmann 风格预处理，默认参数为：

```text
mu3 = 0.1
mu4 = 1.0
max_iter = 5000
tol = 1e-4
```

明确需要保留未预处理观测场时使用 `--no-preprocess`。有效鬼点宽度为
`nghost × 2**(refine_level-1)`，必须与实际 AMR 层级一致。

主要输出及放置位置：

| 预处理输出 | AMRVAC 目标位置 | 用途 |
| --- | --- | --- |
| `pfss_boundary.dat` | `pfss/potential_boundary/pfss_boundary.dat` | 全球径向 PFSS 底边界 |
| `cutout.dat` | `mf/nlfff_boundary/cutout.dat` | 局地三分量磁场边界 |
| `parameters.txt` | 人工核对各阶段参数 | `mod_usr.t` 和 `amrvac.par` 建议值 |
| `parameters.json` | 保留在事件结果目录 | 输入、网格、预处理和输出元数据 |

复制命令示例：

```bash
cp /data/rbsl_boundary/pfss_boundary.dat pfss/potential_boundary/pfss_boundary.dat
cp /data/rbsl_boundary/cutout.dat mf/nlfff_boundary/cutout.dat
```

## 三、生成 RBSL 磁绳轴线

`Preprocess/prepare_rbsl_path.py` 可在 AIA/HMI 图像上交互取点，也可以复用已有
`pts.npy`。它将二维日面路径转换为含高度的三维球面轴线，并写出 AMRVAC 使用
的 6000 点 `path.dat`。

```bash
python Preprocess/prepare_rbsl_path.py \
  --aia /data/aia_193.fits \
  --hmi /data/hmi_br.fits \
  --output-dir /data/rbsl_path \
  --submap-arcsec 400 800 100 400 \
  --lon-range 13 61 --lat-range 8 32 \
  --height-mm 50 \
  --rope-radius-mm 25 --rope-twist 5 \
  --write-3d
```

如果已经保存取点结果，可用 `--points /path/to/pts.npy` 跳过交互。输出完成后：

```bash
cp /data/rbsl_path/path.dat pfss/rbsl/path.dat
```

`path.dat` 是无文件头的 Fortran/Python 二进制流，依次保存 6000 个点的
`float64 lon_deg, lat_deg, radius_rsun`；写出时经度加 180°，AMRVAC 模板按
相同约定读取。`mhd_internal_e` 也复用该轴线放置冷致密物质。

## 四、同步模拟参数

以 `Preprocess` 生成的 `parameters.txt` 和 `parameters.json` 为准，同步：

- `pfss/mod_usr.t` 中的 `lat0/lat1/lon0/lon1`；
- 四个阶段 `amrvac.par` 中的 `domain_nx*`、`xprobmin/max*`、block 尺寸和
  refinement；
- `restart_from_file`、输出编号、保存间隔和终止条件；
- `mod_usr.t` 中的磁场、长度、密度和温度归一化参数。

不要把 README 或历史 `amrvac.par` 中的示例数值直接用于新事件。部分模板仍
保留 namelist 结束符 `/` 后面的历史参考行；真正生效的值必须放在对应 namelist
内部。

## 五、编译和运行顺序

### 1. PFSS

`pfss/` 读取 `potential_boundary/pfss_boundary.dat` 和 `rbsl/path.dat`，生成
去除局部 RBSL 径向贡献后的边界及势场 restart。

```bash
cd pfss
export AMRVAC_DIR=/path/to/amrvac
$AMRVAC_DIR/setup.pl -d=3
make -j 8
bsub < job
```

确认 `pfss/data/lfffNNNN.dat` 正常生成后再进入下一阶段。

### 2. 磁摩擦 MF

`mf/` 从 PFSS DAT restart，读取 `nlfff_boundary/cutout.dat`，进行磁摩擦松弛：

```bash
cd ../mf
$AMRVAC_DIR/setup.pl -d=3
make -j 8
bsub < job
```

不要只根据文件名选择结果，应结合磁摩擦日志确认真正完成松弛的最后一个
`nlfffNNNN.dat`。

### 3. MHD

根据科学问题选择一个分支：

- `mhd_0beta/`：零 beta 动力学演化；
- `mhd_internal_e/`：包含内能、重力以及沿 `pfss/rbsl/path.dat` 加载的冷致密
  日珥物质。

两个分支都从选定的 MF 结果 restart，并统一读取
`../mf/nlfff_boundary/cutout.dat`。

## 六、AIA 球面 VTK 和多视角检查

将 AIA 图像采样到球面网格并输出 ParaView VTS：

```bash
python Preprocess/aia_to_spherical_vtk.py \
  --aia /data/aia_304.fits \
  --output /data/aia_304_surface.vts \
  --lon-range 13 61 --lat-range 8 32 \
  --resolution-degree 0.1 --radius 1.001 --quicklook
```

将同一条 RBSL 轴线投影到 Solar Orbiter/FSI 和 STEREO/EUVI：

```bash
python Preprocess/plot_multiview_path.py \
  --path-dat pfss/rbsl/path.dat \
  --fsi /data/solo_fsi.fits \
  --euvi /data/stereo_euvi.fts \
  --output-dir /data/multiview
```

## 运行前检查清单

1. 四个 HMI 矢量分段来自同一观测时刻，形状和 WCS 一致；
2. synoptic 图与局部 patch 的坐标、数组方向和磁场符号正确；
3. `parameters.json` 中的 patch 索引、网格数、经纬范围和 refinement 一致；
4. `pfss_boundary.dat`、`cutout.dat`、`path.dat` 已复制到模板期望的位置；
5. PFSS → MF → MHD 的 restart 文件真实存在且网格完全相同；
6. 修改 `mod_usr.t` 后已经重新执行 `setup.pl` 和 `make`；
7. 正式计算前检查快速图、磁通平衡、预处理前后诊断量和路径投影。

Python 各入口的全部参数可用 `python Preprocess/<script>.py --help` 查看。更详细
的预处理说明见 `Preprocess/README.md`，交互路径说明见
`Preprocess/prepare_rbsl_path_README.md`。

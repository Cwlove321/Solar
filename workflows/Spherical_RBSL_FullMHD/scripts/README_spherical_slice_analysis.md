# 球坐标 AMR 统一切片分析工具

正式入口：

```text
spherical_slice_analysis.py
```

它从 MPI-AMRVAC DAT 读取原始变量，通过 simesh 按真实 AMR/拉伸坐标采样到均匀切片或薄层，先保存可复用的 `NPZ`，然后从缓存绘图。必须把生成 DAT 时使用的 `amrvac.par` 一并传入；程序不会在缺少拉伸参数时静默假设径向均匀。

## 1. 主要切片模式

固定 `phi`、固定 `theta`、固定 `r` 是主要入口：

```bash
python spherical_slice_analysis.py extract r14_0001.dat amrvac.par phi_0001.npz \
  --mode phi --value 204.9573717 --nu 420 --nv 360

python spherical_slice_analysis.py extract r14_0001.dat amrvac.par theta_0001.npz \
  --mode theta --value 69.9256251 --nu 480 --nv 360

python spherical_slice_analysis.py extract r14_0001.dat amrvac.par r_0001.npz \
  --mode r --value 1.05 --nu 480 --nv 360
```

- `phi/theta --value` 使用度；
- `r --value` 使用太阳半径；
- 固定 `phi` 输出坐标为纬度–高度；
- 固定 `theta` 输出坐标为经度–高度；
- 固定 `r` 输出坐标为经度–纬度。

### 1.1 原生 AMR 网格切片

`plot-amr-grid` 不把 DAT 重采样到均匀网格，而是直接使用 simesh 根据 DAT 树和
`amrvac.par` 重建的真实拉伸节点，绘制与固定 `phi` 或固定 `theta` 平面相交的
叶块轮廓及块内单元线：

```bash
python spherical_slice_analysis.py plot-amr-grid \
  r15_0004.dat amrvac.par amr_phi_0004.png \
  --mode phi --value 204.95960924035637 \
  --min-height-mm 0 --max-height-mm 200 \
  --marker-u 20.07324443847174 --marker-v 39.971263344737324 \
  --marker-label "initial rope-axis midpoint"

python spherical_slice_analysis.py plot-amr-grid \
  r15_0004.dat amrvac.par amr_theta_0004.png \
  --mode theta --value 75 \
  --min-height-mm 0 --max-height-mm 200
```

- 固定 `phi` 图的横轴是纬度，固定 `theta` 图的横轴是经度；纵轴都是高度；
- 彩色轮廓和浅色填充表示原生 AMR 叶块，细线表示每个叶块内部的真实单元；
- `--u-min/--u-max` 可裁剪横轴，`--hide-cell-lines` 可只保留叶块；
- `--marker-u/--marker-v` 必须成对使用，可标记磁绳轴线或其它参考位置；
- 若算例时间单位不同于工作流默认值，用 `--unit-time-s` 传入日志中的精确时间单位；
- PNG 同目录会生成 JSON，记录 DAT/参数文件、坐标来源、时间、全域及切面各层叶块数、
  各层块宽和输出 SHA256；
- 该命令只读取 DAT，不修改模拟、AMR 树或物理参数。

## 2. 自由角度平面

自由平面采用球坐标中心、Cartesian 法向和向上方向：

```bash
python spherical_slice_analysis.py extract r14_0001.dat amrvac.par plane_0001.npz \
  --mode plane \
  --plane-origin-spherical 1.0575,20.07,204.96 \
  --plane-normal 0.2,-0.9,0.38 \
  --plane-up 0,0,1 \
  --plane-width-mm 200 --plane-height-mm 160 \
  --nu 480 --nv 400
```

也可用 `--plane-origin-cartesian=x,y,z`，单位为太阳半径。程序把均匀平面点转换成 `(r,theta,phi)` 后交给 simesh 采样，并把球坐标矢量转换到平面基底。

## 3. 三维导数后端

温度、密度、压力、速度、磁场、阿尔芬速度、等离子体 beta、辐射损失和背景加热只需单面。

`J=curl(B)`、`div(B)`、`div(v)`、`J/B` 和欧姆耗散有两种后端。

### 3.1 `local-fit`：已有 DAT 的推荐方法

程序先确定每个切片点所属的真实 AMR leaf 和单元，再按该单元的
`dr/dtheta/dphi` 在三个坐标方向做多尺度对称采样。默认偏移为当地单元宽度的
`1.5、2.0、2.5` 倍，通过局部最小二乘拟合一阶导数，然后代入完整球坐标旋度和
散度公式：

```bash
python spherical_slice_analysis.py extract r14_0001.dat amrvac.par phi_fit_0001.npz \
  --mode phi --value 204.9573717 \
  --slab-points 1 \
  --derivative-backend local-fit \
  --fit-scales 1.5,2.0,2.5 \
  --fit-residual-max 0.35 \
  --fit-sensitivity-max 0.35 \
  --nu 420 --nv 360
```

这种方法不需要预先构造均匀三维薄层，能够降低“插值后再差分”放大的 AMR 块面
条纹。缓存同时包含：

```text
owner_leaf_id
owner_amr_level
cross_leaf_stencil
cross_level_stencil
stencil_outside_domain
derivative_fit_residual
derivative_scale_sensitivity
derivative_reliable_mask
derivative_strict_reliable_mask
velocity_derivative_reliable_mask
velocity_derivative_strict_reliable_mask
```

用于定量图时优先选择：

```text
current_magnitude_reliable_code
j_over_b_reliable_per_Mm
ohmic_heating_reliable_cgs
divb_reliable_code
divv_reliable_code
```

默认的 `reliable` 字段屏蔽超出物理边界、拟合残差过大或对拟合尺度敏感的点，
但不会仅因模板跨越 AMR level 就屏蔽。实测跨层点与同层点的拟合误差分布相近，
若一刀切去掉跨层点，反而会在图中留下块状空洞。跨层位置仍由
`cross_level_stencil` 明确标记。

需要最保守的块边界排查时，可使用 `strict_reliable` 掩膜或对应字段；它们在默认
质量判据之外还会排除跨 AMR level 的点。默认 `reliable` 适合连续的科学图像，
`strict_reliable` 适合判断某个结构是否恰好压在层级接口上。

### 3.2 `uniform-slab`：快速预览和旧结果复现

旧方法要求奇数个薄层面：

```bash
python spherical_slice_analysis.py extract r14_0001.dat amrvac.par phi_slab_0001.npz \
  --mode phi --value 204.9573717 \
  --derivative-backend uniform-slab \
  --slab-points 5 --slab-half-width 0.2 \
  --nu 420 --nv 360
```

`--slab-half-width` 的单位：

- 固定 `phi/theta`：度；
- 固定 `r` 和自由平面：Mm。

固定球坐标切片使用完整球坐标旋度和散度公式；自由平面在正交 Cartesian 薄层
中计算三维导数。`uniform-slab` 便于快速预览，但块接口两侧的插值 stencil
不保证一阶导数连续，所以不能作为默认定量结果。

无论选择哪种后端，缓存中的 `divb_code` 都是基于单元中心场的后处理诊断，不等价
于 AMRVAC CT 面心磁场的 `get_divb()`。严格 CT 散度必须在 AMRVAC 原生网格上计算。

不需要导数量时使用：

```text
--derivative-backend none
```

`--derivative-backend auto` 是兼容模式：`slab-points>=3` 时选择
`uniform-slab`，否则不计算导数。

## 4. 绘图

```bash
python spherical_slice_analysis.py plot phi_0001.npz temperature_density.png \
  --fields temperature_K number_density_cm3 \
  --overlay magnetic \
  --limit temperature_K:3e3:4e6 \
  --limit number_density_cm3:2e8:1e14

python spherical_slice_analysis.py plot phi_fit_0001.npz diagnostics.png \
  --fields pressure_cgs speed_km_s alfven_speed_classical_km_s \
           plasma_beta radiative_loss_model_cgs background_heating_cgs \
           current_magnitude_reliable_code j_over_b_reliable_per_Mm \
           ohmic_heating_reliable_cgs divb_reliable_code divv_reliable_code \
  --overlay velocity
```

主要可用字段会写进每个缓存的 `available_fields`，同时记录 DAT、参数文件、simesh 坐标来源、采样统计、单位、辐射参数和导数口径。

## 5. 物理口径

- DAT 中 `e` 按本工作流约定是内能；`saveprim=.true.` 只影响 VTU，不改变 DAT 中的守恒变量。
- 半相对论模拟不会使用简单的 `m/rho`，程序按 AMRVAC 半相对论守恒动量反演速度。
- 同时给出经典阿尔芬速度和约化光速限制后的诊断速度，二者不混用。
- `radiative_loss_model_cgs` 重现 JCcorona 表、截止温度和径向阻尼。DAT 不保存动态 TRAC 截止温度及时间积分源项账本，因此它是当前状态上的模型估计，不宣称等于求解器本步实际积分值。
- 当前背景加热按源码使用 `Q0 exp(-h/60 Mm)`。
- 欧姆耗散按参数文件中的常数 `mhd_eta` 计算；若案例启用了自定义空间变电阻，应先扩展工具读取相同的 `eta(x)`。

## 6. 与旧脚本对照

固定 `phi` 缓存可以直接与已验证脚本逐点比较：

```bash
python spherical_slice_analysis.py validate-legacy \
  phi_0001.npz plot_fixed_phi_temperature_density_fieldlines.py \
  r14_0001.dat amrvac.par legacy_numeric_comparison.json

python spherical_slice_analysis.py compare-images \
  legacy.png unified.png legacy_unified_side_by_side.png
```

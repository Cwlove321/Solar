# 从背景 DAT 估算 RBSL 平衡电流

`calculate_rbsl_balance_from_dat.py` 用未插入磁绳的背景 AMRVAC DAT 和
RBSL 轴线，估算薄环近似下抵消外部束缚场所需的电流。它适用于**均匀网格**。
拉伸径向网格必须使用
`Spherical_RBSL_FullMHD/scripts/calculate_rbsl_balance_from_stretched_dat.py`。

## 输入

- `--dat`：松弛后、插绳前的背景 DAT。不要对已经含 RBSL 磁场的 DAT 计算，
  否则束缚场会混入磁绳自场。
- `--path-dat`：轴线，支持三列 ASCII `lon_deg lat_deg r/Rsun`，也支持成熟
  RBSL 工作流写出的无表头小端 `float64` 三元组。
- `--minor-radius-mm`：RBSL 小半径 `a`（Mm）。
- `--solar-radius-mm`：代码长度单位（默认 695.5 Mm）。
- `--unit-b-gauss`：可选，但需要输出等效轴向磁通 Mx 时必须提供。

脚本会从 DAT 的 `b1,b2,b3`（也接受 `br,btheta,bphi` 等别名）直接在每个轴线
点做线性插值。它输出 JSON 汇总、完整轴线 CSV 和三联图。

## 薄环力平衡近似的原理

把磁绳近似成大半径为 `R`、小半径为 `a` 的细圆环电流通道，且 `a << R`。在
局部轴线处，背景场沿磁绳截面法向的分量 `B_strap` 是把环向内拉回的束缚场；
磁绳自身的环形电流则产生向外的 hoop force（自膨胀力）。采用本 RBSL 代码的
无量纲约定，每单位轴长的两项主导力可写成

```math
f_{hoop} = \frac{I^2}{R/L_0}
\left[\ln\left(\frac{8R}{a}\right)-1\right],
\qquad
f_{strap} = I B_{strap}.
```

这里 `L0` 是长度代码单位，`I`、`B_strap` 均是代码单位。令外向 hoop force 与
内向束缚力相等，并取非零电流支，得到

```math
I_{balance} =
\frac{(R/L_0)B_{strap}}{\ln(8R/a)-1}.
```

因此，该诊断问的是：在当前背景场下，要让一根给定几何的细磁环不立刻整体向外
膨胀或向内收缩，RBSL 电流应处于什么量级。它不是由 DAT 直接积分 `J` 得到的
电流，而是由外部束缚场反推的**平衡所需电流**。

这个近似要求：

- `a/R` 足够小，且 `ln(8R/a)-1 > 0`；
- 背景场在一个磁绳截面尺度 `a` 上变化缓慢；
- 轴线接近单一圆弧，故一个 `R` 能代表其曲率；
- 压力梯度、重力、扭转磁压、线系端点和三维非局部 Lorentz 力相对次要。

日冕 RBSL 磁绳通常不完全满足这些条件。因此结果应被用来确定插入磁通/电流的
合理扫描范围，而不是宣称获得了严格 MHD 平衡。特别是 `B_strap(s)` 变号或沿轴线
变化强烈时，单一 `I` 不可能逐点平衡，必须在插绳后进行 CT 连续性和 MHD 松弛测试。

## 数值计算方法

1. 将轴线转成三维物理坐标，并在最佳拟合平面上做最小二乘圆拟合，得到圆心和
   大半径 `R`。若已知应采用的 `R`，用 `--major-radius-mm` 覆盖拟合值；JSON 中
   会同时保留拟合残差，便于判断薄环近似是否可信。
2. 由轴线切向 `t` 与指向拟合圆心的横向单位向量 `n` 定义
   `b = t x n`，并投影背景场得到

   ```math
   B_{strap}(s) = -\mathbf{B}_{bg}(s)\cdot\mathbf{b}(s).
   ```

3. 逐点使用上面的薄环力平衡近似

   ```math
   I_{balance}(s) =
   \frac{(R/L_0) B_{strap}(s)}{\ln(8R/a)-1},
   ```

   其中 `L0=solar_radius_mm`，`B` 和 `I` 均为代码单位。默认只对轴线弧长
   10%--90% 做加权平均，避开端点及其边界效应。
4. 给出 `--unit-b-gauss` 后，将电流等效为 RBSL 轴向磁通幅值：

   ```math
   F = |I| B_0 L_0^2 \frac{12\pi(a/L_0)}{5\sqrt{2}}.
   ```

   这个关系不是单纯由量纲分析得到的，而是由当前 RBSL 模型选定的磁绳截面
   结构决定的。RBSL 把磁绳场写成两部分矢势之和：

   ```math
   \mathbf{B}_{\mathrm{FR}}
   = \nabla\times\mathbf{A}_I+\nabla\times\mathbf{A}_F.
   ```

   `A_I` 与沿轴线的净电流成正比，产生磁绳的环向磁场；`A_F` 与轴向磁通成
   正比，产生磁绳的轴向磁场。二者原则上可以独立缩放，但任意比例通常会使
   直磁绳截面内存在径向 Lorentz 力。

   Titov et al. (2018) 的 RBSL 正则化核以直、圆柱形磁绳为局部基准。令归一化
   截面半径为 `rho=r_perp/a`，在 `0<=rho<=1` 内选取抛物线轴向电流密度

   ```math
   j_{\mathrm{ax}}(\rho)=2(1-\rho^2),
   ```

   其中 `j_ax` 以 `I_phys/(pi a^2)` 归一化；同时要求轴向磁场在 `rho=1` 处
   消失。将这个电流剖面代入直圆柱的无力条件

   ```math
   \mathbf{J}\times\mathbf{B}=0
   ```

   后，环向场与轴向场的相对幅值不再自由，而满足

   ```math
   F=\pm\frac{3}{5\sqrt{2}}\,\mu I_{\mathrm{phys}}a.
   ```

   正负号选择两种螺度。系数 `3/(5 sqrt(2))` 因而来自上述抛物线电流剖面、
   边界条件和对应的 RBSL 正则化核；换成另一种截面电流分布，这个数值系数也会
   改变。

   本项目代码没有直接使用有量纲电流 `I_phys`，而定义无量纲 RBSL 电流

   ```math
   I=\frac{\mu I_{\mathrm{phys}}}{4\pi B_0L_0},
   \qquad
   F_{\mathrm{code}}=\frac{F}{B_0L_0^2},
   \qquad
   a_{\mathrm{code}}=\frac{a}{L_0}.
   ```

   把这些归一化代入上一式，就得到代码实际采用的关系

   ```math
   F_{\mathrm{code}}
   =|I|\frac{12\pi a_{\mathrm{code}}}{5\sqrt{2}},
   ```

   恢复物理磁通单位后即为本节开头的公式。`12 pi` 来自物理电流与代码电流
   之间的 `4 pi` 归一化，再乘上直圆柱无力关系中的系数 `3`。

   其中，`I` 是 RBSL 代码单位下的无量纲轴向电流，`B_0` 是磁场单位，`L_0`
   是长度单位，`a` 是磁绳小半径。计算 Mx 时，`B_0` 应使用 G，`L_0` 和 `a`
   应使用相同的长度单位，并将 `L_0` 换算成 cm；因此 `B_0 L_0^2` 的单位就是
   `G cm^2 = Mx`，而 `a/L_0` 是无量纲量。

   这里使用 `|I|`，所以 `F` 表示非负的**轴向磁通幅值**。电流正负号由磁绳
   轴线方向和螺度约定决定，不包含在 `F` 的幅值中。反过来，给定待插入磁通
   幅值时，对应的电流幅值为

   ```math
   |I| = \frac{F}{B_0 L_0^2}
   \frac{5\sqrt{2}}{12\pi(a/L_0)}.
   ```

   这里的“等效”表示把薄环平衡诊断得到的 `I_balance` 按同一 RBSL 截面模型
   换算成可与 `F_flx` 直接比较的轴向磁通参数；它不是从背景 DAT 的某个截面
   直接积分 `B_axial` 得到的独立磁通测量。还要注意，以上关系只保证 RBSL
   **内部直圆柱基准**近似无力，并不保证弯曲磁绳与外部背景场、压力和重力达到
   **全局平衡**；后者才是本脚本用 `B_strap` 估算 `I_balance` 的目的。
   `--current-flux-mx` 会同时输出待插入磁通对应的带符号电流及其相对平衡电流
   的比值。

   公式来源：Titov et al. (2018), *Regularized Biot-Savart Laws for Modeling
   Magnetic Flux Ropes*, ApJL 852, L21, Equation (12),
   DOI: `10.3847/2041-8213/aaa3da`。

这是局部薄环诊断，不是完整 Grad-Shafranov 平衡。若 `B_strap` 在轴线上变号、
`I_balance(s)` 离散很大，或拟合圆残差接近 `a`，不存在一个单一电流能在整个
磁绳上实现平衡。

## 用法

```bash
cd /Users/xiaoxingyun/Desktop/solar_code
PY=/opt/anaconda3/bin/python

$PY workflows/Spherical_RBSL/scripts/calculate_rbsl_balance_from_dat.py \
  --dat /path/to/relax/data/relax_0006.dat \
  --path-dat /path/to/rbsl/path.dat \
  --minor-radius-mm 15 \
  --solar-radius-mm 695.5 \
  --unit-b-gauss 1.9977 \
  --current-flux-mx 1.3e20 \
  --helicity negative \
  --output-dir /path/to/relax/python/balance_current
```

输出目录包含：

- `rbsl_balance_summary.json`：全局平衡电流、等效磁通、圆拟合质量和 DAT 元数据；
- `rbsl_balance_profile.csv`：沿轴线的三分量磁场、束缚场、电流和可选磁通；
- `rbsl_balance_profile.png`：沿轴线的诊断图。

## 解释与检查

- 轴线方向会改变带符号的 `B_strap` 和 `I_balance`；比较正负螺度时必须使用与
  插绳相同的 `path.dat` 方向。
- `matching_helicity_fraction` 应在中央段接近 1。若不是，先检查轴线方向、
  边界磁场符号和局部磁场拓扑，而不是直接调整磁通。
- 每次显著改变待扣除的 `Br_RBSL` 或 PFSS/PDFI 背景后，都应重新松弛并再次计算。
- DAT 不保存拉伸规则；即便文件名相同，也不能用本脚本读取几何拉伸 DAT。

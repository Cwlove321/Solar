# 拉伸网格 DAT 的 RBSL 平衡电流

`calculate_rbsl_balance_from_stretched_dat.py` 是全热力学工作流的平衡电流
诊断入口。它的物理公式、输出文件和参数含义与
[薄环推导与磁通换算](balance_current_derivation.md)
相同；区别仅在于它从**真实拉伸网格 DAT**采样磁场。

## 薄环平衡的物理含义

该脚本求的不是 DAT 内部已存在的 `J`，而是给定背景束缚场时，细圆环 RBSL 磁绳
所需的平衡电流。将磁绳视为大半径 `R`、小半径 `a` 且 `a << R` 的电流环：环电流
的自场产生向外 hoop force，背景场在截面法向的 `B_strap` 产生向内 Lorentz 力。
在本代码归一化下，两者相等给出

```math
I_{balance} =
\frac{(R/L_0)B_{strap}}{\ln(8R/a)-1}.
```

脚本沿轴线计算局部 `B_strap(s)`，然后在中央 10%--90% 轴段做弧长加权平均。
它假设背景场在尺度 `a` 上平滑、轴线接近单一圆弧，并忽略压力、重力、线系端点、
扭转和完整三维非局部力。因而输出是插入磁通的物理量级约束，不保证插绳后的严格
MHD 平衡；`B_strap` 变号或局部平衡电流离散很大时尤应如此。

完整推导、磁通-电流换算和符号约定见
[薄环推导与磁通换算](balance_current_derivation.md)。

## 为什么必须提供 `amrvac.par`

MPI-AMRVAC DAT 本身只保存域范围、整数网格大小和 AMR 树，不保存
`stretch_dim`、`qstretch_baselevel` 或每个块的实际节点坐标。若把拉伸径向
网格当成 `linspace(rmin,rmax,N)`，轴线会被关联到错误的物理半径，进而污染
`B_strap`、平衡电流和等效磁通。

该脚本强制要求 `--grid-config`。它把原始 `amrvac.par` 传给新版 simesh，由
simesh 按每个 AMR 层的真实几何节点重建 `(r,theta,phi)`，再对轴线点进行物理
坐标线性插值。若该 par 文件没有启用拉伸，脚本会停止而不是静默按均匀网格运行。

## 用法

从工作流根目录运行。`python` 表示已安装分析依赖的解释器；安装步骤见
[依赖说明](../DEPENDENCIES.md)。请替换案例路径、背景帧号和物理参数。

```bash
python scripts/calculate_rbsl_balance_from_stretched_dat.py \
  --dat runs/my_case/02_relax/data/relax_0006.dat \
  --grid-config runs/my_case/02_relax/amrvac.par \
  --path-dat runs/my_case/03_insert/path.dat \
  --minor-radius-mm 15 \
  --solar-radius-mm 695.5 \
  --unit-b-gauss 1.9977 \
  --current-flux-mx 1.29e20 \
  --helicity negative \
  --output-dir runs/my_case/02_relax/python/balance_current
```

`--grid-config` 必须是产生该 DAT 的同一运行版本的 `amrvac.par`。特别是要检查：

- `stretch_dim(1)` 为实际的 `uni` 或 `symm` 设置；
- `qstretch_baselevel(1)` 与运行一致；
- 域范围、`domain_nx*`、块大小和 AMR 级数没有被后续模板改写。

建议首先在 `0000` 或已核对过 VTU 的一个快照上运行，并检查 JSON 中：

- `grid_mode` 是 `stretched`；
- `grid_config` 指向正确的 par 文件；
- `axis_circle_fit_rms_mm` 远小于磁绳小半径；
- `valid_central_points` 覆盖中央 10%--90% 轴线；
- `matching_helicity_fraction` 与待插入磁绳的螺度符号一致。

该计算仍是薄环平衡诊断，不能取代插绳后的一步 CT 连续性测试或完整 MHD 稳定性
测试。它的作用是把插入磁通的扫描范围限制在与松弛背景磁场相容的量级。

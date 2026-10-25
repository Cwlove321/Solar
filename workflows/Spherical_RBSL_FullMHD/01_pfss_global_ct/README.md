# 全球磁图 PFSS-CT

该目录是 `Spherical_RBSL_FullMHD` 的可选第 1 阶段。它使用完整日面径向磁图和
AMRVAC 原有 `mod_pfss` 球谐展开，不调用局部 PDFI_SS。输出 DAT 可直接作为
`02_relax` 的背景场初态。

## 算法

1. 读取 `potential_boundary/pfss_boundary.dat`，格式为两个 `int32`
   (`nphi, ntheta`)，随后是 theta、phi 坐标和 `Br(nphi,ntheta)` 的 `float64`。
2. 在局部模拟角网格上根据 `path.dat`、`Bflux` 和 `Rmfr` 计算 RBSL 的底边界
   `Br_RBSL`。
3. 使用边界产品 `parameters.json` 中的零基索引 `x0/y0`，执行
   `Br_background = Br_observed - Br_RBSL`。全球 theta 存储方向与局部数组相反，
   代码显式使用 `br_rbsl(m+1-itheta,iphi)`。
4. 从扣通量后的完整日面磁图计算球谐系数。默认 `lmax=360`、源面
   `R_s=2.5 Rsun`，与旧全球 PFSS 基线一致。
5. 在均匀真实半径壳层上调用 `pfss_vector_potential()`，以四点三次插值把矢势
   映射到目标拉伸网格的 CT 边，再由 AMRVAC 离散旋度生成面心磁场。
6. `special_refine_grid()` 只保护物理底边界为最高层，当前为L4；其他块保留标准
   AMR误差判据给出的细化/粗化请求，并受树相容性约束。2026-10-22已取消非底部
   固定L3的测试限制，回调与其他四阶段一致。`itfixgrid=0` 仍冻结演化期的树；
   `usr_before_main_loop` 同时输出面心CT散度和中心磁场散度。

此源码修改需要重新编译全球PFSS并重新生成初态；不改变已有DAT的网格树，也不
自动开启演化期间的动态AMR。

这种做法的关键不是直接插值 `B`，而是插值矢势 `A` 后在目标网格离散计算
`B=curl(A)`。因此 CT 面磁通在拉伸网格上保持自洽。

## 准备输入

直接复制当前事件的全球磁图和二进制轴线：

```bash
cp /path/to/boundary_product/pfss_boundary.dat potential_boundary/
cp /path/to/rbsl_path/path.dat ./path.dat
```

`path.dat` 由 `Preprocess` 写成小端 `float64` 二进制文件，记录连续排列为：

```text
lon_1, lat_1, radius_1,
lon_2, lat_2, radius_2, ...
```

`mod_usr.t` 使用 little-endian unformatted stream 直接读取，并要求文件大小严格等于
`npoints * 3 * 8` 字节。当前 `npoints=6000` 时应为 `144000` 字节；旧的 ASCII
`path.dat` 不再适用于本目录。

新事件需要从边界产品 `parameters.json` 手工同步以下参数到 `amrvac.par`：

- `xprobmin/max2/3`；
- `domain_nx2/3`；
- `global_patch_x0/y0`；
- `global_expected_nphi/ntheta`；
- `Bflux`、`Rmfr`、`npoints`、螺度和轴线方向。

当前 rbsl14 的正确全球输入是 `3600 x 1800`、`0.1 deg` 磁图。旧对照目录中的
`1800 x 900`、`0.2 deg` 磁图不能代替它。

## 初始大气配置（2026-10-22 更新）

本阶段与 `01_pfss`、`02_relax`、`03_insert`、`04_fullmhd` 及其 profile 共用同一
`&atmosphere_list` 接口。默认是解析过渡区大气（`analytic_tr`：双曲正切低层＋恒
导热通量上层），并按 `n_H(10 Mm)=6.2e8 cm^-3` 归一化。旧压力锚点的插值分母
已改为位置差，不再使用压力差。

在工作流根目录运行：

```bash
/opt/anaconda3/bin/python scripts/configure_case.py case.json --atmosphere-only
/opt/anaconda3/bin/python scripts/validate_workflow.py
```

只同步大气时，JSON 仅需包含 `atmosphere`，不要求 `boundary_dir`、`path_dat` 或
局部 PDFI 输入；此模式也不会复制磁图、轴线或改动本阶段的全球网格／patch 参数。
未提供 `atmosphere` 时采用新的默认值。正常完整配置同样会同步本阶段的大气，但
不会代替上节列出的全球 PFSS 事件参数配置。

支持 `model="analytic_tr"` 或 `"isothermal"`，以及独立的
`normalization="number_density"` 或 `"pressure"`。压力模式需显式指定参考高度和
参考压力；修正后的结果不兼容旧错误插值的大气数值。`isothermal` 仅指初始温度
剖面恒定，不关闭能量方程或强制后续温度不变。没有内置 VAL-C、C7 或外部大气表。
完整参数映射、JSON 示例和物理单位见[工作流说明 §3.3.2](../README.md#332-选择大气模型与归一化)。

本次源码更新后必须重新编译；此后仅改已支持的大气 namelist 参数无需再次编译。
换大气需要协调初态与所有后续阶段的底边界：**旧 DAT 内部不会因修改参数而自动
更新**。不要只更换后续 MHD 的大气参数、再原位重启旧 DAT 来声称完成大气切换。
初始静力构造不等于热平衡，应核对启动日志的模型、归一化与实际锚点，并检查首帧
温度／密度／压力及底边界，再进行短时测试。本节所述大气更新本身不修改辐射、
加热或磁绳；其后经用户确认取消了非底部固定L3规则，见前述算法说明。
这些修改不表示已完成新的全MHD生产验证。

## 编译与运行

在 guoguo 上：

```bash
export AMRVAC_DIR=/share/home/guoyang1/codes/amrvac
export GLOBAL_FMHD_LIBDIR=/share/home/guoyang1/apps/global_fmhd_deps/chenpf_libs
$AMRVAC_DIR/setup.pl -d=3
make ARCH=intel -j 8
bsub < job
```

模板使用 640 核，且 `span[ptile=64]` 保证整节点分配。核数只能改成 64 的倍数。
模板仍保留较大的 `max_blocks=70000` 和 `autoconvert=.false.`，不在线生成VTU。
`max_blocks` 是每个MPI进程的块数容量，不是全域叶块数；旧的“全域L4有55,296块、
约9 GiB/帧”是另一种网格假设，不能作为当前仅保护底边界最高层的输出大小预测。
实际层级、叶块数及DAT大小应以初始化日志与输出文件为准。

## 输出与检查

主要输出：

```text
data/pfss_global_ct_0000.dat
potential_boundary/pfss_boundary_subtracted.dat
potential_boundary/pfss_boundary_subtracted.dat.coef
python/global_boundary_subtraction.json
python/initial_divb_diagnostics.json
```

`global_boundary_subtraction.json` 记录全球/局部尺寸、patch 索引、输入和输出 Br
范围以及 RBSL Br 范围。`initial_divb_diagnostics.json` 分别报告：

- CT 面心场的 `max/RMS divB` 和归一化散度；
- 单元中心场的 `max/RMS divB` 和归一化散度。

比较全球 PFSS-CT 与 PDFI-CT 时应优先比较面心指标。中心场散度还包含面到中心
重建和诊断离散误差。

## 参数约束

- `global_patch_x0/y0` 是 Python/产品中的零基索引，不是 Fortran 一基索引。
- `path.dat` 必须是交错排列的三个小端 `float64`，不能使用格式化文本。
- `do_rbsl=.true.` 表示在 PFSS 外推前扣除 RBSL 底边界通量，并不在 PFSS DAT 中
  插入磁绳。
- 更改磁图、`Bflux`、`Rmfr`、路径、`lmax` 后必须重建扣通量磁图和球谐系数。
  代码在每次生成扣通量磁图时主动删除同名旧 `.coef`，避免缓存污染。
- 第一轮验证完成前保持 `time_max=0`，不要直接提交后续长时间 MHD。

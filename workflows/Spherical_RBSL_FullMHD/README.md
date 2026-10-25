# 球坐标 PDFI/PFSS、RBSL 与全热力学 MHD 工作流

本目录提供基于 `rbsl14` 模型整理的球坐标源码工作流及可交付案例生成入口。它不是旧
`Spherical_RBSL` 的替代品，而是一个独立的生产工作流：使用局部高分辨率磁图，
通过 PDFI_SS/FISHPACK 构造与 constrained transport（CT）一致的势场，在不含磁绳
的背景场上先做 MHD 松弛，随后以 CT 面心形式插入 RBSL 磁绳和冷致密物质，最后
运行半相对论、含内能的 MHD 演化。

归档基线日期为 2026-07-24。模板数值来自 2024-10-07 22:00 UT 事件的 `rbsl14`，
但所有事件参数都必须重新配置，不能直接把模板数值用于新事件。新用户通过
`case.delivery.example.json`和`workflow.py new`生成独立案例；旧`case.json`
原位接口仅作兼容保留。

2026-10-22 的大气更新修正了压力锚点的空间插值，并将默认大气归一化恢复为成功
RBSL15 使用的 `n_H(10 Mm)=6.2e8 cm^-3`；同时增加可选等温初始剖面和显式压力
锚点。**本次只更新大气构造与配置接口，不把整个模板改成 RBSL15**：磁绳、质量
加载、AMR、辐射与背景加热的既有默认值均保持不变。大气配置方法见 §3.3.2，旧
算例迁移注意事项见 §3.3.3。本次源码更新后必须重新编译；不能沿用旧可执行文件。

### 逐参数说明索引（2026-10-22）

每份 `.par` 都配有同名 `.par.readme`，按原文件顺序逐项说明当前值、含义、单位、
生效条件和容易混淆之处，共覆盖当前 7 份文件的 682 项赋值。数组在原文件中的
一次赋值算一项；不是只说明各阶段之间的差异。

| 配置 | 逐参数说明 | 项数 |
| --- | --- | --- |
| 局部 PDFI/PFSS 初始化 | [01_pfss/amrvac.par.readme](01_pfss/amrvac.par.readme) | 88 |
| 全球球谐 PFSS 初始化 | [01_pfss_global_ct/amrvac.par.readme](01_pfss_global_ct/amrvac.par.readme) | 91 |
| 背景 MHD 松弛 | [02_relax/amrvac.par.readme](02_relax/amrvac.par.readme) | 95 |
| 插绳与质量加载 | [03_insert/amrvac.par.readme](03_insert/amrvac.par.readme) | 105 |
| 默认 full-MHD 续算 | [04_fullmhd/amrvac.par.readme](04_fullmhd/amrvac.par.readme) | 101 |
| 有辐射、有加热 profile | [amrvac_full_thermo.par.readme](04_fullmhd/profiles/amrvac_full_thermo.par.readme) | 101 |
| 无辐射、有加热 profile | [amrvac_norad_with_heating.par.readme](04_fullmhd/profiles/amrvac_norad_with_heating.par.readme) | 101 |

本次按用户要求从 `04_fullmhd` 及两份 profile 移除低层辐射衰减三项
`rad_damp`、`rad_damp_height`、`rad_damp_scale`；其他辐射设置及物理、网格、
时间参数保持不变。当前核验的共享 AMRVAC 源码默认 `rad_damp=false`，因此
有辐射 profile 不再施加这类人为低层抑制；默认和无辐射 profile 仍关闭辐射。
没有改动底层求解器源码、已生成 DAT 或已提交模拟。

说明文件附有对应 `.par` 的 SHA256。之后运行 `configure_case.py`、切换 profile
或手动修改参数时，应同步说明中的当前值；说明文件本身不参与程序运行。

### 可交付版本入口（2026-10-22）

新用户请先阅读[QUICKSTART.md](QUICKSTART.md)，依赖安装与源码来源见
[DEPENDENCIES.md](DEPENDENCIES.md)，验收范围与结果见[VALIDATION.md](VALIDATION.md)。
推荐入口是根目录`workflow.py`：

```bash
python workflow.py new --boundary /data/boundary --axis /data/axis/path.dat \
  --settings case.delivery.example.json --output /data/runs/my_event
python workflow.py check /data/runs/my_event --config-only
python workflow.py machine /data/runs/my_event --profile machine.json
python workflow.py check /data/runs/my_event
python workflow.py build /data/runs/my_event
python workflow.py run /data/runs/my_event
python workflow.py status /data/runs/my_event
```

先确认科学参数，并按接收方机器修改`machine.example.json`另存为`machine.json`。
新入口在独立目录生成案例，不覆盖模板或原始观测产品。该目录带输入快照、
SHA256清单、解析后的配置、机器配置、各阶段参数及逐参数说明。

交付版增加了以下衔接：

- 从实际预处理数组和metadata推导角域、角向根网格及轴线点数，并核对FITS单位/
  方向、全球patch、有限值、轴线范围、参考时刻及网格块整除。
- 同步局部/全球PFSS和后续阶段的网格、磁通、半径、旋向及轴线；选择无绳链时
  自动关闭PFSS扣除磁绳Br并跳过插绳阶段。
- 由`case.delivery.example.json`显式选择固定/动态AMR、04阶段热物理方案、质量加载与
  各阶段时间。底边界最高层保护保留，没有事件专用box或非底部固定L3限制。
- 根据机器JSON生成编译链接与启动文件，支持local和LSF；核数、节点布局、编译器
  和库路径由接收机器决定。Slurm等其他调度器尚未实现。
- 阶段完成后按DAT真实代码时间与结构/状态验收交接；检查所有状态和CT面数据的
  有限性、密度/内能正性及故障日志，不依赖固定末帧编号或调度器DONE。
- 保存已验收状态，校验配置、输入、可执行文件与DAT后才能跳过成功阶段；失败或
  未登记输出不会被自动覆盖，避免重复插绳/加载。机器配置在执行后冻结。

`vendor/amrvac`包含当前工作流所用共享求解器的源码快照和本地修改，
`vendor/simesh`提供可选分析源码。版本、来源、许可证和逐文件哈希见
`vendor/SOURCE_PROVENANCE.md`、`vendor/dependencies.lock.json`；它不是历史RBSL15
某个独立生产二进制的逐位归档。PDFI_SS、FISHPACK、CFITSIO仍须在接收机器安装，
源码接口/命令/库文件预检通过之后，真实构建是否成功以`build`结果为准。

交付验收使用合成输入、微型v5 DAT和受控假求解器验证配置、故障阻断与交接；
这不替代接收者机器上的真实MPI编译及新事件物理验证。测试命令：

```bash
python -m pip install -r requirements-test.txt
python -m pytest -q scripts
python scripts/verify_dependencies.py
```

以下章节保留物理模型与旧手工流程参考。旧`configure_case.py`仍是原模板的原位
配置接口，旧shell脚本仍包含集群/帧号约定；新用户按QUICKSTART使用统一入口。

## 1. 工作流结构

```text
Spherical_RBSL_FullMHD/
├── Preprocess/                    # HMI/AIA、边界、轴线和 PDFI FITS 预处理
├── 01_pfss/                       # PDFI_SS 局部势场，先扣除 RBSL 的底边界 Br
│   ├── Bphot/Brll0000.fits        # 无鬼点，数组顺序为 Python (phi, theta)
│   ├── path.dat                   # RBSL 轴线，用于计算待扣除的 Br_RBSL
│   ├── mod_usr.t
│   ├── amrvac.par
│   ├── local.make
│   └── job
├── 01_pfss_global_ct/             # 可选：全球磁图球谐 PFSS，矢势映射到拉伸 CT 边
│   ├── potential_boundary/        # 完整日面 pfss_boundary.dat
│   ├── path.dat                   # 小端 float64 二进制 RBSL 轴线
│   ├── mod_usr.t
│   ├── amrvac.par
│   ├── local.make
│   └── job
├── 02_relax/                      # 不含磁绳的背景 MHD 松弛
├── 03_insert/                     # 插入磁绳，可选暗条增密并演化到交接时刻
├── 04_fullmhd/                    # continuation，可接 02 或 03 的输出
│   └── profiles/                  # 可选热物理配置
├── scripts/
│   ├── configure_case.py          # 同步事件参数及五阶段/profile 的大气参数
│   ├── build_all.sh               # 统一编译四阶段
│   ├── link_restarts.sh           # 建立阶段间 DAT 软链接
│   ├── submit_chain.sh            # 可选 LSF 依赖链
│   ├── upload_vtu_to_yangyang.sh  # 本机流式中转 VTU
│   ├── spherical_slice_analysis.py # 固定球坐标/自由平面统一切片分析
│   ├── README_spherical_slice_analysis.md
│   └── validate_workflow.py       # 静态一致性检查
├── case.example.json
└── README.md
```

默认依赖关系为：

```text
HMI vector/synoptic + AIA
             |
             v
        Preprocess
             |
      Brll0000.fits + path.dat
             |
             v
   01_pfss: Br_bg = Br_obs - Br_RBSL
             |
             v
   02_relax: background MHD relaxation
        |                         |
        | no flux rope            | insert flux rope
        |                         v
        |              03_insert: curl(A_RBSL)
        |              optional instant/Liakh density
        |              optional finite evolution
        |                         |
        +-------------------------+
                     |
                     v
   04_fullmhd: continuation and eruption evolution
```

`01_pfss_global_ct` 可替换上图的 `01_pfss`。它从完整日面磁图计算旧球谐 PFSS，
再把均匀真实半径壳层上的矢势插值到拉伸网格 CT 边；详细用法见
[`01_pfss_global_ct/README.md`](01_pfss_global_ct/README.md)。

## 2. 相比原成熟代码的重大修改

### 2.1 PDFI 势场与 CT 离散一致性

1. 当前的 PFSS 不是完整 PDFI 电场反演，而是调用 PDFI_SS/FISHPACK，在局部球面
   楔形区域内从单时刻底边界 `Br` 求势场。
2. PDFI 返回网格边上的 `A_theta` 和 `A_phi`，磁场由离散环量除以面面积得到，
   与 AMRVAC 的 CT 面心磁场使用同一几何离散。
3. 拉伸径向网格不再用最近邻或简单线性壳层映射，而是对均匀 `log(r)` 的 PDFI
   壳层做四点三次重构。高径向层数下必须做 `p=288/576` 收敛测试，因为测试曾发现
   高空衰减速度会随径向离散显著改变。
4. 构造保持观测底边界的插绳初态时，先计算 `Br_RBSL`，再用
   `Br_background = Br_observed - Br_RBSL` 做势场外推。最后插入完整 RBSL 后，
   `Br_total` 才近似回到观测值。

### 2.2 restart 插绳的关键 CT 修复

原始 restart 插绳分支把 `curl(A_RBSL)` 加到临时面心数组 `wstmp`，却只用它更新
单元中心磁场，没有写回真正由 CT 推进和 DAT 保存的 `block%ws`。结果是：

- `0000` 可视化中能看到磁绳；
- 第一个 CT 时间步从旧背景 `block%ws` 重建中心场；
- `0001` 中磁绳瞬间消失。

归档代码在 `initonegrid_usr()` 中增加：

```fortran
block%ws(ixGs^S,1:nws) = wstmp(ixGs^S,1:nws)
```

并保证 `b_from_vector_potential_S()` 使用加法：

```fortran
ws(ixC^S,idim1) = ws(ixC^S,idim1) + circ(ixC^S,idim1)
```

而不是用磁绳场覆盖松弛后的背景场。`usr_improve_initial_condition` 也被显式接管，
避免默认初始化器再次用单独的 RBSL 矢势重建总磁场。一步测试已确认 `0000 -> 0001`
的闭合磁绳截面和三个磁场分量连续。

### 2.3 两种暗条增密方法

#### 2.3.1 插绳时一次性增密

`03_insert` 从松弛 DAT 初始化并插入磁绳时，可在轴线附近同时加载暗条物质。
设到轴线的最短距离为 `d`，
`xi=d/a`，核心比例为 `f_c`，密度窗函数为

```math
W(\xi)=
\begin{cases}
1, & 0\le\xi\le f_c,\\
\frac12\left[1+\cos\frac{\pi(\xi-f_c)}{1-f_c}\right],
& f_c<\xi<1,\\
0, & \xi\ge1.
\end{cases}
```

加载后的密度为

```math
\rho'=\rho\left[1+(\chi-1)W(\xi)\right],
```

其中默认 `chi=20`、`f_c=0.6`，并且高度低于 5 Mm 不加载。插入时保持压强和内能
不变，因此核心初始温度约变为环境温度的 `1/chi`。若打开光学薄辐射后发生过快
冷却，应先把 `chi` 降到 10 到 15，而不是任意改变磁绳磁通。

#### 2.3.2 Liakh 时间依赖质量源

`03_insert` 另包含 Liakh, Luna & Khomenko (2023), A&A 673, A154 式 (8)
的质量源，并把论文的二维直磁绳分布推广到观测得到的三维弯曲 RBSL 轴线：

```math
S_\rho=\frac{\chi\rho_0(r)}{t_{\rm load}}
\exp\left[-\left(\frac{d_h}{\sigma_h}\right)^4
-\left(\frac{d_r}{\sigma_r}\right)^4\right],\qquad
t_{\rm start}\le t<t_{\rm start}+t_{\rm load}.
```

其中 `rho0(r)` 取初始化静力大气表 `ra(r)`，不取随时间变化的当前密度，因而是
线性累计注入而不是指数增长；`d_h` 是单元方向与最近轴线方向之间的球面弧长，
`d_r` 是单元与该轴线点的径向距离。代码在启动时把最近轴线半径和夹角预计算到
最高层角网格，Runge-Kutta 子步只做二维查表，不在每个三维单元内搜索全部轴线点。

默认参数复现 rbsl14 已跑通的 100 s 加载设置：

```text
liakh_mass_loading        = .true.
liakh_density_contrast    = 30
liakh_t_start_s           = 0 s
liakh_t_load_s            = 100 s
liakh_sigma_horizontal_Mm = 3.6 Mm
liakh_sigma_radial_Mm     = 4.8 Mm
liakh_min_height_Mm       = 5 Mm
```

源项只增加守恒密度，不增加动量或总能量，等价于注入静止冷物质；压强不会随注入
同步增加。时间参数以物理秒给出，代码内部用启动日志中的 `unit_time` 转为无量纲时间。
严格积分后轴心累计增加 `chi*rho0`，所以若初始密度恰为 `rho0`，`chi=30` 最终是
`31 rho0`；若要求最终严格为背景的 30 倍，应设 `liakh_density_contrast=29`。

两种方法不能同时打开。使用 Liakh 源时，`03_insert` 设置
`load_filament_density=.false.`：初始化只插入磁绳，随后由源项逐步加密；使用瞬时
方法时则设置 `load_filament_density=.true.`、`liakh_mass_loading=.false.`。
`03_insert` 必须有与 PFSS 扣通量相同的二进制 `path.dat`。Liakh 方法与辐射开关
相互独立；rbsl14 的验证算例关闭辐射、保留热传导、TRAC、背景加热、重力和电阻率。

### 2.4 参数化热物理和电阻率

- 经验背景加热由硬编码改成 `background_heating_cgs` namelist 参数。
- 当前实际启用的背景加热只有指数项

```math
H(r)=H_0\exp\left[-\frac{r-r_{\min}}{H_s}\right],\qquad H_s=60\ {\rm Mm}.
```

  原论文中依赖 `B`、密度和磁力线曲率的另一项仍保留在源码注释中，但本工作流
  没有启用，不能在论文中声称已经使用完整的双分支经验加热公式。
- 显式磁扩散率改为论文值 `eta=200 km^2 s^-1`，代码单位为
  `2.4693930682273452e-6`。源码中存在可选的异常电阻函数，但回调没有注册，
  默认计算使用空间常数 `mhd_eta`。
- 光学薄辐射、低层辐射阻尼、热传导、TRAC 和背景加热都由独立参数控制。

### 2.5 拉伸网格、动态 AMR 和边界外推

- 保留球坐标径向拉伸网格，并为拉伸网格单独计算一侧三阶等梯度外推系数，
  不再把均匀网格系数直接用于非均匀 `dr`。
- 代码支持动态 AMR：`refine_criterion=3`，仅密度参与 Lohner 误差，重网格周期为
  每 3 步；接触 `r_min` 的底边界块始终保持最高层。**但当前默认 `.par` 使用
  `itfixgrid=0` 冻结网格树**，不能据上述误差与周期参数就认定运行中会重网格。
  本次大气更新没有同步 RBSL15 的动态 AMR 或强制 L4 盒子。
- 当前没有事件专用的经纬高L4盒，也没有历史故障区专用盒；底边界最高层保护
  不属于这种坐标盒。全球PFSS的非底部固定L3测试限制亦已取消，五阶段使用相同
  的底边界保护回调，详见前面的交付审计。
- 四个 job 强制 MPI 核数为 64 的倍数，并使用 `span[ptile=64]`，不允许把同一
  64 核节点拆成两个 32 核作业。

#### 2.5.1 径向拉伸的实际定义和分辨率

四个阶段使用同一组 `&meshlist` 参数，因此 PFSS、松弛、插绳和 full MHD 的物理
坐标完全对应。计算域是球坐标楔形：

```math
r\in[1.0,1.5]R_\odot,\qquad
\theta\in[0.16,0.24]\ {\rm u.a.},\qquad
\phi\in[0.4616666667,0.6216666667]\ {\rm u.a.}.
```

这里 `u.a.`（unit angle，单位角度）是 SAMRVAC 的角坐标单位：

```math
1\ {\rm u.a.}=360^\circ=2\pi\ {\rm rad}.
```

基础网格为

```math
(N_r,N_\theta,N_\phi)=(72,36,72),\qquad
(n_{\rm block,r},n_{\rm block,\theta},n_{\rm block,\phi})=(12,12,12).
```

只有径向开启单边等比拉伸：`stretch_dim(1)='uni'`，
`qstretch_baselevel(1)=1.041541581503301`；两个角向均匀。设

```math
q_1=1.041541581503301,
\Delta r_{1,0}=(r_{\max}-r_{\min})\frac{1-q_1}{1-q_1^{N_r}},
```

则一级（L1）的径向节点与单元宽度为

```math
r_i=r_{\min}+\Delta r_{1,0}\frac{q_1^i-1}{q_1-1},\qquad
\Delta r_{1,i}=\Delta r_{1,0}q_1^i,\qquad i=0,\ldots,N_r-1.
```

这不是“高度坐标均匀、仅绘图时拉伸”，而是求解器实际使用的 \(r\) 坐标：单元从
底部向顶部逐渐变宽。取 \(R_\odot=695.5\ {\rm Mm}\)，L1 的径向单元宽度约从

```math
\Delta r_{1,0}=0.814\ {\rm Mm}
\quad\text{增至}\quad
\Delta r_{1,71}=14.65\ {\rm Mm}.
```

AMR 的第 \(\ell\) 层并非简单复制 L1 的几何级数，而是使用

```math
q_\ell=q_1^{1/2^{\ell-1}},\qquad
\Delta r_{\ell,0}=
\frac{\Delta r_{\ell-1,0}}{1+q_\ell},
\qquad N_{r,\ell}=72\times2^{\ell-1}.
```

因此每个父单元恰好由两个子单元填满，粗细网格交界处的物理坐标连续。最高的 L4
具有 \(q_4=1.00510070\)，底部和顶部的局部径向宽度约为

```math
\Delta r_{4,0}=0.100\ {\rm Mm},\qquad
\Delta r_{4,575}=1.864\ {\rm Mm}.
```

角向不拉伸；以下数值首先以单位角度给出：

```math
\Delta\theta_{\rm L1}=\Delta\phi_{\rm L1}
=\frac{0.08}{36}\ {\rm u.a.}
=\frac{0.16}{72}\ {\rm u.a.}
=2.222222\times10^{-3}\ {\rm u.a.}=0.8^\circ,
```

```math
\Delta\theta_{\rm L4}=\Delta\phi_{\rm L4}
=\frac{2.222222\times10^{-3}}{2^3}\ {\rm u.a.}
=2.777778\times10^{-4}\ {\rm u.a.}=0.1^\circ.
```

将单位角度转换为弧度后，物理横向尺度为

```math
\Delta s_\theta=2\pi r\Delta\theta,
\qquad
\Delta s_\phi=2\pi r\sin(2\pi\theta)\Delta\phi.
```

底边界块被强制维持 L4，故其角向间距始终为 \(0.1^\circ\)，最小径向间距约为
\(0.1\ {\rm Mm}\)；其余区域由密度 Lohner 指标决定是否达到 L4。

**RBSL14 边界磁图与 L4 网格一致。** RBSL14 的基准边界产品
`rbsl_2200_lon-13.8_43.8_lat3.6_32.4_res01_v1/parameters.json` 确实给出
`resolution_degree=0.1`，其 `Brll0000.fits` 为

```math
(N_{\rm lat},N_{\rm lon})=(288,576),\qquad
\Delta\mathrm{lat}=\Delta\mathrm{lon}=0.1^\circ.
```

这表示输入磁图覆盖 \(28.8^\circ\times57.6^\circ\)，即
\(0.08\times0.16\ {\rm u.a.}\)，恰好对应当前 `xprob` 的角宽。代码以

```math
(m,n)=(36\times2^3,72\times2^3)=(288,576)
```

在 `read_b()` 中将 FITS 数组直接读入 `brll0(m,n)`；数组尺寸与 L4 角网格逐点相同，
无需重采样。`amrvac.par` 的 L4 间距为

```math
\Delta\theta=\frac{0.24-0.16}{288}=2.777778\times10^{-4}\ {\rm u.a.}=0.1^\circ,
\qquad
\Delta\phi=\frac{0.6216666667-0.4616666667}{576}=2.777778\times10^{-4}\ {\rm u.a.}=0.1^\circ.
```

因此，当前归档中 \(0.1^\circ\) 的 RBSL14 底边界磁图与 L4 数值底边界逐点对应，
两者的角向分辨率一致。

**后处理约束。** MPI-AMRVAC 的 DAT 文件头不记录 `stretch_dim` 或
`qstretch_baselevel`。读取 DAT、构造固定高度切片、插值磁绳轴线或计算体积积分时，
必须同时提供对应的 `amrvac.par`，或显式使用上述映射重建每个 leaf block 的真实节点。
若按均匀 \((r_{\max}-r_{\min})/N_r\) 重建，会将数组索引对应到错误高度，尤其会污染
低层的温度/密度切片和 AMR 交界附近的诊断。

### 2.6 归一化初始化顺序

`unit_length`、`unit_numberdensity` 和 `unit_temperature` 必须在 `mhd_activate()`
之前设置，否则半相对论约化光速和派生单位会使用错误尺度。当前基准为：

```text
unit_length        = 6.955e10 cm
unit_numberdensity = 1.0e9 cm^-3
unit_temperature   = 1.0e6 K
```

由当前组分定义得到的典型派生量约为：

```text
unit_density       = 2.341668e-15 g cm^-3
unit_pressure      = 0.317538 erg cm^-3
unit_magneticfield = 1.997574 G
unit_velocity      = 1.164488e7 cm s^-1
unit_time          = 5.9726e3 s
```

原代码中 `unit_time=85.87 s` 是旧长度单位遗留注释，已在本归档修正。正式换算时
始终以 AMRVAC 启动日志实际打印的单位为准。

## 3. 控制方程

### 3.1 半相对论 MHD

质量守恒：

```math
\frac{\partial\rho}{\partial t}+\nabla\cdot(\rho\mathbf v)=0.
```

采用 Boris 型半相对论修正的动量方程。写成电磁动量守恒形式为

```math
\frac{\partial}{\partial t}
\left(\rho\mathbf v+\frac{\mathbf E\times\mathbf B}{\mu_0c_{\rm red}^{2}}\right)
+\nabla\cdot\left[
\rho\mathbf v\mathbf v+p\mathbf I
+\frac{B^2}{2\mu_0}\mathbf I-\frac{\mathbf B\mathbf B}{\mu_0}
+\frac{E^2}{2\mu_0c_{\rm red}^{2}}\mathbf I
-\frac{\mathbf E\mathbf E}{\mu_0c_{\rm red}^{2}}
\right]
=\rho\mathbf g+\mathbf S_{\rm Boris},
```

为便于逐项阅读和截图，同一方程可写成带物理含义标注的形式：

```math
\frac{\partial}{\partial t}\left[
\underbrace{\rho\mathbf v}_{\text{流体动量}}
+\underbrace{\frac{\mathbf E\times\mathbf B}{\mu_0c_{\rm red}^{2}}}_{\text{电磁惯性（Boris 修正）}}
\right]
+\nabla\cdot\left[
\underbrace{\rho\mathbf v\mathbf v}_{\text{惯性动量通量}}
+\underbrace{p\mathbf I}_{\text{气压通量}}
+\underbrace{\frac{B^2}{2\mu_0}\mathbf I}_{\text{磁压通量}}
-\underbrace{\frac{\mathbf B\mathbf B}{\mu_0}}_{\text{磁张力通量}}
+\underbrace{\frac{E^2}{2\mu_0c_{\rm red}^{2}}\mathbf I}_{\text{电场各向同性应力}}
-\underbrace{\frac{\mathbf E\mathbf E}{\mu_0c_{\rm red}^{2}}}_{\text{电场各向异性应力}}
\right]
=\underbrace{\rho\mathbf g}_{\text{重力}}
+\underbrace{\mathbf S_{\rm Boris}}_{\text{Boris 源项}}.
```

English labels for the same equation are:

```math
\frac{\partial}{\partial t}\left[
\underbrace{\rho\mathbf v}_{\text{fluid momentum}}
+\underbrace{\frac{\mathbf E\times\mathbf B}{\mu_0c_{\rm red}^{2}}}_{\text{electromagnetic inertia (Boris correction)}}
\right]
+\nabla\cdot\left[
\underbrace{\rho\mathbf v\mathbf v}_{\text{inertial momentum flux}}
+\underbrace{p\mathbf I}_{\text{gas-pressure flux}}
+\underbrace{\frac{B^2}{2\mu_0}\mathbf I}_{\text{magnetic-pressure flux}}
-\underbrace{\frac{\mathbf B\mathbf B}{\mu_0}}_{\text{magnetic-tension flux}}
+\underbrace{\frac{E^2}{2\mu_0c_{\rm red}^{2}}\mathbf I}_{\text{isotropic electric-field stress}}
-\underbrace{\frac{\mathbf E\mathbf E}{\mu_0c_{\rm red}^{2}}}_{\text{anisotropic electric-field stress}}
\right]
=\underbrace{\rho\mathbf g}_{\text{gravity}}
+\underbrace{\mathbf S_{\rm Boris}}_{\text{Boris source term}}.
```

其中

```math
\mathbf E=-\mathbf v\times\mathbf B,\qquad
c_{\rm red}=6000\ {\rm km\ s^{-1}}.
```

以上记号中的 Boris 源项为

```math
\mathbf S_{\rm Boris}=\frac{1}{\mu_0}
\left(\frac{1}{c_0^2}-\frac{1}{c_{\rm red}^2}\right)
\left[\mathbf E(\nabla\cdot\mathbf E)
-\mathbf E\times(\nabla\times\mathbf E)\right],
```

其中 `c0` 是真实光速。AMRVAC 内部的 Boris 修正保证约化光速系统与原 MHD 极限
一致。半相对论限制的是阿尔芬特征速度及相应时间步，不会限制
由未平衡洛伦兹力产生的真实流速，因此它不能代替磁场离散一致性和背景松弛。

磁感应方程为

```math
\frac{\partial\mathbf B}{\partial t}
=\nabla\times\left(\mathbf v\times\mathbf B-\eta\mathbf J\right),\qquad
\mathbf J=\frac{1}{\mu_0}\nabla\times\mathbf B,\qquad
\nabla\cdot\mathbf B=0.
```

逐项标注的中文形式为：

```math
\underbrace{\frac{\partial\mathbf B}{\partial t}}_{\text{磁场的时间演化}}
=\nabla\times\left[
\underbrace{\mathbf v\times\mathbf B}_{\text{理想感应（磁通冻结）}}
-\underbrace{\eta\mathbf J}_{\text{电阻扩散}}
\right],\qquad
\underbrace{\mathbf J}_{\text{电流密度}}
=\underbrace{\frac{1}{\mu_0}\nabla\times\mathbf B}_{\text{由磁场旋度计算}},
\qquad
\underbrace{\nabla\cdot\mathbf B}_{\text{无散磁场约束}}=0.
```

The same annotated induction equation in English is:

```math
\underbrace{\frac{\partial\mathbf B}{\partial t}}_{\text{magnetic-field time evolution}}
=\nabla\times\left[
\underbrace{\mathbf v\times\mathbf B}_{\text{ideal induction (flux freezing)}}
-\underbrace{\eta\mathbf J}_{\text{resistive diffusion}}
\right],\qquad
\underbrace{\mathbf J}_{\text{current density}}
=\underbrace{\frac{1}{\mu_0}\nabla\times\mathbf B}_{\text{computed from magnetic curl}},
\qquad
\underbrace{\nabla\cdot\mathbf B}_{\text{solenoidal magnetic-field constraint}}=0.
```

数值上使用交错网格 CT 更新面中心磁场，使离散磁通守恒达到机器精度。CT 只能
控制离散散度，不能保证 `J x B` 的离散误差为零。

### 3.2 内能方程和状态方程

启用 `mhd_energy=.true.` 和 `mhd_internal_e=.true.`，使用双能量策略。内能密度
`epsilon=p/(gamma-1)` 满足

```math
\frac{\partial\epsilon}{\partial t}
+\nabla\cdot(\epsilon\mathbf v)
=-p\nabla\cdot\mathbf v
+\nabla\cdot\left(\boldsymbol\kappa\cdot\nabla T\right)
-n_en_H\Lambda(T)
+H(r)
+\eta J^2.
```

逐项标注的中文形式为：

```math
\underbrace{\frac{\partial\epsilon}{\partial t}}_{\text{内能的时间变化}}
+\underbrace{\nabla\cdot(\epsilon\mathbf v)}_{\text{内能的平流输运}}
=\underbrace{-p\nabla\cdot\mathbf v}_{\text{压缩加热或膨胀冷却}}
+\underbrace{\nabla\cdot\left(\boldsymbol\kappa\cdot\nabla T\right)}_{\text{热传导的净热输入}}
-\underbrace{n_en_H\Lambda(T)}_{\text{光学薄辐射损失}}
+\underbrace{H(r)}_{\text{经验背景加热}}
+\underbrace{\eta J^2}_{\text{焦耳加热}}.
```

The same annotated internal-energy equation in English is:

```math
\underbrace{\frac{\partial\epsilon}{\partial t}}_{\text{internal-energy time change}}
+\underbrace{\nabla\cdot(\epsilon\mathbf v)}_{\text{advective transport of internal energy}}
=\underbrace{-p\nabla\cdot\mathbf v}_{\text{compressional heating or expansion cooling}}
+\underbrace{\nabla\cdot\left(\boldsymbol\kappa\cdot\nabla T\right)}_{\text{net conductive heat input}}
-\underbrace{n_en_H\Lambda(T)}_{\text{optically thin radiative loss}}
+\underbrace{H(r)}_{\text{empirical background heating}}
+\underbrace{\eta J^2}_{\text{Joule heating}}.
```

其中

```math
p=(\gamma-1)\epsilon,\qquad \gamma=5/3,
```

温度由 AMRVAC 的完全电离等离子体单位关系从 `p/rho` 得到。各项分别代表压缩
加热或膨胀冷却、场向热传导、光学薄辐射损失、经验背景加热和焦耳加热。
用于守恒推进的总能量基准为

```math
e_{\rm tot}=\frac{p}{\gamma-1}+\frac12\rho v^2+\frac{B^2}{2\mu_0},
```

半相对论形式还包含相应的电磁惯性修正。`mhd_internal_e` 在强磁场、低 beta 区域
提供内能辅助方程，避免用两个大能量之差恢复压强时丢失精度。

### 3.3 重力

```math
\mathbf g(r)=-\frac{GM_\odot}{r^2}\hat{\mathbf e}_r
=g_\odot\left(\frac{R_\odot}{r}\right)^2\hat{\mathbf e}_r,\qquad
g_\odot=-2.74\times10^4\ {\rm cm\ s^{-2}}.
```

初始大气通过一维静力方程

```math
\frac{dp}{dr}=\rho g_r
```

与给定温度分布共同构造。底边界的 \(\rho\) 和 \(p\) 固定为初始化表中的值。

#### 3.3.1 初始大气的实际构造方法

当 `mhd_energy=.true.` 时，`initglobaldata_usr()` 调用 `inithdstatic()`。它不是在每个
三维网格点独立指定密度和压强，而是先沿径向建立一张一维大气表，再把该表映射到三维
球坐标网格。定义相对底边界的几何高度

```math
h=r-r_{\min}.
```

一维表额外向底部和顶部各延伸 `nghostcells` 个最高层级网格宽度，以便同时给物理区和
底边界鬼点提供值。默认 `atmosphere_model='analytic_tr'` 时，表中温度由以下分段
函数给出；默认参数为

```math
T_{\rm pho}=8000\ {\rm K},\qquad
h_{\rm tr}=2\ {\rm Mm},\qquad
w_{\rm tr}=0.25\ {\rm Mm},
```

```math
T_{\rm tr}=1.6\times10^{5}\ {\rm K},\qquad
F_c=2\times10^{5}\ {\rm erg\ cm^{-2}\ s^{-1}},\qquad
h_{\rm gbot}=-0.2\ {\rm Mm}:
```

```math
T(h)=
\begin{cases}
T_{\rm pho}+\dfrac{T_{\rm top}-T_{\rm pho}}{2}
\left[\tanh\!\left(\dfrac{h-h_{\rm tr}+f_{\rm tr}}{w_{\rm tr}}\right)+1\right],
& h_{\rm gbot}<h\le h_{\rm tr},\\
\left[T_{\rm tr}^{7/2}+\dfrac{7F_c}{2\kappa_0}(h-h_{\rm tr})\right]^{2/7},
& h>h_{\rm tr},\\
T(h_{\rm gbot}), & h\le h_{\rm gbot},
\end{cases}
```

其余温度与导热系数参数为

```math
T_{\rm top}=1.5\ {\rm MK},\qquad
\kappa_0=8\times10^{-7}\ {\rm erg\ cm^{-1}\ s^{-1}\ K^{-7/2}},
```

并且

```math
f_{\rm tr}=w_{\rm tr}\operatorname{atanh}
\left[\frac{2(T_{\rm tr}-T_{\rm pho})}{T_{\rm top}-T_{\rm pho}}-1\right].
```

该选择使双曲正切层在 \(h=h_{\rm tr}\) 处恰为 \(T_{\rm tr}\)，从而与上方的恒定场向
导热通量解连续。因此，\(T_{\rm top}\) 是低层双曲正切段的幅度参数；高空温度实际由
导热通量积分决定，并非简单地在整个日冕中固定为 \(1.5\ {\rm MK}\)。

在该温度表上，重力按高度变化：

```math
g_r(h)=g_\odot\left(\frac{R_\odot}{R_\odot+h}\right)^2.
```

代码在其归一化单位下使用状态关系 \(p=\rho T\)，故先积分

```math
p(h)=C\exp\!\left[\int_{h_0}^{h}\frac{g_r(h')}{T(h')}\,dh'\right],\qquad
\rho(h)=\frac{p(h)}{T(h)}.
```

\(g_r<0\)，所以压强随高度增加而降低。积分先以任意正的临时底值开始，再将整张
\(p(h)\)、\(\rho(h)\) 表同比例缩放。默认使用数密度归一化：

```math
n_H(h_{\rm ref}=10\ {\rm Mm})=6.2\times10^8\ {\rm cm^{-3}}.
```

先在一维表上对参考高度做线性空间插值，再计算

```math
\mathrm{ratio}=\frac{n_{H,\rm ref}}{n_H(h_{\rm ref})},\qquad
\rho\leftarrow\mathrm{ratio}\,\rho,\quad p\leftarrow\mathrm{ratio}\,p.
```

这里 `reference_nh_cm3` 是氢数密度，不是质量密度的 cgs 数值；代码通过
`unit_numberdensity` 转换。也可显式选择压力归一化，此时使用
`ratio=p_ref/p(h_ref)`，密度和压力仍同比例缩放，温度不因归一化改变。

2026-10-22 已修正压力锚点的插值公式为

```math
p(h_{\rm ref})=p_{j-1}+
\frac{h_{\rm ref}-h_{j-1}}{h_j-h_{j-1}}(p_j-p_{j-1}).
```

分母使用位置差 `ya(j)-ya(j-1)`，不再使用错误的压力差。修正后的 `pressure`
选项不用于数值复现旧错误分母所产生的高密大气。运行时应核对启动日志的
`Atmosphere model`、`Atmosphere normalization`、参考高度及实际锚点，并检查
`Rho at fixed bottom boundary` 和 `P at fixed bottom boundary`。

最后，对每个三维单元，代码以余弦权重插值表中的相邻点：

```math
q(h)=q_n+\frac{1-\cos(\pi\delta h/\Delta h)}{2}(q_{n+1}-q_n),\qquad q\in\{\rho,p\}.
```

这里 \(\delta h\) 是单元高度相对表格第 \(n\) 点的偏移。底边界速度固定为零，而底部鬼点的
\(\rho\)、\(p\) 由同一一维表在 \(h<0\) 处线性插值得到并固定写入；因此底边界热力学状态不是
零梯度外推，也不是每一步重新求静力平衡。`02_relax` 从该大气开始松弛，后续 restart
阶段继续使用同一构造的底边界表。

#### 3.3.2 选择大气模型与归一化

大气模型和归一化方式彼此独立：

- `model="analytic_tr"`（默认）：上述双曲正切过渡区＋恒导热通量解析上层剖面。
- `model="isothermal"`：初始一维表的温度取 `isothermal_k`，其余仍按球坐标重力积分
  静力平衡。**这不关闭能量方程，也不把后续演化强制为等温 MHD。**
- `normalization="number_density"`（默认）：在指定高度固定氢数密度。
- `normalization="pressure"`：在指定高度固定压力；通过 JSON 配置时必须显式给出参考高度和参考压力。

没有内置 VAL-C、C7 或外部大气表读取接口。两种初始模型都只保证构造上的静力平衡，
不自动满足辐射、加热和热传导的热平衡。

`case.json` 的 `atmosphere` 对象与各阶段 `&atmosphere_list` 的映射如下。所有数值
使用表中物理单位；脚本负责生成 Fortran namelist，运行时再转换成代码单位。

| JSON 键 | Fortran 参数 | 默认值与物理单位 |
| --- | --- | --- |
| `model` | `atmosphere_model` | `analytic_tr` |
| `normalization` | `atmosphere_normalization` | `number_density` |
| `reference_height_mm` | `atmosphere_reference_height_Mm` | 10 Mm，相对底边界 |
| `reference_nh_cm3` | `atmosphere_reference_nH_cgs` | `6.2e8` cm⁻³ |
| `reference_pressure_cgs` | `atmosphere_reference_pressure_cgs` | 138.365 dyn cm⁻²，仅压力归一化使用 |
| `chromosphere_k` | `atmosphere_chromosphere_K` | 8000 K |
| `tanh_top_k` | `atmosphere_tanh_top_K` | `1.5e6` K |
| `transition_height_mm` | `atmosphere_transition_height_Mm` | 2 Mm |
| `transition_width_mm` | `atmosphere_transition_width_Mm` | 0.25 Mm |
| `transition_k` | `atmosphere_transition_K` | `1.6e5` K |
| `conductive_flux_cgs` | `atmosphere_conductive_flux_cgs` | `2e5` erg cm⁻² s⁻¹ |
| `floor_height_mm` | `atmosphere_floor_height_Mm` | −0.2 Mm，低于此高度温度保持常数 |
| `isothermal_k` | `atmosphere_isothermal_K` | `1e6` K，仅等温初始模型使用 |

默认解析模型示例（可只保留需要覆盖的字段，其他字段填入成功大气的默认值）：

```json
{
  "atmosphere": {
    "model": "analytic_tr",
    "normalization": "number_density",
    "reference_height_mm": 10.0,
    "reference_nh_cm3": 6.2e8
  }
}
```

等温初态示例：

```json
{
  "atmosphere": {
    "model": "isothermal",
    "isothermal_k": 1.0e6,
    "normalization": "number_density",
    "reference_height_mm": 10.0,
    "reference_nh_cm3": 6.2e8
  }
}
```

若确实需要底部压力归一化，应主动选择并给出锚点；以下只是配置示例，不是推荐默认，
也不保证与当前背景加热平衡：

```json
{
  "atmosphere": {
    "model": "analytic_tr",
    "normalization": "pressure",
    "reference_height_mm": 0.0,
    "reference_pressure_cgs": 138.365
  }
}
```

可将上述任一对象保存为 `atmosphere.json`，仅同步大气而不要求磁图或 `path.dat`：

```bash
/opt/anaconda3/bin/python scripts/configure_case.py atmosphere.json --atmosphere-only
/opt/anaconda3/bin/python scripts/validate_workflow.py
```

该模式同步 `01_pfss`、`01_pfss_global_ct`、`02_relax`、`03_insert`、`04_fullmhd`
及所有 `04_fullmhd/profiles/*.par`，不改变网格、磁场、暗条加载或热物理开关，不复制
边界或路径文件，也不改写 `case.resolved.json`；应保留本次输入 JSON 作为配置记录。
正常完整配置也会同步上述大气参数。旧 `case.json` 若没有 `atmosphere`，将采用新的
默认密度归一化，不会继续保留旧压力归一化。

#### 3.3.3 旧算例迁移与验证边界

1. 本次更改了各阶段 `mod_usr.t`，必须重新编译需要运行的阶段。完成这次编译后，
   在已有接口内仅改 `&atmosphere_list` 参数无需再次编译。
2. 不要覆盖已完成算例的参数或 DAT。应在独立案例目录保存新旧配置，并从匹配的新
   大气重新生成初态；后续阶段必须使用同一套大气参数和底边界。
3. **改变大气参数不会重置 restart DAT 内部的密度、压力或温度。** 启动时仍会重建
   静力表并设置底边界，因此只替换续算程序／参数、沿用旧 DAT，可能造成“内部旧
   大气＋底边界新大气”的不一致；不能把这种续算称为已完成大气切换。
4. 使用已松弛 DAT 时，应匹配其原始大气归一化和底边界，而不是强行套新默认值。
   辐射、加热与热传导仍需独立检查；大气模型选择不会自动切换热物理 profile。
5. 在启动日志核对模型、归一化、参考高度及实际锚点，并检查首帧 DAT 的温度、
   密度、压力和底边界连续性，再做短时 MHD 测试。配置／一维构造检查不能替代
   完整 AMRVAC 编译和演化验证；本次更新不声明已完成新的全 MHD 生产验证。

#### 3.3.4 本地回归检查

在工作流根目录运行：

```bash
/opt/anaconda3/bin/python -B -m unittest discover -s scripts -p 'test_atmosphere*.py' -v
/opt/anaconda3/bin/python -B scripts/validate_workflow.py
```

配置测试检查五个阶段和两份备用配置的一致性、非法输入不写文件，以及其他物理参数
不被改变。数值测试需要 `gfortran` 和 NumPy：它直接提取五份 `mod_usr.t` 中的大气
例程并编译，而不是仅测试 Python 重写公式；没有编译器时会明确跳过数值测试。
测试覆盖默认解析大气、等温初态、压力锚点的空间插值、文件顺序覆盖和非法参数。

2026-10-22 的更新还在相同单位、64 个基级径向单元和 L4 条件下，与成功 RBSL15
Stage90 保存的 `inithdstatic()` 原例程进行了独立编译比较：8000 点高度、密度、压力
及两层底边界密度／压力共 24004 个值逐位一致。此结果只证明默认一维大气的回归，
不表示工作流的磁场、网格、热物理开关或完整三维演化已经与 RBSL15 相同。

### 3.4 场向热传导和 TRAC

热传导张量为

```math
\boldsymbol\kappa=\kappa_\parallel\mathbf b\mathbf b
+\kappa_\perp(\mathbf I-\mathbf b\mathbf b),\qquad \mathbf b=\frac{\mathbf B}{|\mathbf B|}.
```

默认关闭垂直热传导，使用 AMRVAC 的 Spitzer 系数

```math
\kappa_\parallel=8\times10^{-7}T^{5/2}
\quad {\rm erg\ cm^{-1}\ s^{-1}\ K^{-1}}.
```

`tc_saturate=.true.` 使用 AMRVAC 的平均自由程饱和修正，避免经典热流在陡峭温度
梯度处无限增大。`mhd_trac_type=2` 在未分辨过渡区内提高有效热传导并补偿热源/热汇，
以近似保持积分能量交换。TRAC 不是额外加热项，也不能修复不合理的边界鬼点温度。

### 3.5 光学薄辐射和背景加热

可选辐射项为

```math
Q_{\rm rad}=-n_en_H\Lambda(T)f_{\rm damp}(h).
```

#### 3.5.1 \(\Lambda(T)\) 本身：`JCcorona`

\(\Lambda(T)\) 是单位电子—氢粒子数密度乘积对应的光学薄辐射损失系数，cgs
单位为

```math
[\Lambda]={\rm erg\ cm^3\ s^{-1}}.
```

它只描述给定温度下等离子体把热能辐射出去的效率，不是单位体积的实际能量损失率。
真正进入内能方程的是

```math
\left|Q_{\rm rad}\right|
=n_en_H\Lambda(T)f_{\rm damp}(h)
\quad [{\rm erg\ cm^{-3}\ s^{-1}}].
```

因此必须区分“局部 \(\Lambda(T)\) 小”和“总冷却效应小”。当前完全电离、He/H=0.1
的单位关系给出 \(n_e\simeq1.2n_H\)，所以在温度相同且无高度阻尼时
\(\left|Q_{\rm rad}\right|\propto n_H^2\)。密度增加 \(f\) 倍会使单位体积辐射损失
增加约 \(f^2\) 倍；固定温度下，辐射冷却时间近似为

```math
t_{\rm rad}\sim
\frac{p/(\gamma-1)}{n_en_H\Lambda(T)}
\propto \frac{1}{n_H\Lambda(T)}.
```

所以高密物质即使位于较小的 \(\Lambda\) 区间，也可能具有不可忽略的体积积分辐射
损失。反过来，如果 \(\Lambda\) 在某个低温区突然变小，已经从高温侧冷却下来的物质会
在那里继续降温变慢，形成温度和质量堆积；这不等于该小 \(\Lambda\) 本身产生了强冷却。

AMRVAC 中的 `JCcorona` 是 Colgan--Feldman 2008 日冕光学薄损失曲线。它隐含生成
原始表时采用的元素丰度和电离平衡假设；当前 MHD 演化不会实时求解非平衡电离并重新
计算该曲线。源码 `src/physics/mod_radloss_tables.t` 保存了 45 个锚点：

```text
t_JCcorona = log10(T / K)
l_JCcorona = log10(Lambda / [erg cm^3 s^-1])
```

原始锚点从 \(\log_{10}T=4.0\) 延伸到 \(7.95351\)，即约
\(10^4--8.98\times10^7\ {\rm K}\)。几个代表点为：

| \(T\) [K] | \(\log_{10}\Lambda\) | \(\Lambda\) [erg cm³ s⁻¹] | 含义 |
| ---: | ---: | ---: | --- |
| \(1.0000\times10^4\) | -200.18883 | \(6.47\times10^{-201}\) | 表底部的数值关断，不应视为实测物理值 |
| \(1.3877\times10^4\) | -100.78630 | \(1.64\times10^{-101}\) | 仍近似无光学薄损失 |
| \(1.6594\times10^4\) | -30.60384 | \(2.49\times10^{-31}\) | 低温陡峭肩部的下侧 |
| \(1.9843\times10^4\) | -22.68481 | \(2.07\times10^{-23}\) | 在很窄温区内快速上升 |
| \(2.3729\times10^4\) | -21.76445 | \(1.72\times10^{-22}\) | 已进入显著辐射区 |
| \(9.9479\times10^4\) | -21.12642 | \(7.47\times10^{-22}\) | 约 \(10^5\) K 附近的高损失区 |
| \(1.1604\times10^6\) | -21.34104 | \(4.56\times10^{-22}\) | 典型日冕温度 |
| \(8.9848\times10^7\) | -22.51490 | \(3.06\times10^{-23}\) | 原始表的高温端 |

其中约 \(1.7\times10^4--2.4\times10^4\ {\rm K}\) 不是一个平缓的“小
\(\Lambda\)”平台，而是跨越多个数量级的极陡关断/开启肩部。该区间的低温一侧
\(\Lambda\) 绝对值很小，但 \(d\log\Lambda/d\log T\) 很大，因此单元跨过少量温度
就会得到完全不同的辐射损失率。这会增强温度演化对热传导、压缩、平流和离散温度
误差的敏感性。

初始化时，`mod_radiative_cooling.t` 不直接用 45 点表逐步查值，而是按以下顺序处理：

1. 在 \(\log_{10}T\) 上把原始表重采样为 `ncool` 个点，当前为 4000 点；
2. 通常采用二阶 Lagrange 插值；在表边缘或相邻斜率出现大跳变时退回线性插值，
   避免二阶插值产生明显过冲；
3. 把 \(\log_{10}T,\log_{10}\Lambda\) 转回线性量并变换为代码单位；
4. 演化时先按 \(\log_{10}T\) 定位相邻内部表格点，再在线性 \(T\) 和
   \(\Lambda\) 上作快速线性插值；
5. 高于原始表最高温度时，按轫致辐射近似
   \(\Lambda\propto\sqrt{T}\) 延拓。

`ncool` 只控制内部表的温度分辨率，不会改变原始 `JCcorona` 的物理形状。rbsl14
workflow27 的复现诊断表明，在 \(1.7\times10^4--2.4\times10^4\ {\rm K}\) 内，
上述正式插值和简单 log--log 线性插值的最大差异约为 `0.0029 dex`，没有发现能解释
异常冷却的插值过冲。因此当前问题应区分为：

- 早期从数十万 K 降温时的实际辐射能量损失；
- 密度平方乘子对单位体积及体积积分损失的放大；
- 进入低温陡峭肩部后的温度/质量停滞；
- TRAC、场向热传导、压缩功和平流对后期能量平衡的补偿。

不能仅依据低温端某一个 \(\Lambda\) 数值很小，就判断辐射冷却对整体演化不重要。

`JCcorona` 表自身最低温度约为 `10000 K`；这与用户指定的冷却适用下限和可选的
硬温度地板不是同一个概念。当前扩展源码把三个温度概念明确拆开：

| 名称 | 含义 | 是否改变原本更冷的单元 |
| --- | --- | --- |
| `tcoolmin` | 冷却表自身的最低温度，由所选 `coolcurve` 决定，不是 namelist 参数 | 否；表外不计算该曲线 |
| `tcool_cutoff` | 光学薄冷却的用户适用下限；实际值取 `max(tcool_cutoff,tcoolmin)` | 否；低于它的单元不施加光学薄冷却 |
| `tlow` | 仅供可选硬温度地板使用的目标温度 | 只有 `Tfix=.true.` 时才会补加内能 |

温度参数采用代码单位：

```math
T_{\rm physical}=T_{\rm code}\,T_{\rm unit}.
```

本工作流 `unit_temperature=10^6 K`，所以 `tcool_cutoff=0.02` 表示
`20000 K`。推荐的 20000 K 冷却截断测试为

```fortran
&rc_list
  coolcurve       = 'JCcorona'
  ncool           = 4000
  tcool_cutoff    = 0.02d0
  Tfix            = .false.
  tlow            = 0.008d0
  rc_split        = .true.
/
```

它表示：初始温度低于 `20000 K` 的光球或低色球单元不会被光学薄冷却项修改，
也不会被程序加热到 `20000 K`；从高温侧冷却下来的单元最多由该冷却项降到
`20000 K`。此时 `tlow=0.008` 不生效，因为 `Tfix=.false.`。

`Tfix=.true.` 的含义完全不同：每次辐射源更新后，程序检查全域温度，对
`T<tlow` 的单元直接增加内能，使其回到 `tlow`。因此若设置
`Tfix=.true., tlow=0.02`，就会把真实约 `6000 K` 的光球也抬到至少
`20000 K`，不能用来表达“冷却曲线只适用于 20000 K 以上”。

其余辐射参数含义如下：

- `coolcurve`：选择 \(\Lambda(T)\) 表；`JCcorona` 是当前日冕光学薄损失曲线。
- `ncool`：把原始冷却表插值成内部表时的温度采样数；增加它提高表分辨率，也增加
  每个 MPI rank 的初始化开销。
- `rc_split`：是否把辐射作为算子分裂源项更新；它不改变冷却曲线本身。
- `rad_newton`：额外的光学厚 Newton 冷却近似，与这里讨论的光学薄
  `tcool_cutoff` 分开；当前关闭。
- `mhd_trac` 和 TRAC 的 `Tcoff`：用于未分辨过渡区的导热/热源重分配，
  `Tcoff` 不是 `tcool_cutoff`。

2026-10-22 已从所有 `.par` 删除 `rad_damp`、`rad_damp_height`、
`rad_damp_scale`，不再规定低层高度衰减。共享依赖
`/share/home/guoyang1/codes/amrvac/src/physics/mod_radiative_cooling.t`
的默认值为 `rad_damp=false`；更换底层版本或叠加其他参数文件时要核对默认值与
覆盖值，不能仅凭本模板未写这一项就推断所有外部版本均关闭。

版本核验说明：上述共享源码实际 `cool_exact` 推进含球坐标径向衰减分支，
而 `getvar_cooling` 诊断函数没有对应分支。此前“球坐标下不生效”的初步判断
不适用于实际推进。现在通过移除三项并使用关闭的默认值取消该衰减，不是因为
其一定不受支持；辐射总开关、JCcorona、20000 K 截止和 Tfix 等设置没有删除。

光学薄辐射只移除本单元的能量，不会像辐射 MHD 那样吸收、传播或加热周围单元。

背景加热为

```math
H(r)=H_0\exp[-(r-r_{\min})/(60\ {\rm Mm})].
```

无辐射时若仍取 `H0=1e-4 erg cm^-3 s^-1`，体积分内能会持续增长。默认
`04_fullmhd/amrvac.par` 因此使用“无辐射、无背景加热”。无辐射、有背景加热的
热物理对照可使用 `profiles/amrvac_norad_with_heating.par`；完整热物理则使用
`profiles/amrvac_full_thermo.par`，并应先做短时稳定性测试。两份 profile 都保留
`itfixgrid=0`，切换热物理 profile 不会自动开启动态 AMR，也不保证复现某个历史算例。

### 3.6 PDFI/PFSS 方程

势场满足

```math
\nabla\times\mathbf B=0,\qquad
\nabla\cdot\mathbf B=0,\qquad
\mathbf B=-\nabla\Phi,\qquad
\nabla^2\Phi=0.
```

底边界给定观测 `Br`，外径向边界使用 source-surface 约束。当前局部楔形区域设置
`bcm=bcn=3`，即角向侧壁无穿透法向磁通。它与全日面球谐 PFSS 的物理方程相同，
但局部侧墙代替了区域外磁场的影响，所以高空零点、分界面和连通性可以不同。

### 3.7 RBSL 磁绳

设磁绳轴线为闭合曲线 `R(l)`，小半径为 `a`，

```math
\mathbf r=\frac{\mathbf x-\mathbf R(l)}{a},\qquad r=|\mathbf r|.
```

RBSL 矢势写成

```math
\mathbf A_{\rm RBSL}=\mathbf A_I+\mathbf A_F,
```

```math
\mathbf A_I=\frac{I}{a}\oint K_I(r)\,d\mathbf R,\qquad
\mathbf A_F=\frac{F}{4\pi a^2}\oint K_F(r)
\left(d\mathbf R\times\mathbf r\right).
```

代码采用的正则核为

```math
K_I(r)=
\begin{cases}
\frac{2}{\pi}\left[\frac{\arcsin r}{r}
+\frac{5-2r^2}{3}\sqrt{1-r^2}\right], & r<1,\\
r^{-1}, & r\ge1,
\end{cases}
```

```math
K_F(r)=
\begin{cases}
\frac{2}{\pi r^2}\left(\frac{\arcsin r}{r}-\sqrt{1-r^2}\right)
+\frac{2}{\pi}\sqrt{1-r^2}\\
\quad+\frac{5-2r^2}{2\sqrt6}
\left[1-\frac{2}{\pi}\arcsin\left(\frac{1+2r^2}{5-2r^2}\right)\right],
& r<1,\\
r^{-3}, & r\ge1.
\end{cases}
```

平衡关系使用

```math
I=\pm\frac{5\sqrt2}{12\pi}\frac{F}{a},
```

正负号由 `I_Helix` 决定。增大 `a` 而保持 `F` 不变会降低电流和平均轴向场，通常
减小扭转和初始洛伦兹力；增大 `F` 则同时增强电流和磁压。给定观测轴线曲率后，
仍应把由外部垂直场估计的平衡电流与上述 `I(F,a)` 沿轴线比较，不能只用单点值。

## 4. 数值方法

默认设置为：

| 项目 | 设置 |
| --- | --- |
| 坐标 | 三维球坐标 `(r, theta, phi)` |
| 时间推进 | 三步显式时间推进 |
| Riemann 求解 | HLL |
| 限制器 | van Leer |
| 磁场散度 | 交错网格 CT |
| 半相对论 | `c_red=6000 km/s` |
| 显式磁扩散 | `eta=200 km2/s` |
| AMR | Lohner 误差，密度权重 1，三步重网格 |
| 最大层级 | 4 |
| block | `12 x 12 x 12` |
| Courant 数 | 0.8 |
| 小量下限 | `rho,p = 1e-3` 代码单位 |
| 热传导 | 场向、饱和、STS `dtpar=0.2` |

`refine_threshold=0.2,0.4,0.6,0.8` 分别控制各层阈值；Fortran 表达式
`derefine_ratio=20*1.d0` 表示 20 个数组元素均为 1，并不是退化比等于 20。
`specialthreshold()` 在四个角向侧边界和顶边界附近增加容差，降低无科学目标区域的
过度加密。底边界块由 `special_refine_grid()` 强制保持最高层。

## 5. 边界条件

### 5.1 底边界 `r_min`

- 所有速度在鬼点固定为零，满足 line-tied 非流动边界。
- `rho` 和 `p` 固定为初始静力大气表。
- 切向面心磁场使用适用于拉伸网格的一侧等梯度外推。
- 法向面心磁场通过离散 `div(B)=0` 约束反算。
- 边界切向电场为零，因此底边界 `Br` 保持不变。

### 5.2 顶边界 `r_max`

- `r^2 rho` 连续；压强按顶层温度做等温外推。
- `v_r` 零梯度但禁止流入，`r v_theta` 和 `r v_phi` 连续。
- 切向磁场做几何/拉伸加权外推，法向分量由离散无散条件补齐。

### 5.3 四个角向边界

- 当前成熟代码在所有速度鬼点取零。
- `rho` 和 `p` 使用边界内侧值的镜像等值外推，不分别做 equal-gradient，避免
  `p/rho` 在鬼点产生十倍以上温度过冲。
- 切向磁场使用零梯度外推，法向面心磁场由离散无散条件反算。

边界条件修改后必须同时检查原始变量与守恒变量转换范围。只看 `saveprim` 的 VTU
不足以判断鬼点是否正确，因为 VTU 通常不保存鬼点。

## 6. 输入数据准备

### 6.1 Python 环境

```bash
cd Spherical_RBSL_FullMHD
python -m venv .venv
source .venv/bin/activate
python -m pip install -r Preprocess/requirements.txt
```

本机开发环境可直接使用 `/opt/anaconda3/bin/python`。

### 6.2 生成边界和路径

沿用 `Preprocess/prepare_rbsl_boundary.py` 与 `prepare_rbsl_path.py`。它们生成：

```text
cutout.dat
pfss_boundary.dat
parameters.json
path.dat
```

`path.dat` 在所有阶段统一为无文件头的小端 `float64` 二进制流，每个轴线点连续存储
`(longitude_deg+180, latitude_deg, radius_Rsun)` 三个数。默认 6000 点对应
`6000 * 3 * 8 = 144000` 字节。Fortran 使用 unformatted stream 直接读取，不再
支持旧的三列 ASCII `path.dat`。

再把局地无鬼点 `Br` 转成 PDFI 输入：

```bash
python Preprocess/make_pdfi_brll.py /path/to/rbsl_boundary_product
```

脚本会同时验证 `cutout.dat` 的物理层 `Br` 与插入全局 `pfss_boundary.dat` 的 patch
逐点完全相同。输出 `Brll0000.fits` 的 Python 数组顺序为 `(phi,theta)`；CFITSIO
在 Fortran 中读取为 `Brll(theta,phi)`。不要手工转置或额外加入鬼点。

### 6.3 旧模板原位配置接口（新用户见QUICKSTART）

```bash
cp case.example.json case.json
```

编辑 `boundary_dir`、`path_dat`、网格范围、RBSL 半径/磁通/螺度、暗条加载方式及
`atmosphere`（见 §3.3.2），
然后运行：

```bash
python scripts/configure_case.py case.json
python scripts/validate_workflow.py
```

`configure_case.py` 会：

1. 复制 `Brll0000.fits` 到 `01_pfss/Bphot/`；
2. 验证二进制 `path.dat` 大小，并把同一文件放到 PFSS 扣通量和 `03_insert`
   初始化插绳/质量源阶段；
3. 同步四阶段和两个 full-MHD profile 的网格范围、层级及 RBSL 参数，并同步
   `03_insert` 的暗条加载和演化时间；大气参数另同步到包含可选全球 PFSS 的全部
   五阶段及所有 full-MHD profile；全球 PFSS 的磁图、网格及 patch 参数仍按其
   README 单独配置；
4. 根据 `filament.method` 只打开 `instant`、`liakh_source` 或 `none` 中的一种模式；
5. 写出 `case.resolved.json` 作为可追溯记录。

只更新大气时使用 `--atmosphere-only`（见 §3.3.2）；该模式不会执行上述输入复制、
其他事件参数同步或 `case.resolved.json` 写入。未提供 `atmosphere` 对象会使用新的
默认解析大气及 10 Mm 数密度归一化。

## 7. 旧集群手工编译参考

服务器依赖：MPI-AMRVAC 3.3、Intel oneAPI Fortran、PDFI_SS、FISHPACK、CFITSIO。

```bash
export AMRVAC_DIR=/share/home/guoyang1/codes/amrvac
export GLOBAL_FMHD_LIBDIR=/share/home/guoyang1/apps/global_fmhd_deps/chenpf_libs
bash scripts/build_all.sh
```

若第 1 阶段使用全球球谐 PFSS-CT：

```bash
PFSS_STAGE=01_pfss_global_ct bash scripts/build_all.sh
```

必须使用 `make ARCH=intel`。默认 `ARCH=default` 曾导致 Intel `mpif90` 与 GNU
自由格式参数混用，使 `.f` 被按 fixed-form 解析并出现大量伪编译错误。

## 8. 四阶段运行方法

### 8.1 01_pfss

该阶段不推进时间，读取 `Bphot/Brll0000.fits`，计算并扣除 `Br_RBSL`，然后用
PDFI_SS 生成背景势场的边矢势和 CT 面心磁场。

```bash
cd 01_pfss
bsub < job
```

必须检查：输入 `Br` 范围、`Br_RBSL` 范围、扣除后范围、`mflux`、DAT/VTU 是否
有限，以及底边界 `Br_total` 在未来插绳后能否回到观测值。

### 8.2 02_relax

默认从 `pfss_0000.dat` 开始，不含磁绳，使用 CT、半相对论、重力、场向热传导和
TRAC；关闭辐射和背景加热。默认运行到 `t=1`，每 `0.05` 保存一次。

```bash
ln -sfn ../01_pfss/data/pfss_0000.dat restart.dat
bsub < job
```

不能只以最大动量进入平台判断收敛。至少要求连续两个输出间隔内：

- 体积分质量、内能、磁能和 RMS 速度变化均小于约 1%；
- 净边界质量/能量通量接近稳定；
- 最大速度、最大密度和最小时间步不再单调恶化；
- 温度/密度切片没有持续向外扩展的边界伪结构。

### 8.3 03_insert

把选定的无磁绳松弛 DAT 链接为 `restart.dat`。`firstprocess=.true.`、`iprob=2` 和
`do_rbsl=.true.` 使 `initonegrid_usr()` 在启动时只插入一次磁绳，并把新增磁场写回
CT 面场。默认采用 Liakh 时间源，运行到物理时间 110 s，每 10 s 保存一个 DAT，
不自动生成 VTU：

```bash
cd ../03_insert
ln -sfn ../02_relax/data/relax_0020.dat restart.dat
bsub < job
```

默认 `case.example.json` 使用 Liakh 时间源。在 `case.json` 中设置：

```json
"filament": {
  "method": "liakh_source",
  "density_factor": 20.0,
  "core_fraction": 0.6,
  "minimum_height_mm": 5.0,
  "liakh_mass_loading": {
    "density_contrast": 30.0,
    "start_time_s": 0.0,
    "load_time_s": 100.0,
    "sigma_horizontal_mm": 3.6,
    "sigma_radial_mm": 4.8,
    "minimum_height_mm": 5.0
  }
}
```

`method` 可取 `instant`、`liakh_source` 或 `none`。脚本会同步
`03_insert/amrvac.par` 并复制 `path.dat`。`insert_stage.time_max_code` 和
`save_interval_code` 控制该阶段演化到何时交给 `04_fullmhd`；`time_max` 必须覆盖完整
加载窗，且要根据 `reset_time` 判断 `liakh_t_start_s` 对应的是重置后的时钟还是
restart 原时钟。若只需生成插绳初态，可把 `time_max_code` 设为 0。

输出后仍必须做 `0000 -> 0001` CT 连续性检查。`03_insert` 自己从已有输出 continuation
时必须把 `firstprocess=.false.`、`do_rbsl=.false.`，否则会再次插入磁绳。

### 8.4 04_fullmhd

`04_fullmhd` 只做 continuation，固定使用 `firstprocess=.false.`、`do_rbsl=.false.`。
需要磁绳时接收 `03_insert` 选定时刻的 DAT；不需要磁绳时可直接接收 `02_relax`：

```bash
cd ../04_fullmhd
# 有磁绳
ln -sfn ../03_insert/data/rbsl_insert_load_0011.dat restart.dat
# 或无磁绳
ln -sfn ../02_relax/data/relax_0020.dat restart.dat
bsub < job
```

无磁绳路径中的 `02_relax` 必须来自 `01_pfss do_rbsl=.false.` 的未扣通量 PFSS。
若仍使用为插绳准备的 `Br_background=Br_observed-Br_RBSL`，却跳过 `03_insert`，底边界
就不会由磁绳补回到观测磁通，不能作为严格的无磁绳对照。

切换配置时先备份当前 `amrvac.par`：

```bash
cp profiles/amrvac_full_thermo.par amrvac.par
# 或
cp profiles/amrvac_norad_with_heating.par amrvac.par
```

切换 profile 后重新运行 `configure_case.py`，确保事件参数同步。

### 8.5 旧固定帧号链接和LSF依赖链（兼容入口）

默认输出编号固定时可执行：

```bash
bash scripts/link_restarts.sh
bash scripts/submit_chain.sh
```

默认 `FULLMHD_SOURCE=insert`，依赖链为 `01 -> 02 -> 03 -> 04`。无磁绳模拟改为：

```bash
FULLMHD_SOURCE=relax bash scripts/link_restarts.sh
FULLMHD_SOURCE=relax bash scripts/submit_chain.sh
```

此时不提交 `03_insert`，`04_fullmhd` 直接依赖并读取 `02_relax`。
提交前还必须确认该 `02_relax` 使用的是未扣除 RBSL 法向通量的 PFSS。

全球 PFSS-CT 依赖链使用：

```bash
PFSS_STAGE=01_pfss_global_ct bash scripts/link_restarts.sh
PFSS_STAGE=01_pfss_global_ct bash scripts/submit_chain.sh
```

首次新事件不推荐直接提交整条链。应逐阶段检查磁图、势场、松弛收敛和 CT 一步测试，
确认后再使用依赖链。

## 9. 热物理 profile

| 文件 | 热传导/TRAC | 辐射 | 背景加热 | 用途 |
| --- | --- | --- | --- | --- |
| `04_fullmhd/amrvac.par` | 开 | 关 | 关 | 推荐的轻量动力学基线 |
| `profiles/amrvac_norad_with_heating.par` | 开 | 关 | `1e-4` | 无辐射、有加热的热物理对照；仍默认冻结树，可能积热 |
| `profiles/amrvac_full_thermo.par` | 开 | 开 | `1e-4` | 完整热物理短测试后再生产 |

辐射不是控制焦耳加热产生的机制。它只在能量方程中移除热能，不能减少磁重联产生的
`eta J^2`；在高密区因为损失近似正比于 `n_e n_H`，反而可能触发快速降温、压力下降、
压缩和进一步增密。背景加热与辐射应作为成对的能量收支进行标定。

## 10. 关键测试得到的知识

### 10.1 PFSS 与磁场初始化

- PDFI 局部势场和旧全日面球谐 PFSS 在低层强场主体上可高度相关，但高空弱场、
  零点和分界面会因局部侧墙、磁图分辨率和径向离散显著不同。
- 在一组干净对照中，全域矢量相关系数约 0.988、相对 L2 差约 19%；100 G 以上
  区域的方向差中位数约 4 度。高于 100 Mm 后相对差迅速增大。
- PDFI 径向壳层由 288 增到 576 并配合 `log(r)` 三次重构时，曾出现高空场衰减
  明显变慢。网格更细不自动代表解更正确，必须做固定边界的径向收敛测试。
- 把 PFSS 面心磁场直接保存为后续 `cutout_pfss.dat` 能减小 restart 边界梯度和强场
  区域的非物理爆发，但不能单独消除所有离散洛伦兹力。

### 10.2 强场、拉伸网格和 AMR

- 底层 2000 到 3000 G 强场对离散 `J x B` 极敏感。即使连续 PFSS 是无流势场，
  映射、限制、边界填充和粗细网格交换仍会产生数值电流。
- `J_num` 不能只用 `Delta B/Delta r` 的直觉判断。网格变细时，如果不同层重构并非
  对同一光滑场高阶一致，`Delta B` 不一定按 `Delta r` 同阶下降，故 `J_num` 可能增大。
- 早期崩溃点多次落在 AMR 粗细交界。固定强磁场区为 L4 能显著缓解，但会增大 DAT
  和 VTU。动态 AMR 必须监视每次 regrid 前后的 `dt`、块数和极值位置。
- 半相对论只限制阿尔芬特征速度，不能阻止未平衡洛伦兹力把等离子体加速到非物理
  实际速度。降低 Courant 数通常只是延迟失稳。

### 10.3 边界和热物理

- 四侧边界分别对 `rho`、`p` 做 equal-gradient 曾使鬼点 `p/rho` 产生 10 到 20 倍
  温度过冲，热传导随后把异常高温带入活动区并压低时间步。当前改为 `rho,p` 等值
  外推。
- 打开辐射的早期测试出现低温/高温交替条带；关闭热传导、辐射和 TRAC 后条带消失。
  分项测试表明辐射是主要放大源，热传导负责沿磁场传播和重分配。
- 无辐射的 TC+TRAC 松弛仍形成低温高密舌状结构。与无热传导组比较，舌状掩膜
  Jaccard 约 0.996，说明它不是热传导产生的，而是 PFSS 初态离散力不平衡驱动低熵
  物质沿磁场上升。
- 舌状区径向速度多数向上，速度与三维磁场夹角中位数约 10 度；磁场引导输运，
  但二维舌头轮廓不等同于单条磁力线。

### 10.4 电阻率和能量预算

- `eta=200 km2/s` 与早期 `eta约8100 km2/s` 对照中，`t=0.05` 的累计焦耳热只约为
  高电阻组的 14%，但不是简单的 1/40，因为低电阻保留了更大的电流密度。
- 低电阻组早期压缩功 `integral(-p div v)dV` 由负转正，局部密度和压缩增强。因此
  内能仍可能比高电阻组高，不能只根据 `eta` 判断温度趋势。
- 推荐输出并积分：内能、磁能、动能、`integral eta J^2 dV`、
  `integral -p div(v)dV`，同时记录边界能流，才能区分焦耳加热、压缩加热和输运。

### 10.5 收敛和输出判断

- “动量平台”不等于热力学收敛。曾有区域动能下降约 1%，但同一时段质量增加
  4.8%、内能增加 11.4%、磁能增加 1.3%。
- DAT 日志中的最后一个 `dt=0` 可能只是精确命中 `time_max`，不是崩溃。
- `saveprim` 只决定输出变量表示，不改变内部守恒变量推进，也不包含边界鬼点。
- 每次插绳后至少比较 `0000/0001` 的磁场、温度、密度和磁力线投影，不能只看
  单独的 `0000.vtu`。

## 11. 运行监控与验收

每个长作业至少记录：

```text
it, time, dt
L1/L2/L3/L4 active blocks
max/min rho, p, T, |B|, |v|
min local MHD timestep and its (r,theta,phi,level)
cs, vA, cf, beta, J/B at the limiting cell
volume-integrated internal, magnetic, kinetic energy
integral eta J^2 dV
integral -p div(v) dV
```

验收条件：

1. `stderr` 为空，无 NaN/Infinity；
2. `dt` 没有跨多个数量级持续下降；
3. regrid 前后无极值瞬间跳变；
4. `0000 -> 0001` 磁绳保持；
5. 底边界 `Br` 与目标磁图的差异可解释；
6. 温度下限占比、最大密度和总质量没有非物理单调增长；
7. 能量变化能由焦耳、压缩、热传导、辐射、加热和边界通量闭合解释。

## 12. VTU 传输

从本机流式中转，不在本机保存临时 VTU：

```bash
scripts/upload_vtu_to_yangyang.sh \
  '/share/home/guoyang1/zhouchang/CASE/RUN/data/*vtu'
```

脚本会自动推断 `/home/guoyang/zhouchang/CASE/RUN`，先写 `.part`，比较源端和目标端
SHA256，校验后原子改名。也可以第二个参数显式指定 yangyang 目录。

## 13. 统一切片分析

[`scripts/spherical_slice_analysis.py`](scripts/spherical_slice_analysis.py)
统一完成“simesh 读取拉伸 DAT → 均匀切片/薄层缓存 → 派生量 → 绘图”。固定
`phi`、固定 `theta`、固定 `r` 是主要模式，同时支持由中心点、Cartesian 法向和
向上方向定义的自由角度平面。温度、密度、压力、速度、磁场、阿尔芬速度、
辐射损失和背景加热可用单面；`J`、`J/B`、欧姆耗散、`div(B)` 和 `div(v)`
应使用 `local-fit` 后端在各点所属 AMR 单元的真实尺度上做多尺度三维拟合。
原来的均匀薄层差分作为 `uniform-slab` 预览后端保留，但不再推荐用于定量图。
`plot-amr-grid` 子命令则不做均匀重采样，直接绘制与固定 `phi` 或固定 `theta`
平面相交的原生 AMR 叶块轮廓及块内单元线，并把全域/切面各层叶块统计写入 JSON；
可选参考点用于标记磁绳轴线位置。

完整命令、字段名、单位和物理口径见
[`scripts/README_spherical_slice_analysis.md`](scripts/README_spherical_slice_analysis.md)。
`local-fit` 会保存 AMR level、跨块/跨层 stencil、拟合残差、尺度敏感性和可靠性
掩膜。默认可靠性掩膜由物理边界、拟合残差和多尺度一致性决定，不会仅因跨 AMR
level 就留下块状空洞；另存的严格掩膜会排除跨层模板，供块边界敏感性检查。
即便如此，缓存中的 `divb_code` 仍不等价于 CT 面心磁场的 `get_divb()`；
TRAC 的动态截止温度也不存储在 DAT 中，因此后处理辐射损失是当前状态上的模型估计。

## 14. 运行前检查清单

1. `Brll0000.fits` 无鬼点，shape、经纬范围和 Fortran 顺序正确；
2. PFSS 和插绳使用完全相同的小端 `float64` 二进制 `path.dat`、`Bflux`、`Rmfr`、螺度和点数；
3. `01_pfss do_rbsl=.true.`，`02_relax firstprocess=.false.`；
4. `03_insert firstprocess=.true.`、`iprob=2`、`do_rbsl=.true.`，只在启动时插入一次
   RBSL；它自己的 continuation 必须把 `firstprocess/do_rbsl` 都改为假；
5. `04_fullmhd firstprocess=.false.`、`do_rbsl=.false.`，只读取 `02` 或 `03` 的 DAT；
6. 暗条瞬时增密与 Liakh 时间源不能同时开启，`03_insert/path.dat` 与 PFSS 扣通量
   使用的轴线逐字节相同；
7. Liakh 的 `03_insert time_max` 覆盖加载窗，`t_start/t_load` 秒数与 `reset_time` 口径一致；
8. 四阶段网格、拉伸参数、block 和最大 AMR 层级完全一致；
9. 编译使用 `ARCH=intel`，PDFI/FISHPACK/CFITSIO 动态库路径正确；
10. job 核数是 64 的倍数且包含 `span[ptile=64]`；
11. 新事件先逐阶段运行，不直接提交整条依赖链；
12. `03_insert` 启动后立即完成 `0000 -> 0001` CT 连续性检查和能量预算基线。

## 15. 参考模型

- Titov et al.：regularized Biot-Savart law（RBSL）磁通绳构造。
- Gombosi et al.：Boris 型半相对论 MHD。
- Xia et al.、Keppens et al.：MPI-AMRVAC 与 AMR MHD。
- Gardiner & Stone：交错网格 constrained transport。
- Spitzer：场向电子热传导。
- Liakh, Luna & Khomenko (2023), A&A 673, A154：有限时间暗条质量源。

本 README 描述的是当前归档代码真实启用的方程和开关。论文写作时应以实际运行
目录的 `amrvac.par`、启动日志和源码哈希为准，不应仅引用某个 profile 的计划设置。

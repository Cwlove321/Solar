# 笛卡尔双极场 RBSL 磁绳模拟

本目录是一套纯理论的三维笛卡尔 MPI-AMRVAC 模板：在解析双极磁场中插入一条
Regularized Biot–Savart Laws（RBSL）磁绳，先用磁摩擦方法松弛磁场，再从松弛
结果启动零 beta MHD 演化。

目录名采用 `Cartesian_RBSL_Bipolar`，因为模型的核心不是单独插入一个磁偶极
场，而是“解析双极背景场 + RBSL 磁绳”。它不读取 HMI、AIA 或外部边界 DAT，
初始磁场和磁绳几何全部在 `mod_usr.t` 中解析构造。

## 目录结构

```text
Cartesian_RBSL_Bipolar/
├── mf/
│   ├── mod_usr.t                 # 双极场、RBSL 磁绳和磁摩擦边界
│   ├── amrvac.par                # MF 网格、AMR、输出和物理参数
│   ├── data/                     # MF 输出目录
│   └── magnetic_geometry.*       # 当前解析几何示意图
├── 0beta/
│   ├── mod_usr.t                 # 从 MF 磁场启动零 beta MHD
│   ├── amrvac.par                # MHD restart、网格和保存参数
│   └── data/                     # MHD 输出目录
└── README.md
```

## 计算流程

```text
解析磁荷产生的 bipolar potential field
                     +
       RBSL 半圆形磁绳轴线与有限截面电流
                     │
                     ▼
                  mf/
             磁摩擦松弛到近似 NLFFF
                     │
          restart: mfNNNN.dat
                     │
                     ▼
                 0beta/
       填充密度并进行零 beta MHD 演化
```

必须先运行 `mf/`，再运行 `0beta/`。两个阶段的物理区域、基础网格、AMR 和
磁场归一化必须一致。

## 解析磁场和磁绳几何

`mf/mod_usr.t` 的标题是 `Magnetic flux rope in a bipolar field`。当前主双极场
由埋在下边界以下的一对正负磁荷构造，RBSL 轴线位于 `y-z` 平面，磁绳中心的
`x` 坐标取计算盒中心。

![解析双极场和 RBSL 磁绳几何](mf/magnetic_geometry.png)

当前主要参数位于 `mf/mod_usr.t::initglobaldata_usr()`：

| 参数 | 当前值 | 含义 |
| --- | ---: | --- |
| `unit_length` | `1e9 cm` | 一个长度代码单位为 10 Mm |
| `q_para` | `+0.75e20/(unit_B unit_L²)` | 主双极场磁荷强度尺度 |
| `d_para` | `1` | 主磁荷埋深，10 Mm |
| `L_para` | `1` | 主磁荷半间距，10 Mm |
| `q1_para` | `-1.0e20/(unit_B unit_L²)` | 第二组解析背景分量强度 |
| `d1_para` | `2.5` | 第二组磁荷埋深，25 Mm |
| `L1_para` | `2.5` | 第二组磁荷半间距，25 Mm |
| `r0` | `2` | 磁绳轴大半径，20 Mm |
| `a0` | `1` | 磁绳小半径，10 Mm |
| `k_nonb` | `0.9` | 控制电流相对平衡估计的倍率 |
| `np` | `200` | 地面以上半圆轴线的离散点数 |

按照 `unit_length=1e9 cm=10 Mm`，`r0=2` 实际对应 20 Mm；代码中的行内
注释已同步修正。修改参数时应始终先换算代码单位，并与几何示意图交叉检查。

轴线参数由以下关系给出：

```text
apex height = r0 - d_para
half footpoint separation = sqrt(r0² - d_para²)
ring center z0 = -d_para
```

电流 `I_cur` 根据主双极场在磁绳位置的垂直分量估计，再乘 `k_nonb`；轴向磁通
`F_flx` 随后由 `I_cur` 和 `a0` 构造。修改 `q_para`、`r0`、`a0` 或
`k_nonb` 会同时改变磁绳电流、磁通和初始非平衡程度，不能只改单个参数而不做
稳定性检查。

## 第一阶段：磁摩擦松弛

`mf/mod_usr.t` 使用 `mod_mf`，在 `Cartesian_3D` 坐标中：

1. 通过解析磁荷构造双极背景场的矢势；
2. 用 RBSL 核函数构造磁绳的电流与轴向磁通；
3. 在 staggered/CT 网格上由总矢势生成磁场；
4. 使用磁摩擦速度逐步降低洛伦兹力；
5. 输出 `divB`、电流和力自由度相关辅助量及分析记录。

当前 `mf/amrvac.par` 使用：

```text
domain_nx = 60 × 60 × 60
domain    = [-10,10] × [-10,10] × [0.3,20.3]
block_nx  = 10 × 10 × 10
refine_max_level = 3
typedivbfix = ct
time_max = 100
dtsave_dat = 25
base_filename = data/mf
```

按当前保存间隔，完整运行通常会产生 `mf0000.dat` 到 `mf0004.dat`，但正式进入
MHD 前仍应根据日志和力自由诊断选择真正完成松弛的最后一帧。

## 第二阶段：零 beta MHD

`0beta/mod_usr.t` 使用 `mod_mhd`。它从 MF 数据读取并保留磁场，重新填充密度，
将速度设为零，然后执行零 beta MHD 演化。当前 `amrvac.par` 中：

```text
restart_from_file = ../mf/data/mf0004.dat
base_filename     = data/eru
domain_nx         = 60 × 60 × 60
domain            = [-10,10] × [-10,10] × [0.3,20.3]
refine_max_level  = 4
mhd_energy        = false
mhd_adiab         = 0
typedivbfix        = ct
```

用户模块还定义了：

```text
t_relax  = 5
t_driven = 20
v_0      = 7.5 km/s（再除以速度单位）
```

底边界驱动窗口为 `t_relax ≤ t ≤ t_relax+t_driven`。修改驱动时长或速度前，应
同时检查底边界速度、电场约束、时间单位和输出时间分辨率。

## 编译和运行

需要 Guo/项目对应版本的 MPI-AMRVAC、MPI 和 Fortran 编译器。每个阶段必须在
自己的目录单独运行 `setup.pl` 和 `make`。

### 1. 编译并运行 MF

```bash
cd /path/to/Cartesian_RBSL_Bipolar/mf
export AMRVAC_DIR=/path/to/amrvac
$AMRVAC_DIR/setup.pl -d=3
make -j 8
mpirun -np 64 ./amrvac -i amrvac.par
```

集群上应使用对应调度系统提交 MPI 作业，不要在登录节点直接运行大规模任务。

### 2. 选择 MF restart

确认以下内容：

- 输出 DAT 完整且没有 MPI/NaN 错误；
- 磁摩擦分析量已经趋于稳定；
- `divB` 和归一化洛伦兹力满足科学要求；
- 选中快照的网格与 `0beta/amrvac.par` 完全一致。

如果最终文件不是 `mf0004.dat`，修改 `0beta/amrvac.par` 中的
`restart_from_file`。

### 3. 编译并运行零 beta MHD

```bash
cd ../0beta
$AMRVAC_DIR/setup.pl -d=3
make -j 8
mpirun -np 64 ./amrvac -i amrvac.par
```

## 修改新实验时的检查清单

1. `mf` 和 `0beta` 的 `domain_nx*`、`xprobmin/max*`、block 大小一致；
2. restart 文件路径从 `0beta/` 目录解析且文件真实存在；
3. `unit_length`、`unit_temperature`、`unit_numberdensity` 在两阶段一致；
4. 磁荷埋深、间距、磁绳大/小半径满足几何关系；
5. `r0 > d_para`，否则足点半间距 `sqrt(r0²-d_para²)` 无实数解；
6. `a0` 相对网格间距有足够分辨率，尤其在最低 AMR 层级；
7. 修改 `mod_usr.t` 后重新运行 `setup.pl` 和编译；
8. 比较初始场、松弛场和 MHD 初态，确认 restart 分量映射正确。

## 当前模板中的历史项

这是一套研究用理论实验模板，`amrvac.par` 中保留了少量位于 namelist 结束符
`/` 后面的旧参数候选值和重复行。它们通常不属于当前生效 namelist，但容易让人
误判。开展新实验时建议：

- 只在对应 `&... /` namelist 内保留唯一生效值；
- 删除无效的拼写，例如历史候选行中的 `.flase.`；
- 为每组物理参数使用独立输出目录；
- 不复用旧实验 DAT 作为新参数的 restart，除非网格和变量定义完全一致。

本目录不需要 Python 预处理或观测数据包；其外部依赖只有 MPI-AMRVAC 的编译和
运行环境。

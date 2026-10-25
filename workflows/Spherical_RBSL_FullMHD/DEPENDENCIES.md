# 接收者安装与依赖版本

本交付包含工作流及匹配参数接口的求解器源码。接收者先安装一次 Python
预处理环境、MPI/Fortran 编译环境和三个外部库，再按 `QUICKSTART.md` 生成案例。
现成的集群二进制、个人磁图、轴线、旧 DAT、口令和运行环境没有随包分发。

## 随包锁定的源码

| 依赖 | 路径与版本 | 用途 |
|---|---|---|
| MPI-AMRVAC | `vendor/amrvac`；上游 HEAD `be70e911d82861bee017e1879141d0565a64733b` 加导出时的本地修改 | 所有 Fortran 阶段 |
| simesh | `vendor/simesh`；来源锁定上游 `e223ceaa7c0c5d72bbcf24f0143b250ad8ed4a0d` 加项目修正 | 可选的平衡电流 DAT 采样分析 |

来源、许可证与修改文件清单见 `vendor/SOURCE_PROVENANCE.md`。实际交付版本由
`vendor/dependencies.lock.json` 中的逐文件 SHA-256 确定：

```bash
python scripts/verify_dependencies.py
```

该命令仅用 Python 标准库，离线运行，不编译或安装。它核对原始源码文件，
额外产生的编译缓存不纳入核对；任何已锁定源文件发生变化都会报告。
AMRVAC 快照包含 `src/`、`arch/`、`tools/`、文字文档、`setup.pl` 与 LICENSE，
没有依赖指向外部目录的符号链接。`lib/` 会在首次编译时创建。

这份 AMRVAC 保留 `tcool_cutoff`、半相对论 MHD、内能方程、TRAC 和 CT。
只写“AMRVAC 3.3”不足以复现这些扩展；不要将普通上游版本或另一个未核对
源码目录替换进来。导出来源是工作流当前配置的共享源码，**不是**历史
RBSL15 `code_production` 的特定 AMR-CT 修复/非有限数保护二进制归档。

## Python 预处理环境

要求 Python 3.11 或更高版本。以下 `python` 指接收者选定的解释器；
原开发机验证使用 `/opt/anaconda3/bin/python`，不要求接收者也采用这个路径。

```bash
python -m venv .venv
.venv/bin/python -m pip install -r requirements-preprocess.txt
```

之后用该解释器调用`workflow.py`，例如`.venv/bin/python workflow.py check CASE`。
预处理依赖包括 NumPy、SciPy、
Astropy、SunPy Map、Matplotlib、Plotly、pyevtk与OpenCV模板匹配。选用
`opencv-python-headless`提供计算功能，交互取点由Matplotlib图形后端承担。`sunpy[map]` 安装 Map 所需
扩展依赖。`tested-environment.txt` 记录本次检查环境的具体版本，不是将
开发机所有包/操作系统组件硬锁给接收者的安装清单。

计算平衡电流或自行制作 GIF 时再安装可选分析依赖：

```bash
.venv/bin/python -m pip install -r requirements-analysis.txt
.venv/bin/python -m pip install --no-build-isolation ./vendor/simesh
.venv/bin/python scripts/calculate_rbsl_balance_from_stretched_dat.py --help
```

simesh 需要 C 编译器与 Cython 编译本机扩展；默认不要求 OpenMP。求解器运行、
案例生成、DAT 交接验收不依赖 simesh。随包平衡电流工具直接使用本工作流内的
实现，不需要相邻的旧 `Spherical_RBSL` 工作流。读拉伸 DAT 时仍须传入
该 DAT 对应的 `--grid-config amrvac.par`。

## Fortran/MPI 与外部库

机器需要 Linux/Unix、GNU make、Perl、rsync、MPI Fortran 包装器（例如
`mpif90`）以及匹配的 Fortran 编译器。源码预处理器和构建规则已随 AMRVAC
提供。Intel 或 GNU 工具链应各自使用一致的 MPI 与库，不混用不同编译器
生成的 `.mod` 或静态库。用户的工作流已在 Intel 集群上运行过；本次未在
GNU 或接收者机器验证全流程数值运行。

当前五份 `mod_usr.t` 都包含 PDFI/FITS 相关例程，因此即使选择全球 PFSS 分支，
链接时仍需以下三个外部库：

| 库 | 用途与接收者操作 | 官方来源 |
|---|---|---|
| PDFI_SS | 提供 `scrbpot_ss`、`mflux_ss`、`ahpot_ss` 等球面势场例程。下载作者源码，按其 README 用当前 Fortran 工具链构建 `libpdfi_ss.a` 或对应共享库；同时保留版本和许可证 | [作者发布的源码快照](https://zenodo.org/records/3711571) |
| FISHPACK | PDFI_SS 的椭圆方程求解依赖。使用 PDFI_SS README 指定的传统接口和精度构建 `libfishpack`，不能盲换成改变了函数接口的现代重写版 | [NCAR 经典库说明](https://ncar-hpc-docs.readthedocs.io/en/latest/environment-and-software/ncar-classic-libraries-for-geophysics/) |
| CFITSIO | 提供 `ftopen` 等 Fortran FITS 接口。按官方文档配置本机前缀，执行 `make`、`make check`、`make install` | [HEASARC 源码与构建说明](https://github.com/HEASARC/cfitsio) |

例如 CFITSIO 可在自己有写权限的安装前缀下构建：

```bash
./configure --prefix=/absolute/path/to/deps
make
make check
make install
```

PDFI_SS/FISHPACK 的版本和构建选项应以作者发布的配套 README 为准，尤其确认
实数精度与外部符号接口。本文没有为未经实际编译的第三方包编造通用构建
参数；接收者可使用管理员已验证的同工具链安装。外部库源码/预编译库没有
自动下载或打包，接收者应按各自许可获取。

将库目录、MPI 启动器、AMRVAC 路径、架构与编译环境写入机器配置，见
`machine.example.json` 与 `QUICKSTART.md`。用随包 AMRVAC 时，
`amrvac_dir` 可指向本工作流的 `vendor/amrvac`。

本机编译的基本形态如下；正式执行由运行脚本根据机器配置生成：

```bash
export AMRVAC_DIR=/absolute/path/to/workflow/vendor/amrvac
cd /absolute/path/to/case/01_pfss
perl "$AMRVAC_DIR/setup.pl" -d=3 -arch=default
make ARCH=default "VACPP=perl $AMRVAC_DIR/src/vacpp.pl"
```

显式使用Perl调用源码预处理器，可避免ZIP解压后丢失可执行权限；统一入口会自动处理。

`arch/default.defs` 使用 GNU 风格选项，Intel 编译器应选匹配的架构定义。
`local.make` 的库目录来自机器配置。只有库由 Intel 编译且需要 Intel runtime
时才添加相应运行库；不能把开发机的 Intel 目录复制到别人的路径中。
共享库运行时需同时让加载器能找到安装目录（例如在机器环境中配置
`LD_LIBRARY_PATH`），仅有链接目录并不能保证运行时加载成功。

独立案例保存其机器配置。把案例搬到另一台机器时，应重新配置求解器、库和
解释器路径；不要继续使用旧机器上解析出的绝对路径。

## 本次验证边界

本次已核对全部随包源文件 SHA-256、求解器必要接口、DAT v5 写入格式与 CT
块长度，检查平衡电流工具可从任意工作目录加载自身实现，并验证磁通/电流
换算往返一致。Python 环境版本已记录。

本次没有安装额外 Python 包，没有在登陆节点编译或运行 MPI，也没有提交新
模拟。Cython 在本次开发机基础环境中未安装，因此未声称已构建可选 simesh
扩展。接收者第一次编译和首次新事件运行仍须由实际机器检查结果确认。

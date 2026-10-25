# 随包依赖源码来源

快照日期：2026-10-22。每个文件的 SHA-256 和完整快照标识见
`dependencies.lock.json`，在工作流根目录执行 `python scripts/verify_dependencies.py`
可离线核对。本文件记录来源，不替代各包的许可证。

## MPI-AMRVAC

- 上游：<https://github.com/amrvac/amrvac>。
- 上游 Git HEAD：`be70e911d82861bee017e1879141d0565a64733b`。
- 来源：guoguo 上实际使用的 `codes/amrvac` 工作目录；它**不是干净的上游提交**，
  应使用随包的全部实际源码，不能仅检出 HEAD 就认为等价。
- 纳入：`src/`、`arch/`、`tools/`、文字 `doc/`、`setup.pl`、原始 `LICENSE` 与 `README.md`。
  排除 `.git`、已编译 `lib`、测试输出、历史备份、演示媒体与缓存。
- 原始许可证：`amrvac/LICENSE`（GNU GPL v3）。保留原源码中的署名与许可头。
- 导出时已修改的 tracked 文件：

```text
src/io/mod_gimli.t
src/mhd/mod_mhd_phys.t
src/particle/mod_particle_advect.t
src/particle/mod_particle_base.t
src/particle/mod_particle_gca.t
src/particle/mod_particle_lorentz.t
src/particle/mod_particle_sample.t
src/physics/mod_pfss.t
src/physics/mod_radiative_cooling.t
```

本次交付没有更改这些求解器文件。关键扩展 `tcool_cutoff`、内能方程、
TRAC 与 CT 都随实际源码保留。`src/io/mod_input_output_helper.t` 明确写 DAT v5；
DAT 带 CC 演化变量及 CT staggered 数据。尚未在接收者机器上编译或执行新事件。

## simesh（可选分析依赖）

- 上游：<https://github.com/Astery0502/simesh>。
- 来源包锁定的上游提交：`e223ceaa7c0c5d72bbcf24f0143b250ad8ed4a0d`。
- 来自项目现有 `Packages/simesh` 源码，含项目的拉伸网格读取修正；
  随包文件哈希是本次交付的权威版本，不能只用 PyPI 同名版本代替。
- 原项目依赖锁中的完整包 tree 标识：
  `41b1b97327b7ca84b7bbfff84692a90e6616a2e2a2c1b81a05cab4e2d5770fcc`。
  当前交付只包含构建必需源码、构建脚本、安装说明和 LICENSE，
  所以裁剪后快照标识不同，见本目录 `dependencies.lock.json`。
- 原始许可证：`simesh/LICENSE`（GNU GPL v3）。
- 没有分发已有 `.so` 扩展、生成的 C 文件、缓存或个人运行数据。

## 平衡电流工具

`scripts/rbsl_balance_dat_core.py` 来自原项目
`workflows/Spherical_RBSL/scripts/rbsl_balance_dat_core.py`，现复制到当前工作流。
只将 simesh 定位改成当前工作流的 `vendor/simesh/src`；公式和采样算法未变。
原文件没有单独的许可证声明，本次不为原作者代码另行声称许可。

# 球坐标有限体积磁螺度

本工作流按照 Moraitis et al. (2018, *Solar Physics*, 293:92) 在有限球坐标楔形区域
中计算相对磁螺度。适用于三维球坐标 MPI-AMRVAC `.dat`，包括径向拉伸和 AMR 网格。

计算顺序为：

1. 使用原始 `amrvac.par` 恢复 DAT 没有保存的径向拉伸参数；
2. 用仓库内的 `simesh` 将 `(B_r,B_\theta,B_\phi)` 线性采样到均匀的
   `(r,\theta,\phi)` 单元中心网格；
3. 对逐分量插值产生的离散散度做有限体积 Hodge 投影；投影保持六个边界的法向
   磁通不变，并将每个均匀网格单元的净磁通降到线性求解器误差；
4. 在六个边界上令 `n·B_p=n·B`，用有限体积全 Neumann Laplace 方程求
   `B_p=grad(Phi)`；
5. 在 `A_r=0` 的 DeVore simple 规范下分别计算 `A` 和 `A_p`；
6. 使用精确球坐标单元体积计算

   ```text
   H_V  = integral (A + A_p) dot (B - B_p) dV
   H_J  = integral (A - A_p) dot (B - B_p) dV
   H_PJ = 2 integral A_p dot (B - B_p) dV
   H_V  = H_J + H_PJ
   ```

## 使用

必须给出与 DAT 同一次计算使用的 `amrvac.par`，否则程序拒绝运行，避免把拉伸径向
网格误当作均匀网格。

```bash
/opt/anaconda3/bin/python \
  workflows/Spherical_Helicity/spherical_helicity.py \
  /path/to/snapshot.dat \
  --grid-config /path/to/amrvac.par \
  --resolution 64,32,64 \
  --references bottom,top \
  --output /path/to/helicity_result
```

在 `guoguo` 的 `zcc` Python 环境中：

```bash
/share/home/guoyang1/apps/miniconda3/envs/zcc/bin/python \
  /share/home/guoyang1/zhouchang/solar_code/workflows/Spherical_Helicity/spherical_helicity.py \
  /share/home/guoyang1/zhouchang/20241009/rbsl15/03_insert_doubleieq_a10_vtu/data/r15_i2eq_0000.dat \
  --grid-config /share/home/guoyang1/zhouchang/20241009/rbsl15/03_insert_doubleieq_a10_vtu/amrvac.par \
  --resolution 64,32,64 \
  --output /share/home/guoyang1/zhouchang/20241009/rbsl15/helicity_test/r15_i2eq_0000
```

输出：

- `spherical_helicity.json`：输入、网格、边界磁通、求解器、散度、规范差异以及
  `H_V/H_J/H_PJ`；
- `spherical_helicity_fields.npz`：均匀网格上的 `B/B_p/A/A_p`、螺度被积函数和
  数值诊断场；不需要三维场时可添加 `--no-save-fields`。

## 单位

默认结果为 DAT 代码单位：磁螺度单位是 `B_unit^2 L_unit^4`，磁通单位是
`B_unit L_unit^2`。同时指定下面两个参数可额外输出 `Mx^2` 和 `Mx`：

```bash
--unit-length-cm 6.957e10 --unit-b-gauss VALUE
```

不要猜测 `unit_b_gauss`；应从该算例的 MPI-AMRVAC 归一化参数中读取。

## 必须检查的质量指标

- `boundary_flux.relative_imbalance`：原磁场六个边界的相对净磁通。默认超过
  `5e-3` 时停止；小的不平衡采用所有边界上最小的统一法向场偏移修正，并明确记录；
- `potential_solver.relative_linear_residual`：势场线性方程残差；
- `solenoidal_projection.relative_field_correction_l2`：均匀化后为了恢复无散条件所需
  的相对场修正；该值过大说明当前重采样分辨率不足；
- `solenoidal_projection.cleaned_flux_residual_l2`：投影后所有单元净磁通的 L2；
- `quality.status`：只有无散投影修正和最佳 DeVore 旋度误差均通过阈值时才是
  `validated`；否则明确标记为 `diagnostic_only`；
- `quality.cell_centered_differential_divergence_l2`：由单元中心数组直接差分得到的
  辅助指标。最终的无散接受指标应以有限体积的 `cleaned_flux_residual_l2` 为准；
- `curl_A_relative_l2_interior` 和 `curl_Ap_relative_l2_interior`：矢势旋度重构误差；
- `quality.gauge_spread`：DeVore 参考面取内边界和外边界所得标量的相对差异；
- `decomposition_closure`：`H_V=H_J+H_PJ` 的代数闭合误差。

若规范差异或旋度重构误差明显，应先提高均匀网格分辨率，并检查输入 DAT 的
`div B`，不能只采用一个看似合理的螺度数值。

实际 `guoguo/rbsl15` 的 PFSS、插绳场和三档分辨率结果见
[`VALIDATION.md`](VALIDATION.md)。

## 文献

- Moraitis, K., Pariat, E., Savcheva, A., & Valori, G. 2018,
  *Computation of Relative Magnetic Helicity in Spherical Coordinates*,
  Solar Physics, 293, 92, DOI: 10.1007/s11207-018-1314-5.
- Finn, J. M., & Antonsen, T. M. 1985, relative-helicity definition.
- Berger 2003; Pariat et al. 2017, gauge-invariant `H_J/H_PJ` decomposition.

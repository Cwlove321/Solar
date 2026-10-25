# solar-qsl

这是 `jet_mhd/fig.ipynb` 所用 CPU/OpenMP FastQSL 程序的可移植 Python
包装。包只保留当前计算所需内容：`qfactor` Python API、Fortran 源码、构建
脚本和可执行文件；不包含未使用的 GPU/CuPy 演示数据。

算法来源：Zhang, P., Chen, J., Liu, R. & Wang, C. (2022),
*FastQSL: A Fast Computation Method for Quasi-separatrix Layers*, ApJ 937, 26,
doi:10.3847/1538-4357/ac8d61。上游代码许可证见 `LICENSE`
（CC BY-NC-SA 4.0）。

## 依赖与构建

- Python >= 3.9
- NumPy >= 1.23
- macOS/Linux 上支持 OpenMP 的 `gfortran`；当前包内可执行文件是
  macOS 15+ ARM64，并依赖
  Homebrew GCC 的 `libgfortran`、`libgomp` 和 `libquadmath`

仓库中已包含本机编译的 ARM64 可执行文件。需要重建时：

```bash
cd /Users/xiaoxingyun/Desktop/solar_code
Packages/qsl/scripts/build_qfactor.sh
```

可选安装（使用指定的基础 Python 环境）：

```bash
/opt/anaconda3/bin/python -m pip install -e /Users/xiaoxingyun/Desktop/solar_code/Packages/qsl
```

不安装时也可把源码目录加入 `PYTHONPATH`：

```bash
export PYTHONPATH=/Users/xiaoxingyun/Desktop/solar_code/Packages/qsl/src:$PYTHONPATH
```

## Python 用法

磁场数组顺序必须是 `(nz, ny, nx)`，与 `fig.ipynb` 的
`read_binary_data()` 返回值一致。均匀网格的坐标单位是像素。

```python
from qsl import qfactor

result = qfactor(
    bx, by, bz,
    xreg=[0, 128],
    yreg=[0, 128],
    zreg=[0, 0],       # 底面（photosphere），不是 y=64 的中间切面
    delta=0.25,
    twistFlag=True,
    odir="qfactor/bottom",
    fstr="bottom_qsl_t0010_delta0.25",
)
q = result["q"]
twist = result["twist"]
```

默认调用包内 `src/qsl/bin/qfactor.x`。若需使用另一份可执行文件，设置
环境变量 `FASTQSL_EXECUTABLE=/absolute/path/to/qfactor.x`。

## `fig.ipynb` 参数对应关系

原图单元使用 `xreg=[0,128], yreg=[64,64], zreg=[0,128]`，得到的是
`y=64` 的 `x-z` 截面，输出键为 `qcs`。底面应使用
`xreg=[0,128], yreg=[0,128], zreg=[0,0]`，输出键改为 `q`；`delta=0.25`
和 `twistFlag=True` 保持不变。

本案例的可复现脚本位于
`/Users/xiaoxingyun/Desktop/jet_mhd/recalculate_bottom_qsl.py`。默认计算快照
10 和 95，并在 `jet_mhd/qfactor/bottom/` 写出两个 NPZ、统计元数据 JSON 和
预览 PNG：

```bash
/opt/anaconda3/bin/python /Users/xiaoxingyun/Desktop/jet_mhd/recalculate_bottom_qsl.py
```

预览图额外需要 Matplotlib；核心 `qsl` 包本身只依赖 NumPy。

# FastQSL2 local integration

这是 FastQSL2 的并行安装，不替换 `../qsl/`。上游版本固定在
`UPSTREAM_COMMIT` 记录的提交，原始说明见 `README.md`。

本地只对 `fastqsl.py` 做了一项运行适配：默认从本目录的
`bin/fastqsl.x` 加载 Fortran/OpenMP 核心，避免依赖作者电脑上的绝对路径。
也可以通过 `FASTQSL2_EXECUTABLE=/absolute/path/to/fastqsl.x` 覆盖。

## 依赖与构建

- Python、NumPy、Matplotlib
- 可选 SciPy（仅用于部分上游 demo）
- 支持 OpenMP 的 gfortran

从 `solar_code` 仓库根目录运行：

```bash
Packages/FastQSL2/scripts/build_fastqsl.sh
```

## Python 导入

该上游项目是单文件模块，不需要安装。脚本应相对于仓库根目录定位
`Packages/FastQSL2`，也可用 `FASTQSL2_DIR` 覆盖。

```python
import os
import sys
from pathlib import Path

fastqsl_dir = Path(os.environ.get("FASTQSL2_DIR", "Packages/FastQSL2")).resolve()
sys.path.insert(0, str(fastqsl_dir))
from fastqsl import fastqsl
```

底面计算示例：

```python
result = fastqsl(
    bx, by, bz,
    xreg=[0, 128], yreg=[0, 128], zreg=[0, 0],
    delta=0.25,
    length_out=True,
    twist_out=True,
    rF_out=True,
    silent=True,
)
```

FastQSL2 的数组顺序是 `(nz, ny, nx)`；二维返回数组顺序是 `(y, x)`。
旧 `qsl` 包的二维结果顺序是 `(x, y)`，新旧比较时应转置旧数组。

## 与旧 qsl 的关键差异

- `qfactor` 政名为 `fastqsl`
- `twistFlag` 政名为 `twist_out`
- `nbridges` 政名为 `nthreads`
- `rboundary` 改为 `10*rbs+rbe`
- `length`、`twist` 和端点坐标默认不输出，需显式启用
- `Bnr` 不再直接输出，可用 `targetB_out` 的端点磁场推导

案例脚本应通过 `--magnetic-dir`、`--output-dir` 或相应环境变量传入数据路径，
不要把个人目录写进 FastQSL2 集成层。

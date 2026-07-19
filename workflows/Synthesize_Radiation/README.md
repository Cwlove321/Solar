# 笛卡尔 AMRVAC 辐射合成

本项目集中存放从笛卡尔 AMRVAC `.dat` 快照生成合成辐射图像的工作流入口。
通用的物理量转换、探测器平面和视线积分辅助函数仍位于
`src/solarutils/sim/radiation.py`，以便测试和其他项目复用。

## 脚本说明

| 脚本 | 用途 |
| --- | --- |
| `synthesize_radiation.py` | 批量处理多个 AMRVAC 快照，支持光学薄合成和带吸收的 DDA 合成。 |
| `synthesize_aia_optically_thick.py` | 对单个 AMRVAC 快照进行 AIA 171/211/304 光学厚辐射传输计算。 |

两个脚本均依赖原作者的 `radsyn_tools`：
<https://github.com/fuwentai/radsyn_tools>。本仓库不上传 `Packages/`，
首次使用时在仓库根目录执行：

```bash
bash bootstrap_dependencies.sh
```

拉取后目录中应包含 `RST_cart/` 和 `AIA_response/`，默认使用：

```text
Packages/radsyn_tools/DDA_method1
```

整个 `solar_code` 移动到其他机器后不需要修改绝对路径。需要使用另一份辐射合成工具时，
仍可通过 `--dda-root` 覆盖默认目录。

## 批量辐射合成

`synthesize_radiation.py` 通过 glob 读取一组 AMRVAC `.dat`，将密度和能量转换为
数密度、压强和温度，再针对指定 AIA 波段执行辐射合成。

```bash
cd /path/to/solar_code
PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_radiation.py \
  --dat-glob '/Volumes/ZC/jet_mhd/*dat' \
  --output-dir '/Volumes/ZC/jet_mhd/radsyn' \
  --start 1 --stop 180 --step 3 \
  --wavelength 171 131 --save-physical --plot
```

可通过 `--normal` 设置视线方向；`--save-physical` 会保存中间物理量立方，
`--plot` 会生成 PNG 快速检查图。

## 单快照 AIA 光学厚合成

`synthesize_aia_optically_thick.py` 直接使用当前 `simesh.amrvac` 接口读取一个
AMRVAC `.dat`，计算发射率、吸收系数、总光学深度和输出强度：

```text
I = integral(exp(-tau_before) * emissivity * ds) + I0 * exp(-tau_total)
```

示例：

```bash
PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_aia_optically_thick.py \
  /path/to/mhd0006.dat \
  --wavelength 171 --normal -1 0 0 \
  --vmin 3 --vmax 10 --label example_171
```

未设置 `--output-dir` 时，输出写入输入 `.dat` 所在目录下的 `radiation/`。
每次运行会保存线性强度、`log10` 强度、总光学深度、PNG 快速检查图和元数据
JSON。没有显式设置 `--resolution` 时，默认使用 AMRVAC 最高有效分辨率。

## 检查命令行参数

```bash
PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_radiation.py --help

PYTHONPATH=src /opt/anaconda3/bin/python \
  workflows/Synthesize_Radiation/synthesize_aia_optically_thick.py --help
```

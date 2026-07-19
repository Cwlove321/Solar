# solar-code

面向太阳观测处理、AMRVAC 数值模拟和科学案例复现的长期维护仓库。

## 目录职责

- `src/solarutils/`：跨项目复用的 Python 函数和稳定 API。
- `workflows/`：可单独交付的预处理、模拟、诊断和辐射合成流程。
- `cases/`：事件数据、事件专用脚本和实际计算产物。
- `Packages/`：本地外部依赖目录，不上传到本仓库；通过原作者 Git 仓库拉取 `simesh` 和 `radsyn_tools`。
- `data/`：仅保存跨事件共享的标定数据；事件数据必须放入 `cases/<event>/data/`。
- `tests/`：通用工具、工作流共享核心和科学 I/O 的回归测试。

## 案例约定

```text
cases/<event>/
├── README.md          # 研究目标和复现步骤
├── config.yaml       # 运行参数和案例相对路径
├── data_manifest.md  # 数据来源、时间、文件与用途
├── data/             # 原始数据和必需输入
├── products/         # 中间、诊断和最终产物
└── scripts/          # 事件专用入口
```

配置和元数据中的本地文件默认使用相对于仓库根目录或案例根目录的路径，
不写入开发者的 `/Users/...` 绝对路径。

## 安装与测试

```bash
cd /path/to/solar_code
bash bootstrap_dependencies.sh
/opt/anaconda3/bin/python -m pip install -e .
/opt/anaconda3/bin/python -m pip install -e Packages/simesh
/opt/anaconda3/bin/python -m pytest
```

`bootstrap_dependencies.sh` 会从原作者仓库拉取：

- `simesh`：<https://github.com/Astery0502/simesh>
- `radsyn_tools`：<https://github.com/fuwentai/radsyn_tools>

已存在的本地 `Packages/simesh` 或 `Packages/radsyn_tools` 不会被覆盖。

大型 FITS、AMRVAC DAT、PSF 和可再生成产物默认不纳入 Git。它们的来源、服务器位置、
文件数和完整性信息应记录在对应案例的 `data_manifest.md`。

# solar-code

面向太阳观测处理、AMRVAC 数值模拟和科学案例复现的长期维护仓库。

## 目录职责

- `src/solarutils/`：跨项目复用的 Python 函数和稳定 API。
- `workflows/`：可单独交付的预处理、模拟、诊断和辐射合成流程。
- `cases/`：事件数据、事件专用脚本和实际计算产物。
- `Packages/`：科研依赖源码。带本地修正的 `simesh`、`radsyn_tools`、`qsl`
  随主仓库保存；大型或独立上游包通过锁文件固定版本。
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

需要运行图像配准、模板匹配、Movie 输出和完整测试时，安装统一的测试附加依赖：

```bash
/opt/anaconda3/bin/python -m pip install -e '.[test]'
```

`bootstrap_dependencies.sh` 会验证 vendored 包的源码哈希，并按
`Packages/dependencies.lock.json` 准备固定版本的 FastQSL2 及本地适配补丁。
它不会用上游最新版覆盖已经验证的 `simesh`、`radsyn_tools` 或 `qsl`。

可随时单独审计依赖来源、提交和源码内容：

```bash
/opt/anaconda3/bin/python scripts/audit_packages.py
```

大型 FITS、AMRVAC DAT、PSF 和可再生成产物默认不纳入 Git。它们的来源、服务器位置、
文件数和完整性信息应记录在对应案例的 `data_manifest.md`。

## 同步到 guoguo

同步脚本默认只显示差异，不写入服务器：

```bash
scripts/sync_guoguo.sh
```

确认清单后，显式执行同步和远端完整测试：

```bash
scripts/sync_guoguo.sh --apply --full-test
```

执行模式会先在远端 `backups/` 创建可恢复归档，再同步源码、做 checksum
复核，并使用 `zcc` Python 运行语法检查和测试。脚本不会删除远端独有文件，
也不会上传案例数据、Git 元数据、缓存或平台相关二进制文件。远端主机、目录和
Python 路径可分别通过 `SOLAR_CODE_REMOTE_HOST`、`SOLAR_CODE_REMOTE_DIR` 和
`SOLAR_CODE_REMOTE_PYTHON` 覆盖。

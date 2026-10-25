# Packages

这里保存工作流需要、并随 `solar_code` 一起分发的独立源码包。

版本与内容校验记录在 `dependencies.lock.json`。其中：

- `simesh/`、`radsyn_tools/`、`qsl/` 是带项目修正的 vendored 源码，纳入
  主仓库管理并使用确定性源码哈希防止静默漂移；
- `FastQSL2/` 保留独立上游 Git 历史，主仓库在 `integrations/FastQSL2/`
  保存固定补丁和集成文件；
- `amrvac/` 是可选的大型独立仓库，只锁定提交，不并入主仓库。

执行 `/opt/anaconda3/bin/python scripts/audit_packages.py` 可核对全部依赖。

- `simesh/`：新版 AMRVAC Python/Cython 网格工具。推荐在仓库根目录执行
  `/opt/anaconda3/bin/python -m pip install -e Packages/simesh`。
- `radsyn_tools/`：随仓库分发的辐射传输/DDA 工具与 AIA 响应表。
  `workflows/Synthesize_Radiation/` 会根据脚本位置自动定位 `DDA_method1`，
  不需要修改 `PYTHONPATH` 或依赖个人 `solar_package` 目录。
- `qsl/`：`jet_mhd/fig.ipynb` 使用的 CPU/OpenMP FastQSL Python 包装及
  Fortran 源码。包内可执行文件使用相对路径定位，不再依赖个人
  `pycodes/solar_package` 目录；安装、重建和底面参数说明见 `qsl/README.md`。
- `FastQSL2/`：FastQSL2 2026 的并行安装，固定上游提交并使用包内
  `bin/fastqsl.x`。它不会覆盖旧 `qsl/`；构建、导入及 API 差异见
  `FastQSL2/LOCAL_README.md`。

工作流代码应通过包的公开 API 调用这里的项目，避免依赖开发者个人目录中的副本。
`simesh` 的完整安装、构建和测试说明见 `simesh/docs/installation.md`。

`radsyn_tools` 目前作为 vendored scientific source/data bundle 使用，不单独安装。
这是因为其上游 notebook 和算法同时依赖 `RST_cart/` 与同级 `AIA_response/`；
工作流已使用相对于源文件的路径读取响应表，从任意工作目录启动都可以正常定位数据。

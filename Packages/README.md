# Packages

这里保存工作流需要、并随 `solar_code` 一起分发的独立源码包。

- `simesh/`：新版 AMRVAC Python/Cython 网格工具。推荐在仓库根目录执行
  `/opt/anaconda3/bin/python -m pip install -e Packages/simesh`。
- `radsyn_tools/`：随仓库分发的辐射传输/DDA 工具与 AIA 响应表。
  `workflows/Synthesize_Radiation/` 会根据脚本位置自动定位 `DDA_method1`，
  不需要修改 `PYTHONPATH` 或依赖个人 `solar_package` 目录。

工作流代码应通过包的公开 API 调用这里的项目，避免依赖开发者个人目录中的副本。
`simesh` 的完整安装、构建和测试说明见 `simesh/docs/installation.md`。

`radsyn_tools` 目前作为 vendored scientific source/data bundle 使用，不单独安装。
这是因为其上游 notebook 和算法同时依赖 `RST_cart/` 与同级 `AIA_response/`；
工作流已使用相对于源文件的路径读取响应表，从任意工作目录启动都可以正常定位数据。

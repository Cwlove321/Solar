# 最近时刻的 AIA / HMI FITS 与对比图

`download_aia_hmi_nearest.py` 是可复用的下载、对齐和绘图工具。目标UTC时间、视场、数据目录、图片目录及JSOC邮箱必须显式传入；事件运行命令、实测时间、下载记录与图片放在对应case中。

调用参数：

```text
--time UTC_TIME
--email JSOC_EMAIL
--wavelength AIA_WAVELENGTH
--aia-dir AIA_DIRECTORY
--hmi-dir HMI_DIRECTORY
--figure-dir FIGURE_DIRECTORY
--extent XMIN XMAX YMIN YMAX
```

使用 `/opt/anaconda3/bin/python download_aia_hmi_nearest.py --help` 查看全部参数。

- 默认查询目标时间前后90秒，按T_OBS选择QUALITY=0且最近的AIA和HMI记录；HMI文件名中的TAI不直接当作UTC。
- 使用DRMS正式FITS导出保留WCS，校验文件大小与FITS结构；数据目录JSON保存记录、时间和SHA256，并支持恢复导出。
- 使用AIA观测者位置、统一太阳半径和太阳表面差分自转，对齐到指定目标时间的同一WCS网格。
- 默认0.5角秒/像素。AIA按曝光时间归一化为DN/s，使用asinh拉伸；HMI默认为灰度线性±800G。
- `--aia-colormap q-d` 选择青、黄、橙、红色表；`--aia-vmin`、`--aia-vmax` 设置色标。`--overlay-threshold` 设置低于阈值透明，`--overlay-alpha` 控制其余像素透明度。
- 使用原始Level-1 WCS，未追加外部指向修正、PSF反卷积或辐射响应退化校正。
- 输出PNG、PDF及图像处理JSON；绘图前检查有效覆盖比例。

下载端直连JSOC官方主机，采用8段并发并检查Content-Range与文件大小。最终观测时间以正式FITS头中的T_OBS为准。

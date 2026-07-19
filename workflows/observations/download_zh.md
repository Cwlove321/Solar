# 观测数据下载方法

下载后的 AIA level-1 到 level-1.5 预处理说明见 [aia_preprocessing_zh.md](aia_preprocessing_zh.md)。
HMI 45s LOS 磁图的 submap derotate workflow 也记录在同一份预处理说明中。

这个说明用于长期维护观测数据下载工作流。当前工具统一通过
`workflows/observations/download_observations.py` 下载常用太阳观测数据，下载后的
FITS 文件默认不提交到 Git，只保存在本地 `data/` 目录中。

## 1. 环境准备

在代码库根目录运行：

```bash
pip install -e ".[aia,solar-orbiter,test]"
```

如果只下载 SDO/AIA、HMI、SHARP 数据，可以不安装 Solar Orbiter 依赖：

```bash
pip install -e ".[aia,test]"
```

数据来源如下：

- SDO/AIA、HMI：通过 SunPy/Fido 调用 JSOC；SHARP CEA：直接用
  `drms.Client().export()` 创建 JSOC 导出请求。JSOC 下载需要注册邮箱。
- STEREO-A：通过 SunPy/Fido 调用 VSO，不需要 JSOC 邮箱。
- Solar Orbiter：通过 `sunpy-soar` 调用 SOAR，可以用 Python 下载；如果 SOAR 服务不可用，再从 Solar Orbiter Archive 网站手动下载。

## 2. 推荐目录结构

普通科研数据建议放在：

```text
data/raw/
├── sdo/
│   ├── aia/
│   ├── hmi/
│   └── sharp/
├── stereo_a/
└── solar_orbiter/
```

用于测试代码正确性的少量数据建议放在：

```text
cases/<YYYYMMDD>/data/<time_range>/
```

例如这次已经下载好的测试数据：

```text
cases/20241008/data/2300_2400/
```

## 3. 基本命令格式

所有下载命令都有同一组常用参数：

```bash
python workflows/observations/download_observations.py <data-source> \
  --start <开始时间> \
  --end <结束时间> \
  --out <输出目录> \
  --fetch
```

常用参数含义：

- `--start`：开始时间，例如 `2024-10-08T23:00:00`。
- `--end`：结束时间，例如 `2024-10-09T00:00:00`。如果你说 `23:00-24:00`，代码里应写成第二天 `00:00:00`。
- `--sample-minutes`：按多少分钟抽样，例如 `10` 表示大约每 10 分钟取一帧。
- `--out`：下载目录。
- `--fetch`：真正下载文件；不加这个参数时只查询，不下载。
- `--overwrite`：覆盖已经存在的文件，默认不覆盖。
- `--email`：JSOC 注册邮箱，只对 SDO/AIA、HMI、SHARP 必需。

建议先不加 `--fetch` 查询一次，看返回结果是否合理；确认数量和时间范围正确后再加
`--fetch` 下载。

## 4. 下载 SDO/AIA

AIA EUV 常用波段包括 `94 131 171 193 211 304 335` Angstrom。

```bash
python workflows/observations/download_observations.py sdo-aia \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --wavelengths 94 131 171 193 211 304 335 \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/sdo/aia \
  --fetch
```

说明：

- 当前默认下载 `image` segment。
- AIA 原始时间采样很密，做测试数据时一定建议加 `--sample-minutes`。
- 文件时间可能与整点差几秒，这是仪器实际观测时间和 Fido 取最近样本造成的，属于正常情况。

## 5. 下载 HMI LOS 磁图

下载 45 秒线-of-sight 磁图：

```bash
python workflows/observations/download_observations.py hmi \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --series hmi_m_45s \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/sdo/hmi_45s \
  --fetch
```

下载 720 秒线-of-sight 磁图：

```bash
python workflows/observations/download_observations.py hmi \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --series hmi_m_720s \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/sdo/hmi_720s \
  --fetch
```

说明：

- `hmi_m_45s` 对应 JSOC series `hmi.M_45s`。
- `hmi_m_720s` 对应 JSOC series `hmi.M_720s`。
- 为了兼容以前的记法，代码中也把 `hmi.75` 映射到了 `hmi.M_720s`。
- 默认 segment 是 `magnetogram`。

## 6. 下载 SHARP 矢量磁场

SHARP CEA 数据需要先确定 `HARPNUM`。如果已经知道目标活动区的 `HARPNUM`，可以直接下载：

```bash
python workflows/observations/download_observations.py sharp \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --harpnum 12001 \
  --sample-minutes 10 \
  --segments Br Bp Bt bitmap conf_disambig \
  --out cases/20241008/data/2300_2400/sdo/sharp_harp12001 \
  --fetch
```

`sharp --fetch` 会直接提交 DRMS export 请求，不走 SunPy/Fido 下载层。例如本次
2015-08-27 05:24 TAI 的 HARP 5885：

```bash
python workflows/observations/download_observations.py sharp \
  --start 2015-08-27T05:23:00 \
  --end 2015-08-27T05:25:00 \
  --email 181830256@smail.nju.edu.cn \
  --harpnum 5885 \
  --segments Br Bt Bp bitmap conf_disambig \
  --out data/raw/sdo/sharp/jsoc_export_drms_hmi.sharp_cea_720s_5885_20150827_0524 \
  --fetch
```

常用 segment：

- `Br`：径向磁场分量。
- `Bp`：纬向磁场分量。
- `Bt`：经向磁场分量。
- `bitmap`：活动区掩膜。
- `conf_disambig`：消歧置信度信息。

说明：

- 不指定 `--harpnum` 时可能会返回很多活动区，测试数据会变得非常大。
- 当前测试数据使用 `HARPNUM=12001`。

## 7. 下载 HMI 全日面 720s 矢量磁场

如果需要全日面的三维矢量磁场，可以下载 JSOC series `hmi.B_720s`。常用 segment 是：

- `field`：磁场强度。
- `inclination`：磁场倾角。
- `azimuth`：磁场方位角。
- `disambig`：180 度方位角消歧信息。

命令：

```bash
python workflows/observations/download_observations.py hmi-vector \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/sdo/hmi_b_720s \
  --fetch
```

也可以用通用 `hmi` 命令等价查询：

```bash
python workflows/observations/download_observations.py hmi \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --series hmi_b_720s \
  --segments field inclination azimuth disambig \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/sdo/hmi_b_720s \
  --fetch
```

## 8. 下载 HMI synoptic-frame 径向磁场

你给出的文件：

```text
hmi.mrdailysynframe_polfil_720s.20241008_120000_TAI.Mr_polfil.fits
```

对应 JSOC series `hmi.mrdailysynframe_polfil_720s`，segment 是 `Mr_polfil`。
可以用下面的命令下载当天 12:00 TAI 附近的文件：

```bash
python workflows/observations/download_observations.py hmi-synoptic \
  --start 2024-10-08T11:50:00 \
  --end 2024-10-08T12:10:00 \
  --email your_jsoc_email@example.com \
  --out cases/20241008/data/2300_2400/sdo/hmi_mrdailysynframe_polfil_720s \
  --fetch
```

如果只想精确查某个时刻，也可以把 `--start` 和 `--end` 缩得更近一些。

Carrington rotation synoptic 文件：

```text
hmi.synoptic_mr_polfil_720s.2289.Mr_polfil
```

对应 JSOC series `hmi.synoptic_mr_polfil_720s`，prime key 是 `CAR_ROT=2289`。
可以这样下载：

```bash
python workflows/observations/download_observations.py hmi-synoptic-cr \
  --car-rotation 2289 \
  --email your_jsoc_email@example.com \
  --out cases/20241008/data/2300_2400/sdo/hmi_synoptic_mr_polfil_720s \
  --fetch
```

如果只给一个时间，工具会自动计算该时间所在的 Carrington rotation：

```bash
python workflows/observations/download_observations.py hmi-synoptic-cr \
  --time 2024-10-08T12:00:00 \
  --email your_jsoc_email@example.com \
  --out cases/20241008/data/2300_2400/sdo/hmi_synoptic_mr_polfil_720s \
  --fetch
```

## 9. 下载 STEREO-A

下载 STEREO-A/SECCHI EUVI 195 Angstrom：

```bash
python workflows/observations/download_observations.py stereo-a \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --instrument SECCHI \
  --detector EUVI \
  --wavelength 195 \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/stereo_a/euvi_195 \
  --fetch
```

说明：

- STEREO-A 通过 VSO 下载，不需要 JSOC 邮箱。
- `--wavelength 195` 会转换为 SunPy 的 Angstrom 波长约束。

## 10. 下载 Solar Orbiter

Solar Orbiter 可以通过 Python 下载，前提是安装了 `sunpy-soar`。

下载 EUI/FSI 174 Angstrom L2 图像：

```bash
python workflows/observations/download_observations.py solar-orbiter \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --instrument EUI \
  --level L2 \
  --product eui-fsi174-image \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/solar_orbiter/eui_fsi174 \
  --fetch
```

说明：

- `--product eui-fsi174-image` 是 SOAR 中的 EUI FSI 174 图像产品。
- SOAR 不支持 Fido 原生的 `Sample` 查询约束，所以本工具会先查询完整结果，再在 Python 中按 `--sample-minutes` 抽样。
- 如果 Python 下载失败，可以从 Solar Orbiter Archive 网站手动下载，并按同样目录结构放入 `cases/<YYYYMMDD>/data/<time_range>/solar_orbiter/` 或 `cases/<YYYYMMDD>/data/`。

## 11. 校验下载的数据

下载完成后，建议立刻打开少量文件确认 FITS 不是坏文件：

```bash
python workflows/observations/validate_downloads.py \
  cases/20241008/data/2300_2400/sdo/aia \
  --max-files 4 \
  --sunpy-map
```

HMI、STEREO-A、Solar Orbiter 也可以用同一个校验脚本：

```bash
python workflows/observations/validate_downloads.py \
  cases/20241008/data/2300_2400/stereo_a/euvi_195 \
  --max-files 2 \
  --sunpy-map

python workflows/observations/validate_downloads.py \
  cases/20241008/data/2300_2400/solar_orbiter/eui_fsi174 \
  --max-files 2 \
  --sunpy-map
```

SHARP CEA 文件也可以先用 Astropy 打开：

```bash
python workflows/observations/validate_downloads.py \
  cases/20241008/data/2300_2400/sdo/sharp_harp12001 \
  --max-files 5
```

## 10. 当前测试数据集

目前已经下载好的测试数据位于：

```text
cases/20241008/data/2300_2400/
```

内容包括：

- SDO/AIA：49 个 FITS。
- HMI 45s：7 个 FITS。
- HMI 720s：6 个 FITS。
- SHARP `HARPNUM=12001`：30 个 FITS。
- STEREO-A EUVI 195 Angstrom：6 个 FITS。
- Solar Orbiter EUI FSI 174 Angstrom L2：6 个 FITS。

总计 104 个数据文件，约 815 MB。这个数据集的目的不是完整科学分析，而是用于测试下载、
读入、Map 构造、裁剪、配准和画图等后续工具是否工作正常。

## 11. 后续维护建议

每个 case 研究建议记录一个 `data_manifest.md`，至少写清楚：

- 事件时间范围。
- 下载命令。
- 数据源、波段、series、segment。
- SHARP 的 `HARPNUM` 和对应 NOAA 活动区编号，如果已知。
- 下载日期。
- 是否做过重采样、裁剪、配准或坐标变换。

这样以后复现实验或写文章时，不需要靠记忆追溯数据来源。

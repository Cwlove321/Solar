# Observation Download Workflow

Observation workflow overview: [overview.md](overview.md).
中文下载说明见 [download_zh.md](download_zh.md).
AIA preprocessing notes: [English](aia_preprocessing.md), [中文](aia_preprocessing_zh.md).

Processed HMI 45s LOS submap workflow:

```bash
PYTHONPATH=src python workflows/observations/process_hmi45_derotate_submaps.py --method pixel --overwrite
PYTHONPATH=src python workflows/observations/plot_processed_hmi45_submaps.py --method pixel --overwrite
```

CHASE Halpha Doppler workflow:

```bash
PYTHONPATH=src SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mpl_config XDG_CACHE_HOME=.mpl_config/cache \
python workflows/observations/process_chase_spectrum.py
```

By default this plots selected velocity channels and calculates a centroid Doppler velocity map with
differential-rotation LOS correction enabled. Use `--no-solar-rotation-correction` to save the raw
centroid Doppler map instead.

This workflow is the first step for observational case studies: query and download raw SDO,
STEREO-A, and Solar Orbiter data into a predictable directory layout.

## Data sources

- SDO/AIA and HMI use JSOC through SunPy/Fido. SHARP CEA uses direct
  `drms.Client().export()` requests. JSOC downloads require a registered
  JSOC email address.
- STEREO-A uses SunPy/Fido through VSO.
- Solar Orbiter can be downloaded from Python with the optional `sunpy-soar` package, which
  registers the Solar Orbiter Archive (SOAR) client with Fido.

## Install

From the repository root:

```bash
pip install -e ".[aia,solar-orbiter,test]"
```

For SDO-only work, `sunpy-soar` is not required:

```bash
pip install -e ".[aia,test]"
```

## Recommended output layout

```text
data/raw/
├── sdo/
│   ├── aia/
│   ├── hmi/
│   └── sharp/
├── stereo_a/
└── solar_orbiter/
```

The repository `.gitignore` excludes FITS and other large data products. Keep these files local and
record important downloads in the relevant case `data_manifest.md`.

Use `--sample-minutes` to thin dense search results. JSOC and VSO use SunPy/Fido sampling directly;
Solar Orbiter/SOAR is sampled after the query because SOAR does not accept the Fido `Sample`
attribute.

## Examples

Search SDO/AIA 171 and 304 Angstrom without downloading:

```bash
python workflows/observations/download_observations.py sdo-aia \
  --start 2024-10-08T22:00:00 \
  --end 2024-10-08T23:00:00 \
  --email your_jsoc_email@example.com \
  --wavelengths 171 304
```

Fetch the results by adding `--fetch`:

```bash
python workflows/observations/download_observations.py sdo-aia \
  --start 2024-10-08T22:00:00 \
  --end 2024-10-08T23:00:00 \
  --email your_jsoc_email@example.com \
  --wavelengths 171 304 \
  --fetch
```

Search HMI 45-second LOS magnetograms:

```bash
python workflows/observations/download_observations.py hmi \
  --start 2024-10-08T22:00:00 \
  --end 2024-10-08T23:00:00 \
  --email your_jsoc_email@example.com \
  --series hmi_m_45s
```

Search HMI 720-second LOS magnetograms:

```bash
python workflows/observations/download_observations.py hmi \
  --start 2024-10-08T22:00:00 \
  --end 2024-10-08T23:00:00 \
  --email your_jsoc_email@example.com \
  --series hmi_m_720s
```

Search full-disk HMI 720-second vector magnetic-field segments:

```bash
python workflows/observations/download_observations.py hmi-vector \
  --start 2024-10-08T22:00:00 \
  --end 2024-10-08T23:00:00 \
  --email your_jsoc_email@example.com \
  --segments field inclination azimuth disambig
```

Search HMI daily synoptic-frame radial-field maps, for example the file
`hmi.mrdailysynframe_polfil_720s.20241008_120000_TAI.Mr_polfil.fits`:

```bash
python workflows/observations/download_observations.py hmi-synoptic \
  --start 2024-10-08T11:50:00 \
  --end 2024-10-08T12:10:00 \
  --email your_jsoc_email@example.com \
  --segments Mr_polfil
```

Search HMI Carrington-rotation synoptic radial-field maps, for example
`hmi.synoptic_mr_polfil_720s.2289.Mr_polfil`:

```bash
python workflows/observations/download_observations.py hmi-synoptic-cr \
  --car-rotation 2289 \
  --email your_jsoc_email@example.com \
  --segments Mr_polfil
```

You can also pass a time and let the workflow infer the integer Carrington rotation:

```bash
python workflows/observations/download_observations.py hmi-synoptic-cr \
  --time 2024-10-08T12:00:00 \
  --email your_jsoc_email@example.com \
  --segments Mr_polfil
```

Search SHARP CEA vector data for a known HARPNUM:

```bash
python workflows/observations/download_observations.py sharp \
  --start 2024-10-08T22:00:00 \
  --end 2024-10-08T23:00:00 \
  --email your_jsoc_email@example.com \
  --harpnum 12345
```

Add `--fetch` to submit a direct DRMS export request and download the selected
segments as FITS files:

```bash
python workflows/observations/download_observations.py sharp \
  --start 2015-08-27T05:23:00 \
  --end 2015-08-27T05:25:00 \
  --email your_jsoc_email@example.com \
  --harpnum 5885 \
  --segments Br Bt Bp bitmap conf_disambig \
  --out data/raw/sdo/sharp/hmi.sharp_cea_720s_5885_20150827_0524 \
  --fetch
```

Search STEREO-A/SECCHI EUVI 195 Angstrom:

```bash
python workflows/observations/download_observations.py stereo-a \
  --start 2024-10-08T22:00:00 \
  --end 2024-10-08T23:00:00 \
  --wavelength 195
```

Search Solar Orbiter EUI data through SOAR:

```bash
python workflows/observations/download_observations.py solar-orbiter \
  --start 2024-10-08T22:00:00 \
  --end 2024-10-08T23:00:00 \
  --instrument EUI \
  --product EUI-FSI174-IMAGE \
  --level L2
```

If `sunpy-soar` is not installed or SOAR is unavailable, download manually from the Solar Orbiter
Archive website and place files under `data/raw/solar_orbiter/<case-name>/`.

## Validate downloaded files

After downloading, open a few files to verify they are readable:

```bash
python workflows/observations/validate_downloads.py data/raw/sdo/aia --sunpy-map
```

For Solar Orbiter or STEREO FITS files, the same validator works:

```bash
python workflows/observations/validate_downloads.py data/raw/solar_orbiter --max-files 5
```

For `.cdf` files, install `cdflib` first:

```bash
pip install cdflib
python workflows/observations/validate_downloads.py data/raw/solar_orbiter
```

## Test Dataset: 2024-10-08 23:00-24:00

The local smoke-test data for the first observation workflow is stored under:

```text
cases/20241008/data/2300_2400/
```

It covers `2024-10-08T23:00:00` to `2024-10-09T00:00:00` with roughly 10-minute sampling. The
current dataset contains SDO/AIA EUV images, HMI 45-second and 720-second LOS magnetograms,
SHARP CEA vector-field segments for `HARPNUM=12001`, STEREO-A/SECCHI EUVI 195 Angstrom, and
Solar Orbiter/EUI FSI 174 Angstrom L2 images.

Example commands:

```bash
python workflows/observations/download_observations.py sdo-aia \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --wavelengths 94 131 171 193 211 304 335 \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/sdo/aia \
  --fetch

python workflows/observations/download_observations.py hmi \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --series hmi_m_45s \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/sdo/hmi_45s \
  --fetch

python workflows/observations/download_observations.py hmi \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --series hmi_m_720s \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/sdo/hmi_720s \
  --fetch

python workflows/observations/download_observations.py sharp \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --email your_jsoc_email@example.com \
  --harpnum 12001 \
  --sample-minutes 10 \
  --segments Br Bp Bt bitmap conf_disambig \
  --out cases/20241008/data/2300_2400/sdo/sharp_harp12001 \
  --fetch

python workflows/observations/download_observations.py stereo-a \
  --start 2024-10-08T23:00:00 \
  --end 2024-10-09T00:00:00 \
  --instrument SECCHI \
  --detector EUVI \
  --wavelength 195 \
  --sample-minutes 10 \
  --out cases/20241008/data/2300_2400/stereo_a/euvi_195 \
  --fetch

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

## References

- SunPy Fido acquiring data guide: https://docs.sunpy.org/en/stable/tutorial/acquiring_data/index.html
- SunPy JSOC client behavior: https://docs.sunpy.org/en/latest/_modules/sunpy/net/jsoc/jsoc.html
- SunPy wavelength searches: https://docs.sunpy.org/en/stable/api/sunpy.net.attrs.Wavelength.html
- sunpy-soar SOAR plugin: https://docs.sunpy.org/projects/soar/en/stable/api.html
- sunpy-soar query examples: https://docs.sunpy.org/projects/soar/en/latest/dev_guide/working.html

## 自定义 SHARP 类 CEA 与 AIA 投影

磁场生成与 AIA 共用网格、中心及参考时间的工作流，见
[README_custom_cea.md](README_custom_cea.md)。支持角秒范围选择、官方 JSOC 数值内核、
完整分辨率 FITS/NPY、AIA 曝光归一化及差分自转；bin2 改为显式选项。

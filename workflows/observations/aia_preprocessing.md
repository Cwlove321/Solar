# AIA Preprocessing

English notes for the maintained SDO/AIA level-1 to level-1.5 preprocessing
workflow. 中文说明见 [aia_preprocessing_zh.md](aia_preprocessing_zh.md).

## Purpose

AIA level-1 FITS files are useful for quick inspection. Science workflows often
need corrected pointing, standard plate scale, standard orientation, and sometimes
PSF deconvolution. The maintained helper is:

```python
from solarutils.obs.maps import aia_prep
```

The processing order is:

1. `update_pointing`: correct AIA pointing metadata.
2. optional `respike`: restore despiked pixels if spike positions and values are supplied.
3. optional `aiapy.psf.deconvolve`: deconvolve the AIA point-spread function.
4. `register`: rotate and resample onto the standard AIA level-1.5 grid.

PSF deconvolution is intentionally optional and defaults to off in `solarutils`.
Enable it explicitly when needed.

## Install

```bash
pip install -e ".[aia,test]"
```

## PSF Cache

The AIA PSF is a full-disk `4096 x 4096` array. Recomputing it is slow, so the
test workflow caches common EUV channels locally:

```text
data/calibration/sdo/aia/psf/
├── aia_0094_psf.npz
├── aia_0131_psf.npz
├── aia_0171_psf.npz
├── aia_0193_psf.npz
├── aia_0211_psf.npz
├── aia_0304_psf.npz
└── aia_0335_psf.npz
```

Generate or refresh the cache:

```bash
PYTHONPATH=src python tests/visual/cache_aia_psf_examples.py
```

Use cached PSFs in code:

```python
from solarutils.obs.maps import aia_psf_path, cache_aia_psfs, load_aia_psf

cache_aia_psfs([171, 193])
psf_path = aia_psf_path(193)
psf = load_aia_psf(193)
```

## Basic Usage

Registration without PSF deconvolution:

```python
import astropy.units as u
from sunpy.map import Map
from aiapy.calibrate.utils import get_pointing_table

from solarutils.obs.maps import aia_prep

aia_map = Map("aia.lev1_euv_12s.2024-10-08T230006Z.193.image_lev1.fits")
pointing_table = get_pointing_table(
    "JSOC",
    time_range=(aia_map.date - 12 * u.h, aia_map.date + 12 * u.h),
)

prepared = aia_prep(aia_map, pointing_table=pointing_table)
```

Registration with PSF deconvolution:

```python
from solarutils.obs.maps import aia_prep, aia_psf_path

prepared = aia_prep(
    aia_map,
    pointing_table=pointing_table,
    deconvolve_psf=True,
    psf=aia_psf_path(193),
    deconvolve_iterations=25,
)
```

For a lightweight smoke test:

```python
prepared = aia_prep(
    aia_map,
    correct_pointing=False,
    deconvolve_psf=False,
    register_map=True,
)
```

## Batch Recommendation

When preparing many AIA maps from the same time range:

- Fetch the pointing table once and reuse it.
- Use cached PSF paths instead of recalculating PSFs.
- Use `deconvolve_iterations=25` for formal processing.
- Use smaller values such as `5` only for quick visual tests.

```python
from pathlib import Path

import astropy.units as u
from sunpy.map import Map
from aiapy.calibrate.utils import get_pointing_table

from solarutils.obs.maps import aia_prep, aia_psf_path

files = sorted(Path("cases/20241008/data/2300_2400/sdo/aia").glob("*193*.fits"))
first_map = Map(files[0])
pointing_table = get_pointing_table(
    "JSOC",
    time_range=(first_map.date - 12 * u.h, first_map.date + 12 * u.h),
)
psf = aia_psf_path(193)

prepared_maps = [
    aia_prep(
        Map(path),
        pointing_table=pointing_table,
        deconvolve_psf=True,
        psf=psf,
    )
    for path in files
]
```

## Local Submap Deconvolution

Local PSF deconvolution can be much faster, but it is an approximation. A direct
submap deconvolution with the full `4096 x 4096` PSF fails because the image and
PSF FFT shapes do not match. Use `crop_psf_to_shape` to center-crop and normalize
the PSF to the submap shape.

The target submap shape must have even dimensions because `aiapy.psf.deconvolve`
uses `rfft2/irfft2`.

```python
from sunpy.map import Map
from aiapy.psf import deconvolve

from solarutils.obs.maps import crop_psf_to_shape, get_submap, load_aia_psf

submap = get_submap(aia_map, (100, 800, 150, 500))
submap = Map(
    submap.data[: submap.data.shape[0] // 2 * 2, : submap.data.shape[1] // 2 * 2],
    submap.meta.copy(),
)
psf = crop_psf_to_shape(load_aia_psf(193), submap.data.shape)
local_deconvolved = deconvolve(submap, psf=psf, iterations=5)
```

The recommended helper crops a larger padded region, deconvolves that local
region, and returns only the central target extent:

```python
from solarutils.obs.maps import deconvolve_submap_with_padding, aia_psf_path

local = deconvolve_submap_with_padding(
    aia_map,
    extent=(100, 800, 150, 500),
    psf=aia_psf_path(193),
    padding=200,
    iterations=25,
)
```

`psf` is optional; if omitted, it is calculated from the map wavelength. The
default padding is `200 arcsec`.

A test on the AIA 193 submap `(100, 800, 150, 500)` compared padded local
deconvolution against full-frame deconvolution followed by cropping:

```text
padding=0      mean |diff|=3.703   max=127.384
padding=100    mean |diff|=0.757   max=10.851
padding=200    mean |diff|=0.436   max=6.569
padding=400    mean |diff|=0.193   max=1.722
```

Use local deconvolution for previews and experiments. For final quantitative
analysis, prefer full-frame deconvolution, or use a padded local region and crop
the padding away afterward.

Visual comparison:

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/local_psf_deconvolution_compare.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/local_psf_padding_compare.png
```

## Derotation

For time sequences, older case notebooks commonly used a `de_rotate` helper to
crop each map to the solar-rotated center of a reference submap. The maintained
equivalent is:

```python
from solarutils.obs.maps import derotate_to_reference, aia_psf_path

aligned = derotate_to_reference(
    path,
    reference_map,
    aia_prep=True,
    deconvolve_psf=True,
    psf=aia_psf_path(193),
    deconvolve_iterations=25,
    method="pixel",
)
```

`method="pixel"` is the default because it rotates every reference pixel to the
source time and then interpolates the source map. This is more accurate for large
submaps and longer time gaps, while the local benchmark showed only a modest
speed cost. `method="center"` rotates only the reference submap center and is
still available for quick previews:

```python
aligned = derotate_to_reference(
    path,
    reference_map,
    aia_prep=True,
    deconvolve_psf=True,
    psf=aia_psf_path(193),
    method="center",
)
```

The returned map uses the reference map's FITS metadata for WCS/geometry, then
replaces the observation time keywords with the source map's own time. The data
are divided by the source exposure time and then multiplied by the reference
exposure time, so the output frames are on the same reference-exposure DN scale.
`EXPTIME` is kept consistent with the reference map. The solar-rotation
calculation uses the source map's explicit observer coordinate, which avoids
SunPy's Earth-observer warning.

FITS smoke test:

```bash
PYTHONPATH=src python tests/visual/verify_aia_derotate.py --overwrite
```

Output:

```text
cases/20241008/products/processed/observations/2300_2400/sdo/aia/193/derotate_20241008T231006_to_230006.fits
```

Batch-process the local test AIA dataset:

```bash
PYTHONPATH=src python workflows/observations/process_aia_prep_derotate_submaps.py --overwrite
```

Choose the derotation mode with `--method center|pixel`. The default is now
`pixel`, which writes to a separate product directory. Use `--method center`
only when you want the historical fast approximation:

```bash
PYTHONPATH=src python workflows/observations/process_aia_prep_derotate_submaps.py \
  --method center \
  --overwrite
```

This prepares each AIA map with `aia_prep(deconvolve_psf=True)`, uses the first
frame of each wavelength as the reference submap, derotates every frame to that
reference, and saves FITS files under:

```text
cases/20241008/products/processed/observations/2300_2400/sdo/aia/<wavelength>/prep_psf_derotate_submap/
cases/20241008/products/processed/observations/2300_2400/sdo/aia/<wavelength>/prep_psf_derotate_submap_pixel/
```

The current processed test database contains 49 FITS files, all with shape
`(584, 1168)`, no `BLANK` keyword, and valid `BUNIT='DN'`.

Plot processed FITS quicklooks:

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py --overwrite
```

The plotting workflow also defaults to pixel-mode products. Plot center-mode
products with:

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --derotate-method center \
  --overwrite
```

By default, each wavelength uses the first frame's `1` and `99.5` percentiles
as a common `vmin/vmax` for all frames. Override this with:

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --norm-source first \
  --percentiles 1 99.5 \
  --overwrite
```

Use `--norm-source all` to calculate percentiles from all frames of each
wavelength.

Optional sharpening is applied before percentile calculation and plotting, so
`vmin/vmax` correspond to the displayed sharpened data:

```bash
PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --sharpen none \
  --percentiles 1 99.5 \
  --overwrite

PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --sharpen single \
  --sigma 1 \
  --weight 1.5 \
  --percentiles 1 99.5 \
  --overwrite

PYTHONPATH=src python workflows/observations/plot_processed_aia_submaps.py \
  --sharpen multiscale \
  --sigmas 1 2 4 \
  --weights 0.8 0.35 0.15 \
  --percentiles 1 99.5 \
  --overwrite
```

Sharpened quicklooks are written under `single_sharpen/` or
`multiscale_sharpen/` inside the processed figure directory.

This writes 49 single-frame PNGs and 7 wavelength sequence PNGs under:

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/processed/
```

Create movies from quicklook PNG frames:

```bash
PYTHONPATH=src python workflows/observations/make_processed_aia_movies.py \
  --variant multiscale_sharpen \
  --backend imageio \
  --fps 3 \
  --quality simple \
  --overwrite
```

Use `--quality high` for a clearer, larger MP4. The default `simple` preset is
intended for compact preview files.

The reusable helper is:

```python
from solarutils.plotting import make_movie_from_images
```

`backend='auto'` uses system `ffmpeg` when available and falls back to `imageio`.
The imageio backend pads odd frame dimensions to even dimensions before writing
MP4 files. Both backends support the `simple` and `high` quality presets.

## HMI 45s LOS Submap Derotation

HMI 45-second LOS magnetograms do not use the AIA pointing/PSF/register steps.
The maintained workflow rotates the first HMI map, crops the reference submap,
derotates every later HMI map to that reference, and saves FITS files with the
reference WCS/geometry and each source frame's own time keywords. Use
`--method pixel` for per-pixel differential-rotation remapping, or
`--method center` for the faster center-point approximation.

```bash
PYTHONPATH=src python workflows/observations/process_hmi45_derotate_submaps.py --method pixel --overwrite
PYTHONPATH=src python workflows/observations/plot_processed_hmi45_submaps.py --method pixel --overwrite
```

Outputs:

```text
cases/20241008/products/processed/observations/2300_2400/sdo/hmi_45s/los_derotate_submap_pixel/
cases/20241008/products/figures/observations/2300_2400/sdo/hmi_45s/los_derotate_submap_pixel/
```

The current test output contains 7 FITS files with shape `(695, 1390)`,
`BUNIT='Gauss'`, no `BLANK` keyword, and no `EXPTIME` keyword.

The same unified helper is used by AIA and HMI:

```python
from solarutils.obs.maps import derotate_to_reference

aligned = derotate_to_reference(path, reference_map)
aligned = derotate_to_reference(path, reference_map, method="center")
```

Local benchmark on the `(100, 800, 150, 500)` arcsec test submap:

```text
hmi45 center best=4.645640 s  mean=4.773800 s
hmi45 pixel  best=4.857958 s  mean=5.055192 s
aia193 center best=0.346969 s  mean=0.364392 s
aia193 pixel  best=0.635687 s  mean=0.673612 s
```

Benchmark script/output:

```text
tests/visual/benchmark_derotation_methods.py
data/benchmarks/derotation_methods.txt
```

## Sharpening

Older case notebooks also used a lightweight unsharp-mask step after derotation.
The maintained helper is:

```python
from solarutils.obs.maps import multiscale_sharpen_map, sharpen_map

sharpened = sharpen_map(aligned, sigma=1, weight=1.5)
multiscale = multiscale_sharpen_map(
    aligned,
    sigmas=(1, 2, 4),
    weights=(0.8, 0.35, 0.15),
)
```

The single-scale formula is:

```text
sharpened = data + (data - gaussian_filter(data, sigma)) * weight
```

`multiscale_sharpen_map` applies the same unsharp-mask detail extraction at
several Gaussian scales and adds the weighted details back to the original map.
The default weights are intentionally mild so that the result is useful for
visual inspection without looking over-processed.

Use it as a visual enhancement step, not as a replacement for PSF
deconvolution. The visual smoke test is:

```bash
PYTHONPATH=src python tests/visual/verify_aia_sharpen.py
PYTHONPATH=src python tests/visual/verify_aia_sharpen.py --multiscale \
  --output cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/multiscale_sharpen_compare.png
```

Output:

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/sharpen_compare.png
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/multiscale_sharpen_compare.png
```

## Validation

The maintained `aia_prep` was compared against
`/path/to/193.ipynb`:

```text
update_pointing -> aiapy.psf.deconvolve -> register
```

Using the same AIA 193 FITS file, JSOC pointing table, cached 193 PSF, and 5
deconvolution iterations:

```text
max_abs_difference = 0
```

Comparison figure:

```text
cases/20241008/products/figures/observations/2300_2400/sdo/aia/193/notebook_flow_compare.png
```

## Notes

- Full-disk PSF deconvolution is memory- and CPU-intensive.
- Cached PSFs are local calibration/test data and are not meant to be committed to Git.
- `respike` is optional and skipped unless `(positions, values)` are supplied.
- `register` changes the map to the standard `0.6 arcsec / pix` scale and removes
  the small roll angle in the level-1 map.

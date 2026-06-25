# Observation Workflow Overview

This directory contains reusable command-line workflows for observational data.
Case-specific paths and fixed science choices should live under `cases/<case-id>/`.

## Boundary

- `src/solarutils`: reusable readers, coordinate handling, calibration helpers, measurements, and plotting functions.
- `workflows/observations`: generic scripts that expose reusable functions through command-line arguments.
- `cases/<case-id>`: event-specific configuration, fixed paths, run notes, and case scripts.
- `data/` and external drives: large local FITS files, generated arrays, movies, and temporary figures.

## Current Pipeline

1. Download and validate raw observations.
2. Prepare AIA level-1 data toward level 1.5 when needed.
3. Cut submaps and derotate them to a reference time.
4. Generate quick-look maps, context/submap figures, and movies.
5. Measure derived quantities:
   - HMI magnetic flux in selected regions.
   - AIA time-distance diagrams.
   - CHASE Halpha Doppler channels and centroid Doppler velocity maps.
6. Move event-specific runs into `cases/<case-id>` when a workflow becomes part of a formal case study.

## Naming

- Keep descriptive FITS filenames for processed products when they encode important provenance.
- Add processing switches to derived diagnostic names when ambiguity matters, for example
  `doppler_centroid_rotcorr`.
- Avoid committing generated FITS, NPZ, PNG, MP4, or cache products unless a final case figure is
  intentionally placed under `cases/<case-id>/figures/final/`.

## Environment

Recommended local cache variables:

```bash
PYTHONPATH=src SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mpl_config XDG_CACHE_HOME=.mpl_config/cache
```

These keep SunPy, Astropy, and Matplotlib caches inside ignored project-local directories.

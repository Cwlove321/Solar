# Spherical Simulation Workflows

Workflows for global or shell simulations on spherical grids.

## Scope

- Grids in `(r, theta, phi)`, `(r, latitude, longitude)`, or Carrington longitude/latitude.
- Global coronal/MHD models.
- Synoptic-map driven or magnetogram-driven boundary conditions.

## Expected Data Products

- Spherical shell scalar fields.
- Vector fields in spherical components: `Br`, `Btheta`, `Bphi`, and velocity equivalents.
- Boundary products from HMI synoptic maps and full-disk vector magnetograms.

## First Workflows To Add

- Read HMI synoptic `Mr_polfil` maps as spherical lower boundaries.
- Convert longitude/latitude grids into simulation coordinate arrays.
- Slice spherical data at fixed radius, longitude, or latitude.
- Convert spherical vector components to Cartesian components for comparison plots.
- Export Carrington-map quicklooks and movies.

## RBSL/PFSS Boundary Preparation

`prepare_rbsl_boundary.py` reproduces the spherical-boundary preparation flow from the
2024-10-08 RBSL notebooks:

- plot an HMI context submap and draw the Stonyhurst simulation region;
- extract `Br`, `Btheta`, and `Bphi` on a regular lon/lat grid;
- plot the local vector field and report the spherical flux-balance diagnostic;
- optionally preprocess the vector boundary to reduce force/torque and smooth the field;
- replace the matching `mrdailysynframe_polfil` patch with the processed `Br`;
- write `pfss_boundary.dat`, `cutout.dat`, and the parameter snippets for
  `pfss_mod_usr.t` and `pfss_amrvac.par`.

Example for the 2024-10-08 case data:

```bash
PYTHONPATH=src python workflows/simulations/spherical/prepare_rbsl_boundary.py \
  --field /Users/xiaoxingyun/Desktop/20241008/cme/rbsl/data/hmi.B_720s.20241008_230000_TAI.field.fits \
  --inclination /Users/xiaoxingyun/Desktop/20241008/cme/rbsl/data/hmi.B_720s.20241008_230000_TAI.inclination.fits \
  --azimuth /Users/xiaoxingyun/Desktop/20241008/cme/rbsl/data/hmi.B_720s.20241008_230000_TAI.azimuth.fits \
  --disambig /Users/xiaoxingyun/Desktop/20241008/cme/rbsl/data/hmi.B_720s.20241008_230000_TAI.disambig.fits \
  --synoptic /Users/xiaoxingyun/Desktop/20241008/cme/rbsl/data/hmi.mrdailysynframe_polfil_720s.20241008_120000_TAI.Mr_polfil.fits \
  --output-dir cases/20241008/simulations/rbsl08 \
  --extent-arcsec -400 900 -300 600 \
  --lon-range -7 65 \
  --lat-range 0 40 \
  --resolution-degree 0.05 \
  --nghost 2 \
  --nr 400 \
  --length-r 0.8
```

For a quick dry run of extraction and file formatting, add `--no-preprocess` or lower
`--preprocess-max-iter`.

## AIA To Spherical VTS

`aia_to_spherical_vtk.py` samples one AIA FITS image on a Stonyhurst longitude/latitude
grid and writes a notebook-style Cartesian spherical surface `.vts` for ParaView.

Example:

```bash
SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mplconfig PYTHONPATH=src \
python workflows/simulations/spherical/aia_to_spherical_vtk.py \
  --aia /Users/xiaoxingyun/Desktop/20241008/cme/rbsl/data/aia.lev1_euv_12s.2024-10-08T230007Z.304.image_lev1.fits \
  --output cases/20241008/products/aia_spherical_vtk/aia304_230007.vts \
  --lon-range -7 65 \
  --lat-range 0 40 \
  --resolution-degree 0.2 \
  --nghost 2 \
  --radius 1.001 \
  --quicklook
```

## RBSL Path Preparation

`prepare_rbsl_path.py` prepares the flux-rope axis file `path.dat`.

Workflow:

- show an AIA submap and pick two footpoints, or pass them with `--footpoints-arcsec`
  / `--footpoints-stonyhurst`;
- mark the two footpoints on a longitude/latitude base map;
- click additional axis points on the base map, or pass existing `pts.npy` with `--points`;
- compute the 3D axis from height `Hs` and write binary `path.dat`.

`path.dat` is written as 6000 `float64` triples:

```text
lon_deg_plus_180, lat_deg, radius_rsun
```

Example without interactive clicks:

```bash
SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mplconfig PYTHONPATH=src \
python workflows/simulations/spherical/prepare_rbsl_path.py \
  --aia /Users/xiaoxingyun/Desktop/pycode/solar_package/rbsl_py/data/aia.lev1_euv_12s.2024-10-08T230007Z.304.image_lev1.fits \
  --output-dir cases/20241008/products/rbsl_path_test \
  --submap-arcsec 400 800 100 400 \
  --lon-range -9 67 \
  --lat-range 0 40 \
  --resolution-degree 0.1 \
  --nghost 2 \
  --height-mm 20 \
  --footpoints-arcsec 686 249 455 255 \
  --path-points-stonyhurst 48.89174591 19.16451568 40.0 20.0 30.28951974 20.78224867
```

For an interactive run, omit `--footpoints-*`, `--path-points-stonyhurst`, and `--points`.
Use `--base-map <fits>` to choose a different AIA or HMI FITS background for path picking.

## Multiview Path Overlay

`plot_multiview_path.py` projects a prepared RBSL `path.dat` and its two footpoints
onto Solar Orbiter FSI and STEREO/EUVI images. The projection follows the notebook
`Sun_rot_batch` convention: rotate `(lon, lat, r/Rsun)` into the observer view and
scale by that image's `rsun_obs`.

Example:

```bash
SUNPY_CONFIGDIR=.sunpy_config MPLCONFIGDIR=.mplconfig PYTHONPATH=src \
python workflows/simulations/spherical/plot_multiview_path.py \
  --path-dat cases/20241008/products/rbsl_path_171_interactive/path.dat \
  --points cases/20241008/products/rbsl_path_171_interactive/pts.npy \
  --summary cases/20241008/products/rbsl_path_171_interactive/path_summary.json \
  --fsi /Volumes/Seagate/case/20241009/data/solo/304/solo_L2_eui-fsi304-image_20241008T230020275_V01.fits \
  --euvi data/test/observations/20241008_2300_2400/stereo_a/euvi_195/20241008_230000_n4euA.fts \
  --output-dir cases/20241008/products/rbsl_path_171_interactive/multiview
```

## Reusable Code Location

Place geometry-specific helpers under:

```text
src/solarutils/sim/spherical.py
```

Keep HMI synoptic download and observation readers in the observation layer; only simulation-grid
conversion belongs here.

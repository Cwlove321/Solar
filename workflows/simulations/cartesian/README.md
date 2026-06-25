# Cartesian Simulation Workflows

Workflows for local-box simulations on `(x, y, z)` grids.

## Scope

- Uniform or stretched Cartesian grids.
- Local active-region boxes, jet boxes, and data cubes.
- Boundary-condition preparation and diagnostics for rectangular domains.

## Expected Data Products

- Scalar cubes: density, temperature, pressure, energy terms.
- Vector cubes: magnetic field, velocity, current density.
- Boundary maps: bottom-boundary magnetic field and driving velocity.

## First Workflows To Add

- Read and validate existing binary/VTK Cartesian snapshots.
- Make 2D slices and line cuts.
- Compute curl, divergence, magnetic energy, and flux through rectangular faces.
- Export quicklook PNGs and movies.
- Synthesize optically thin or absorbing EUV/radio/white-light images from AMRVAC snapshots.

## Reusable Code Location

Place geometry-specific helpers under:

```text
src/solarutils/sim/cartesian.py
```

Keep generic readers in `src/solarutils/io` when the file format is not geometry-specific.

## Synthetic Radiation

`synthesize_radiation.py` is the cleaned workflow version of
`Desktop/jet_mhd/辐射合成.ipynb`.

It expects:

- AMRVAC `.dat` snapshots readable by `simesh.amr_loader`.
- The external DDA tools directory containing `RST_cart/` and `AIA_response/`.

Example:

```bash
python workflows/simulations/cartesian/synthesize_radiation.py \
  --dat-glob '/Volumes/ZC/jet_mhd/*dat' \
  --output-dir '/Volumes/ZC/jet_mhd/radsyn' \
  --dda-root '/Users/xiaoxingyun/Desktop/pycode/solar_package/radsyn_tools/DDA_method1' \
  --start 1 --stop 180 --step 3 \
  --wavelength 171 131 \
  --save-physical \
  --plot
```

Notebook defaults are preserved for the slab bounds, detector plane, unit
conversion, log-intensity image output, and `128`-pixel side crop in PNG
quicklooks. Use `--normal 1 0 0` for the side-on line of sight tested near the
end of the notebook.

## AIA Optically Thick From `.dat`

`synthesize_aia_optically_thick.py` is the newer single-snapshot workflow for
AIA 171/211/304 optically thick synthesis. It reads AMRVAC `.dat` files directly
with `simesh.amrvac.read_uniform`; it does not use `simesh.legacy` and does not
read cached physical `.npy/.npz` cubes.

By default it uses the highest effective AMR resolution:

```text
resolution = domain_nx * 2**(levmax - 1)
```

If `--output-dir` is omitted, products are written beside the `.dat` file in
`radiation/`. Each run writes linear intensity, log intensity, total optical
depth, a PNG quicklook, and a metadata JSON file.

Basic example:

```bash
python /Users/xiaoxingyun/Desktop/solar_code/workflows/simulations/cartesian/synthesize_aia_optically_thick.py \
  /Users/xiaoxingyun/Desktop/test1/mhd0006.dat \
  --wavelength 171 \
  --normal -1 0 0 \
  --vmin 3 \
  --vmax 10 \
  --label test1_171_normal_minus_x
```

Common view directions:

```bash
--normal 0 1 0
--normal 0 -1 0
--normal -1 0 0
--normal -1 -1 -0.5
```

Axis-aligned views use fast formal integration. Oblique views use DDA ray
integration through the projected box volume.

The transfer equation is:

```text
I = integral(exp(-tau_before) * emissivity * ds) + I0 * exp(-tau_total)
```

Useful options:

```bash
--wavelength 171          # also 211 or 304
--normal -1 -1 -0.5       # line-of-sight direction
--i0 0                    # background intensity
--z-min 0.2               # discard z < 0.2 Mm before synthesis
--vmin 3 --vmax 10        # PNG log10 color limits only
--resolution 240 240 240  # override native resolution
--kappa-scale 1e10        # diagnostic absorption scaling
--emission-scale 1        # diagnostic emissivity scaling
```

The external RST/DDA helper root defaults to the local installation:

```text
/Users/xiaoxingyun/Desktop/pycode/solar_package/radsyn_tools/DDA_method1
```

Override it with:

```bash
--dda-root /path/to/DDA_method1
```

## SHARP-like CEA Patch From Full-disk HMI Vector Data

`make_cea_patch.py` builds an arbitrary CEA vector magnetic-field patch from
full-disk HMI `field`, `inclination`, `azimuth`, and `disambig` FITS segments.
It follows the SHARP CEA order of operations: form native image-plane vector
components, sample them on the target CEA grid, then transform to
`Br/Btheta/Bphi`.

Example:

```bash
PYTHONPATH=src python workflows/simulations/cartesian/make_cea_patch.py \
  --field /path/to/hmi.B_720s.TIME.field.fits \
  --inclination /path/to/hmi.B_720s.TIME.inclination.fits \
  --azimuth /path/to/hmi.B_720s.TIME.azimuth.fits \
  --disambig /path/to/hmi.B_720s.TIME.disambig.fits \
  --output-dir /path/to/output/cea_patch \
  --center-lon 35.0 \
  --center-lat -20.0 \
  --width-degree 12.0 \
  --height-degree 9.0 \
  --resolution-degree 0.03
```

Outputs include `Br.fits`, `Bt.fits`, `Bp.fits`, `cea_patch.npz`, and
`cartesian_boundary.dat/.npz`. The local Cartesian boundary convention is
`Bx = Bp`, `By = -Bt`, and `Bz = Br`.

For closer comparison with the official SHARP CEA remapping, use
`--sampling-mode sharp`. This samples the native image-plane components on a
fine CEA grid, applies a Gaussian smoothing kernel, and block-reduces to the
requested output resolution. With the usual SHARP CEA spacing, the defaults
approximate the documented `0.01 deg -> 0.03 deg` resize:

```bash
PYTHONPATH=src python workflows/simulations/cartesian/make_cea_patch.py \
  --field /path/to/hmi.B_720s.TIME.field.fits \
  --inclination /path/to/hmi.B_720s.TIME.inclination.fits \
  --azimuth /path/to/hmi.B_720s.TIME.azimuth.fits \
  --disambig /path/to/hmi.B_720s.TIME.disambig.fits \
  --output-dir /path/to/output/cea_patch_sharp_mode \
  --center-lon 49.11080820595 \
  --center-lat -16.6075001 \
  --width-degree 34.3799991978 \
  --height-degree 16.799999608 \
  --resolution-degree 0.0299999993 \
  --sampling-mode sharp \
  --oversample-resolution-degree 0.01 \
  --smooth-sigma-degree 0.01 \
  --smooth-truncate 2
```

If an official SHARP CEA patch is available, compare the generated output with:

```bash
PYTHONPATH=src python workflows/simulations/cartesian/compare_cea_patch_to_sharp.py \
  --patch-dir /path/to/output/cea_patch_sharp_mode \
  --official-dir /path/to/official/sharp/files \
  --output-dir /path/to/output/cea_patch_sharp_mode
```

The comparison writes `official_sharp_comparison_metrics.json`,
`official_vs_make_cea_patch_components.png`, and
`official_vs_make_cea_patch_scatter_strong.png`.

## Guo Cartesian Potential/NLFFF Boundaries From CEA FITS

`prepare_nlfff_boundary.py` converts an existing `Br.fits`, `Bt.fits`, and
`Bp.fits` CEA patch into the binary input files used by Guo's
`example/NLFFF_CAR_UNI` Cartesian potential and magnetofriction examples.

The local Cartesian convention is:

```text
Bx = Bp
By = -Bt
Bz = Br
```

The pixel window is specified with 0-based Python pixels before multigrid
rebinning. `x0`, `y0`, `nx`, and `ny` describe the potential/physical region,
not the ghost-padded NLFFF boundary. The script automatically adds two ghost
cells on each side after multigrid. If the corresponding source pixels lie
outside the FITS frame, they are filled with `--blank-value` (default `0 G`).

Guo's modified multigrid rule is used: `level=1` keeps the input resolution,
`level=2` does 2x2 binning, and `level=3` does 4x4 binning. The physical
`nx` and `ny` must be divisible by `2^(level-1)`. The written
`nlfff_boundary.dat` is four pixels larger than `potential_boundary.dat`
because two ghost pixels are kept on each side.

Example using the 2024-10-08 CEA patch:

```bash
PYTHONPATH=src python workflows/simulations/cartesian/prepare_nlfff_boundary.py \
  --cea-dir /Users/xiaoxingyun/Desktop/solar_code/cases/20241008/cea_patch_center40_lat20_w40_h24 \
  --output-dir /Users/xiaoxingyun/Desktop/solar_code/cases/20241008/nlfff_boundary_center40_lat20_w40_h24_x80_y0_1200x800_lv3 \
  --x0 80 --y0 0 --nx 1200 --ny 800 --level 3 \
  --no-preprocess --vmax 500
```

Main outputs:

```text
potential_boundary.dat  # cropped Bz-only binary for potential extrapolation
nlfff_boundary.dat      # full Bx, By, Bz binary including ghost pixels
allboundaries.dat       # Guo multigrid ASCII intermediate
grid.ini
nd.ini
amrvac_parameters.txt   # domain_nx and xprobmin/xprobmax snippet
parameters.json
boundary_quicklook.png
```

Add `--preprocess` to apply the Wiegelmann-style force/torque preprocessing
after multigrid binning. The defaults match Guo's IDL wrapper:
`mu3=0.001`, `mu4=0.01`, and `max_iter=2000`.

# Simulation Analysis

This workflow collects the AMRVAC analysis tasks that recur across PFSS,
relaxation, RBSL insertion, and full-MHD runs. It is named
`Simulation_Analysis` rather than `Simulation_analysis_demo` because the code
is intended for routine production use.

The workflow grew out of the diagnostics used for `rbsl12`, `rbsl13`, and
`rbsl14`. Case paths, hard-coded output names, and copied physical constants
have been replaced by one JSON configuration.

## Contents

```text
Simulation_Analysis/
|-- README.md
|-- analysis.example.json
|-- simulation_analysis/
|   |-- config.py       # defaults and JSON configuration
|   |-- geometry.py     # spherical AMR and stretched-r geometry
|   |-- physics.py      # state interpretation and derived variables
|   |-- snapshot.py     # AMR slices and extrema
|   |-- rbsl_balance.py # RBSL axis sampling and balance current
|   |-- logs.py         # AMRVAC per-step log parser
|   `-- plotting.py     # shared plotting styles
|-- scripts/
|   |-- inspect_snapshot.py
|   |-- plot_slice.py
|   |-- plot_slice_series.py
|   |-- compare_slices.py
|   |-- mark_extremum_on_bottom_br.py
|   |-- find_magnetic_nulls.py
|   |-- monitor_run.py
|   |-- plot_log_history.py
|   |-- integrate_energy_budget.py
|   |-- calculate_rbsl_balance.py
|   |-- calculate_null_y_qsl_slices.py
|   |-- calculate_fixed_y_qsl_slices.py
|   |-- plot_y0_enclosed_cross_section.py
|   |-- q4_seam_tracking.py
|   |-- make_y0_hybrid_mask_evolution.py
|   `-- make_movie.py
|-- examples/
|   `-- validate_cold_dense_structure.py
`-- tests/
    `-- test_core.py
```

## Requirements

The scripts require Python 3.10 or newer, NumPy, SciPy, Matplotlib, Pillow,
scikit-image, and the current `simesh` AMRVAC reader. QSL calculations also
require the repository copy of `Packages/FastQSL2`.

On this Mac, use the base interpreter directly:

```bash
cd /Users/xiaoxingyun/Desktop/solar_code
/opt/anaconda3/bin/python workflows/Simulation_Analysis/scripts/inspect_snapshot.py --help
```

The scripts automatically find the bundled source at
`Packages/simesh/src`. On `guoguo`, use the default `zcc` Python:

```bash
/share/home/guoyang1/apps/miniconda3/envs/zcc/bin/python script.py ...
```

If `simesh` is neither installed nor in the same repository, point to its
source explicitly:

```bash
export SIMESH_SRC=/path/to/simesh/src
```

`calculate_null_y_qsl_slices.py` resolves `Packages/FastQSL2` relative to the
repository. Its magnetic input, null catalogue, and output directories default
to a sibling `jet_mhd` checkout and can be overridden either with command-line
arguments or the `FASTQSL_MAGNETIC_DIR`, `FASTQSL_NULL_DIR`,
`FASTQSL_OUTPUT_DIR`, and `FASTQSL2_DIR` environment variables.

## Quick Start

1. Copy `analysis.example.json` into the simulation directory as
   `analysis.json`.
2. Change `run_dir`, grid parameters, units, and slice position.
3. Run one snapshot inspection before making plots.

```bash
WORKFLOW=/Users/xiaoxingyun/Desktop/solar_code/workflows/Simulation_Analysis
PY=/opt/anaconda3/bin/python

$PY "$WORKFLOW/scripts/inspect_snapshot.py" \
  --config /path/to/run/analysis.json \
  --variable temperature_K --extreme max
```

Every command writes into `output_dir`, which is relative to `run_dir` unless
an absolute path is given.

## Configuration

The complete template is [analysis.example.json](analysis.example.json). The
important sections are:

### Run selection

```json
{
  "run_dir": "/share/home/.../one_run",
  "data_glob": "data/*_[0-9][0-9][0-9][0-9].dat",
  "output_dir": "python/analysis"
}
```

`data_glob` is a Python glob relative to `run_dir`. Four bracket expressions
are used because Python glob does not understand regular-expression
quantifiers such as `{4}`.

### Radial grid

For the mature stretched grid:

```json
{
  "grid": {
    "radial_type": "geometric",
    "qstretch_base": 1.041541581503301,
    "surface_radius_code": 1.0,
    "radius_unit_Mm": 695.5
  }
}
```

At AMR level `l`, the refinement factor is

\[
f_l=2^{l-1},\qquad q_l=q_1^{1/f_l}.
\]

For `N_r` radial cells at that level,

\[
\Delta r_0=(r_{max}-r_{min})\frac{q_l-1}{q_l^{N_r}-1},
\]

and edge `i` is

\[
r_i=r_{min}+\Delta r_0\frac{q_l^i-1}{q_l-1}.
\]

This must match the stretch used by `mod_usr.t`. Set `radial_type` to
`uniform` for an unstretched run. A wrong `qstretch_base` shifts cells,
corrupts volumes, and makes current/divergence derivatives physically wrong.

### State interpretation

```json
{
  "physics": {
    "energy_kind": "internal",
    "momentum_model": "semirelativistic",
    "reduced_c_km_s": 6000.0,
    "eta_code": 2.4693930682273452e-6,
    "gamma": null
  }
}
```

- `energy_kind=internal`: DAT field `e` is internal-energy density.
- `energy_kind=total`: kinetic and magnetic energies are removed before
  pressure is calculated.
- `momentum_model=classical`: stored momentum is `rho*v`.
- `momentum_model=semirelativistic`: stored momentum contains the Boris
  correction and is inverted before velocity diagnostics are calculated.
- `gamma=null`: read `gamma` from DAT metadata.
- `eta_code`: explicit magnetic diffusivity used only by the energy-budget
  diagnostic. It must match the run.

For the paper value `eta=200 km^2/s` in the current normalization,
`eta_code=2.4693930682273452e-6`. The former high-resistivity test used about
`1e-4`, or roughly `8100 km^2/s`.

### Units

All unit values describe one code unit. The supplied defaults match the
current spherical full-MHD workflow:

```json
{
  "units": {
    "length_cm": 6.955e10,
    "number_density_cm3": 1e9,
    "temperature_K": 1e6,
    "velocity_km_s": 116.448846777562,
    "magnetic_field_G": 1.9977023903979652,
    "time_s": 5972.4772,
    "pressure_erg_cm3": 0.31754927
  }
}
```

Do not reuse these blindly for a case with different normalization.

### Slice position

```json
{
  "slice": {
    "axis": "phi",
    "coordinate": 176.5,
    "angular_unit": "degree",
    "target_level": null,
    "max_height_Mm": 60.0,
    "angular_bounds_degree": [8.0, 22.0],
    "reference_angle_degree": 15.0,
    "equal_physical_scale": true
  }
}
```

- `axis` is the fixed coordinate: `r`, `theta`, or `phi`.
- `target_level=null` means the highest leaf level present in that DAT.
- A lower target level is rejected so that fine cells are not silently lost.
- For fixed `phi`, `angular_bounds_degree` means latitude bounds.
- For fixed `theta`, it means longitude bounds.
- `reference_angle_degree` is the zero of horizontal arc distance.
- `equal_physical_scale=true` gives one plotted Mm the same length in both
  directions. This is the correct choice for judging shape and inclination.

## Commands

### 1. Inspect one DAT

```bash
$PY "$WORKFLOW/scripts/inspect_snapshot.py" \
  --config analysis.json \
  --dat data/run_0006.dat \
  --variable alfven_speed_km_s --extreme max
```

The JSON report includes:

- domain, block size, field names, and leaf counts by AMR level;
- ranges of density, temperature, pressure, magnetic field, velocity,
  sound speed, Alfven speed, and beta;
- the extremum's block, AMR level, local cell index, radius, height,
  latitude, longitude, and local physical context.

Useful variables include `temperature_K`, `number_density_cm3`, `br_G`,
`bmag_G`, `vmag_km_s`, `sound_speed_km_s`, `alfven_speed_km_s`, `beta`, and
`amr_level`.

For semi-relativistic runs, `alfven_speed_km_s` is the classical diagnostic
`B/sqrt(rho)`, while `alfven_reduced_km_s` is the Boris-limited value that
approaches `reduced_c_km_s` in very strong field.

### 2. Plot one fixed-coordinate slice

```bash
$PY "$WORKFLOW/scripts/plot_slice.py" \
  --config analysis.json \
  --dat data/run_0006.dat \
  --axis phi --coordinate 176.5 \
  --fields temperature_K number_density_cm3 br_G amr_level \
  --fieldlines --columns 2 --save-npz
```

The fixed `phi` and fixed `theta` displays use local surface arc distance in
Mm and height in Mm. Magnetic streamlines use the field projected into the
same physical plane:

- fixed `phi`: horizontal component `-B_theta`, vertical component `B_r`;
- fixed `theta`: horizontal component `B_phi`, vertical component `B_r`;
- fixed `r`: horizontal component `B_phi`, vertical component `-B_theta`.

The white/blue curves are projected field lines, not individual 3D field
lines and not plasma trajectories.

Scalar cells are rendered on their actual non-uniform radial coordinates with
`pcolormesh`. Only the projected vector field is interpolated to a uniform
physical plotting grid because Matplotlib `streamplot` requires equally spaced
axes. Thus the stretched-grid cell widths are not visually flattened.

### 3. Plot a time series

Six rows and one column, matching true physical proportions:

```bash
$PY "$WORKFLOW/scripts/plot_slice_series.py" \
  --config analysis.json \
  --axis theta --coordinate 75.15 \
  --first 0 --count 6 --stride 1 --layout 6x1 \
  --fields temperature_K number_density_cm3 \
  --fieldlines
```

One PNG is written per requested field. All times share the same color limits
from the config, so apparent evolution is not caused by automatic rescaling.

### 4. Compare two snapshots

```bash
$PY "$WORKFLOW/scripts/compare_slices.py" \
  --config analysis.json \
  --dat-a data/run_0000.dat --dat-b data/run_0001.dat \
  --axis phi --coordinate 176.5 --field br_G \
  --difference-limit 100
```

The first two panels share one normalization and the third shows `B-A` with a
symmetric color scale. Add `--relative` for `(B-A)/abs(A)`.

### 5. Mark a failure or extreme on bottom Br

```bash
$PY "$WORKFLOW/scripts/mark_extremum_on_bottom_br.py" \
  --config analysis.json \
  --variable temperature_K --extreme max --vmax 3000
```

The script first finds the true 3D AMR leaf-cell extremum, then marks its
longitude and latitude on the first active radial layer. The marker label also
records height and AMR level.

### 6. Monitor a running job

One-time report:

```bash
$PY "$WORKFLOW/scripts/monitor_run.py" --config analysis.json --job-id 282398
```

Continuous report every ten minutes:

```bash
nohup $PY "$WORKFLOW/scripts/monitor_run.py" \
  --config analysis.json --job-id 282398 --watch --interval 600 \
  > python/analysis/monitor.out 2>&1 &
```

Outputs:

- `run_monitor.json`: job state, current iteration/time/dt, DAT list,
  minimum positive dt, first non-finite row, non-empty stderr tails, and
  recent slopes;
- `run_monitor.csv`: all parsed per-step rows.

The default parser matches the current custom AMRVAC line with 19 columns.
If `speciallog` changes, edit `log.columns` in the JSON instead of editing
Python.

### 7. Plot log history

```bash
$PY "$WORKFLOW/scripts/plot_log_history.py" --config analysis.json
```

This makes the three-panel `dt`, momentum-extrema, and density/internal-energy
history used to judge collapse and relaxation plateaus.

### 8. Integrate the energy budget

```bash
$PY "$WORKFLOW/scripts/integrate_energy_budget.py" \
  --config analysis.json --first 0 --count 7 --stride 1
```

Outputs are CSV, JSON, and a three-panel plot containing:

\[
E_{int}=\int e_{int}\,dV,
\qquad
E_B=\int \frac{B^2}{2}\,dV,
\qquad
E_K=\int \frac{\rho v^2}{2}\,dV,
\]

\[
P_J=\int \eta J^2\,dV,
\qquad
P_C=\int -p\nabla\cdot\mathbf{v}\,dV.
\]

The exact spherical cell volume is

\[
\Delta V=\frac{r_{i+1}^3-r_i^3}{3}
(\cos\theta_j-\cos\theta_{j+1})
(\phi_{k+1}-\phi_k).
\]

Current and velocity divergence use the full spherical expressions and
second-order differences on each leaf block with two AMR-exchanged ghost
cells.

Important interpretation:

- `eta*J^2` is a diagnostic heating rate. It contributes to internal energy
  only if the simulation's energy equation actually includes Joule heating.
- The residual
  `Delta E_int - integral(P_J + P_C) dt` is not automatically numerical
  error. It also contains conductive flux, radiation, background heating,
  energy advection, and boundary fluxes.
- DAT files do not store physical-boundary ghost cells. Derivatives in the
  first and last active cells depend on the reconstruction used by `simesh`.
  Treat boundary-local source values more cautiously than volume trends.

### 9. Make a movie

```bash
$PY "$WORKFLOW/scripts/make_movie.py" \
  --frames 'python/analysis/frames/*.png' \
  --output python/analysis/evolution.mp4 --fps 6
```

This requires `ffmpeg` on `PATH` and writes H.264/YUV420p for broad playback
compatibility.

### 10. Find magnetic nulls on a rectilinear grid

First save the three magnetic components on one common rectilinear grid. Array
order must be `(nx, ny, nz)`. Physical coordinate arrays are strongly
recommended, especially when the three cell widths differ:

```python
np.savez_compressed(
    "magnetic_field.npz",
    Bx=bx,
    By=by,
    Bz=bz,
    x=x_coordinates,
    y=y_coordinates,
    z=z_coordinates,
)
```

Run the null search with:

```bash
$PY "$WORKFLOW/scripts/find_magnetic_nulls.py" magnetic_field.npz \
  --output python/analysis/magnetic_nulls.json --save-npz
```

If the archive does not contain `x`, `y`, and `z`, provide uniform spacing and
the coordinate of the first grid point explicitly:

```bash
$PY "$WORKFLOW/scripts/find_magnetic_nulls.py" magnetic_field.npz \
  --origin 0 0 0 --spacing 0.0234375 0.0234375 0.0234375
```

The QSL/FastQSL binary format can be read directly without first duplicating
the magnetic arrays in an NPZ. Its default on-disk array order is `(z,y,x)`:

```bash
$PY "$WORKFLOW/scripts/find_magnetic_nulls.py" nlfff_0010.binary \
  --origin -39.53125 -59.53125 0.46875 \
  --spacing 0.9375 0.9375 0.9375 \
  --output magnetic_nulls_t0010.json --save-npz
```

Here the origin values are cell-center coordinates in Mm for a domain whose
displayed cell-boundary extent is `[-40,80] x [-60,60] x [0,120]` Mm.

The default search applies the eight-corner component-sign condition and then
searches every compatible cell. It also records a solid-angle Poincare degree
computed on a twice-subdivided cell surface, but does not use that diagnostic
as a hard filter. It uses several deterministic starting points for a damped
Newton solve of the trilinear field.
Newton updates are performed in local cell coordinates, while the reported
Jacobian is converted to the supplied physical coordinates. Every accepted
point must pass a magnetic residual test and, by default, have a non-singular
physical-coordinate Jacobian. Add `--keep-degenerate` only when deliberately
investigating non-isolated null lines or surfaces.

The JSON output records both physical and continuous index positions,
residuals, Poincare degree, Jacobian, divergence diagnostic, eigenvalues,
eigenvectors, and real spine/fan basis vectors. Type names state the signs of
the fan and spine eigenvalues explicitly instead of relying on ambiguous
`A`/`B` conventions.

Important limitations:

- Input arrays must represent one common rectilinear grid. Interpolate AMR
  leaf data conservatively to such a grid before using this script; do not
  treat duplicated coarse cells as independent high-resolution samples.
- The optional `--candidate-method degree` filter assumes an isolated,
  non-degenerate null with non-zero net degree in each cell. The default
  `sign` method is safer for a coarse field or a cell containing a close pair.
  Add `--seed-grid 3` for a more exhaustive multi-start check, but refine the
  grid before interpreting any closely spaced null pair.
- A null on a computational boundary is not an isolated volume null. Use
  `--boundary-cells 1` to exclude the outermost cell layer when the boundary
  values are numerical padding or otherwise unsuitable for topology.
- A large `relative_divergence` in the output means the local field is not
  sufficiently solenoidal for reliable spine/fan interpretation.

### 11. Track the y=0 QSL-bounded cross-section

The production command below reuses exact `y=0`, `delta=0.25` FastQSL2
slices for snapshots 10--95 and calculates only missing slices. It then makes
the white-background mask animation with the purple `rboundary==11` outer
boundary and the orange `log10(Q)=4` boundary:

```bash
$PY "$WORKFLOW/scripts/make_y0_hybrid_mask_evolution.py"
```

The first and last masks are anchored by the measured magnetic-axis
intersections. Forward and reverse component tracking automatically identify
the interval where the left and right low-Q regions merge. Within that
interval, `q4_seam_tracking.py` follows the same bottom-rooted Q=4 crossing
with a row-wise dynamic-programming path. This avoids switching to the large
right-hand region when the discrete low-Q components merge and split. The
mask is filled from the current purple left boundary to the tracked QSL and
all enclosed holes are made solid.

The output directory contains:

- one GIF and all 86 rendered PNG frames;
- one compact mask-stack NPZ with `masks`, `snapshots`, `time_minutes`,
  `x_Mm`, and `height_Mm`;
- one NPZ per snapshot with the binary mask, its right boundary, the tracked
  Q4 seam where applicable, and every selected seed in both plotting and
  AMRVAC coordinates;
- CSV and JSON quality reports with areas, component/hole checks, temporal
  IoU, seam gaps, and endpoint-regression results.

For field-line tracing, read `selected_seed_xyz_amrvac`. Its columns are
`(x, y, z_absolute)` and use
`x=x_Mm/10`, `y=0`, and `z_absolute=height_Mm/10+0.1`. The accompanying
`selected_seed_grid_indices_iz_ix`, `selected_seed_log10Q`, and
`selected_seed_rboundary` arrays preserve the grid-level provenance of every
seed.

### 12. Trace the y=0 mask to bottom footpoints and measure flux

After generating the 86 solid masks above and mounting the Seagate data disk,
trace their field lines to the bottom boundary with:

```bash
$PY "$WORKFLOW/scripts/trace_y0_mask_footpoints_flux_evolution.py"
```

The left panel of the resulting GIF reproduces the selected y=0 cross-section.
The right panel is independently selected on a bottom grid with four samples
per native magnetic cell: a bottom sample belongs to the footprint when its
field line crosses the solid y=0 mask at least once. A sample is counted only
once even if its path crosses the plane repeatedly. The paired footprint
requires fresh `rboundary==11` tracing and opposite signs of `Bz` at its two
terminals. Zero-length bottom traces and pixels that reach another boundary
are saved as unpaired diagnostics and are excluded from the paired flux.

The two requested flux measurements are kept distinct:

\[
\Phi_y=\int_{A_{y=0}} B_y\,dx\,dz,
\qquad
\Phi_{z,\mathrm{pair}}=
\frac{\Phi_{z,+}+|\Phi_{z,-}|}{2}.
\]

The CSV also contains signed, positive, negative, and unsigned components for
the full solid y=0 mask, its bottom-to-bottom subset, and the non-closed filled
pixels. Units are maxwells with `B_unit_G=2` and `length_unit_Mm=10`. In this
cell-centred FastQSL cube, traced bottom terminals lie at `z=0.146875`; `z=0.1`
is retained separately as the AMRVAC bottom-height datum.

Default outputs are under
`qfactor/y0_qsl/enclosed_region/evolution_logq4/bottom_footpoint_flux/`:

- `y0_to_bottom_footpoints_t0010_t0095.gif`: white-background evolution GIF;
- `y0_mask_bottom_flux_t0010_t0095.csv`: all per-frame fluxes and QA counts;
- `y0_mask_bottom_footpoints_t0010_t0095_stack.npz`: reusable 86-frame source
  masks, positive/negative bottom masks, crossing counts, bottom `Bz`, and flux
  arrays;
- `data/traced_y0_mask_footpoints_tXXXX.npz`: complete endpoints and per-frame
  provenance;
- `bottom_flux_convergence_t0010_t0095.csv`: 2x, 4x, and 8x bottom-sampling
  checks at snapshots 10, 25, 26, and 95;
- `y0_mask_bottom_footpoints_t0010_t0095_metadata.json`: definitions, units,
  settings, quality summaries, and source paths.

To draw a clean bottom-boundary-only flux evolution figure from the saved CSV:

```bash
$PY "$WORKFLOW/scripts/plot_bottom_footpoint_flux_evolution.py"
```

This writes `bottom_boundary_bz_flux_evolution_t0010_t0095.png` and PDF with
the positive flux, absolute negative flux, and their half unsigned sum.

The quantitative mask is never morphologically closed or filled for display.
If a selected footprint touches the trace-window edge, the script stops rather
than silently clipping its flux.

### 13. Animate the exact-y=0 curl-to-field ratio

Use the same uniform magnetic binaries as the QSL calculation to evaluate

\[
\frac{|\nabla\times\mathbf{B}|}{|\mathbf{B}|}
\]

for snapshots 10--95:

```bash
$PY "$WORKFLOW/scripts/make_y0_j_over_b_evolution.py"
```

The native magnetic grid is cell centred and has an even number of y cells,
so `y=0` lies halfway between `y=-0.046875` and `y=+0.046875`. The script
linearly interpolates B to that exact plane, uses the centred difference
between the two layers for y derivatives, and differentiates the interpolated
plane in x and z. This avoids treating index 63 as y=0.

The animation uses the continuous Matplotlib `rainbow` colormap with one
logarithmic color range derived from all 86 frames and an `x=-15..45 Mm`,
`height=0..60 Mm` view.

To compare the first and last selected times side by side with a shared y axis
and color bar, while overlaying exact-y=0 `log10(Q)=1,2,4` contours:

```bash
$PY "$WORKFLOW/scripts/make_y0_j_over_b_comparison.py"
```

The contour colors and solid/dashed/dash-dot line styles are repeated in the
legend, and a white line halo keeps them readable over the rainbow background.
The same figure marks the local `J/B` maximum in `x=12..18 Mm` and
`height=8..60 Mm` as an approximate current-sheet location. Adjust that
diagnostic window with `--peak-window XMIN XMAX HMIN HMAX` when needed.

The saved data remain on the complete native 128-by-128 x-z plane. Because
MPI-AMRVAC uses normalized `J=curl(B)`, the ratio is stored in both `L0^-1`
and `Mm^-1`, with `L0=10 Mm`; it is not an SI current-density ratio.

Default outputs are under `qfactor/y0_qsl/j_over_b_evolution/`:

- `y0_j_over_b_t0010_t0095.gif`: 86-frame evolution animation;
- `y0_j_over_b_t0010_t0095_stack.npz`: snapshots, times, coordinates,
  B and curl(B) vectors and magnitudes, both ratio units, and validity masks;
- `data/y0_j_over_b_tXXXX.npz`: independently reusable per-snapshot slices;
- `frames/y0_j_over_b_tXXXX.png`: rendered source frames;
- `y0_j_over_b_t0010_t0095_statistics.csv`: full-domain and displayed-window
  percentiles for every snapshot;
- `y0_j_over_b_t0010_t0095_metadata.json`: definitions, interpolation and
  derivative conventions, source paths, display limits, checks, and hashes.

### 14. Calculate the RBSL balance current

Use a background snapshot before the flux rope is inserted. The path file may
be the original binary float64 file or the mature code's ASCII `lon, lat, r`
file. `path_summary.json` supplies the fitted depth and major radius.

For a uniform radial mesh, a DAT can be sampled directly:

```bash
$PY "$WORKFLOW/scripts/calculate_rbsl_balance.py" \
  --config analysis.json \
  --dat data/background_0006.dat \
  --path-dat path.dat \
  --path-summary path_summary.json \
  --minor-radius-mm 15 \
  --current-flux-mx 1.29e20 \
  --helicity negative
```

For a geometrically stretched radial mesh, first let AMRVAC write the true
coordinates and then use `--vtu`:

```bash
mpirun ./amrvac -i amrvac.par -if data/background_0006.dat -convert

$PY "$WORKFLOW/scripts/calculate_rbsl_balance.py" \
  --config analysis.json \
  --vtu data/background_0006.vtu \
  --path-dat path.dat \
  --path-summary path_summary.json \
  --minor-radius-mm 15 \
  --current-flux-mx 1.29e20 \
  --helicity negative
```

The VTU route expects AMRVAC raw-binary point output (`vtuB` or `vtuBmpi`)
with the default `nocartesian=.false.`. It reads the converter's Cartesian
coordinates and vector components directly. Do not use direct DAT balance
sampling for `grid.radial_type="geometric"`: the current `simesh` AMR mesh
reconstruction assumes uniform block spacing and can associate the path with
the wrong physical radius.

At every path sample, the script locates the AMR leaf block and performs
trilinear interpolation of all three magnetic components. The DAT route uses
two exchanged ghost layers; the VTU route uses AMRVAC's true block-corner
coordinates. It projects the external field onto the local rope binormal and
evaluates

\[
I_{balance}(s)=
\frac{R_c B_{strap}(s)}{\ln(8R_c/a)-1},
\qquad
I_{RBSL}=\pm\frac{5\sqrt{2}\,F}{12\pi a}.
\]

The reported global value is the arc-length-weighted mean over the central
10%--90% of the axis. The output directory contains a JSON summary, the full
axis profile as CSV, and a three-panel plot of `B_strap`, normalized current,
and equivalent local flux. `--current-flux-mx` is optional; when supplied it
adds the trial current and its ratio to the inferred balance current.

This is a thin-torus diagnostic, not a guarantee of pointwise equilibrium.
If `B_strap` changes sign or varies strongly along the axis, no single current
can balance the entire rope. Also keep the background self-consistent: after
changing `F_flx`, regenerate the RBSL-subtracted boundary/PFSS, relax it, and
repeat this calculation when the subtraction materially changes the field.

## Cold/Dense Structure Validation

The low-temperature, high-density tongue analysis is kept as an example
because its physical criteria are project-specific:

```bash
$PY "$WORKFLOW/examples/validate_cold_dense_structure.py" \
  --config analysis.json --axis phi --coordinate 176.5 --all \
  --temperature-anomaly -0.10 --density-anomaly 0.10 \
  --min-height 1 --max-height 40
```

At each height, define horizontal anomalies relative to the median:

\[
\delta\log_{10}T=\log_{10}T-
\operatorname{median}_{horizontal}(\log_{10}T),
\]

and similarly for density. The default structure mask is

\[
\delta\log_{10}T\le -0.10,
\qquad
\delta\log_{10}n\ge +0.10.
\]

For every DAT it plots temperature, density, entropy anomaly, radial
velocity, the acute angle between velocity and magnetic field, and AMR level.
The metrics include the fraction of structure pixels near an AMR interface.
This separates three hypotheses:

1. field-guided material redistribution;
2. a thermodynamic source-term artifact;
3. an AMR coarse/fine-interface artifact.

## Implemented Physics

### Pressure and temperature

For an internal-energy DAT:

\[
p=(\gamma-1)e_{int},
\qquad
T=\frac{p}{\rho}T_{unit}.
\]

For a total-energy DAT:

\[
e_{int}=E-\frac{1}{2}\rho v^2-\frac{1}{2}B^2.
\]

The formula assumes the AMRVAC normalization used by the supplied units. If a
model stores pressure directly or uses a different mean molecular weight,
update `physics.py` or the field mapping explicitly.

### Semi-relativistic momentum inversion

For the current Boris-corrected model, split stored momentum into components
parallel and perpendicular to `B`:

\[
\mathbf{m}_{parallel}=
\frac{\mathbf{m}\cdot\mathbf{B}}{B^2}\mathbf{B}.
\]

The velocity is reconstructed as

\[
\mathbf{v}=
\frac{\mathbf{m}_{parallel}}{\rho}+
\frac{\mathbf{m}-\mathbf{m}_{parallel}}
{\rho+B^2/c_{reduced}^2}.
\]

Using `m/rho` for this model overestimates perpendicular velocity in strong
field regions. Set `momentum_model=classical` for ordinary MHD.

### Speeds and beta

\[
c_s=\sqrt{\gamma p/\rho},
\qquad
v_A=\frac{|B|}{\sqrt{\rho}},
\qquad
c_{f,upper}=\sqrt{c_s^2+v_A^2},
\qquad
\beta=\frac{2p}{B^2}.
\]

`c_f,upper` is an angle-independent upper estimate, not the exact directional
fast-mode speed.

For a Boris-reduced model, the reported reduced Alfven diagnostic is

\[
v_{A,reduced}=\frac{v_A}{\sqrt{1+v_A^2/c_{reduced}^2}}.
\]

## Lessons Encoded in the Workflow

1. Always inspect the first DAT before plotting a series. Wrong energy or
   momentum interpretation can produce plausible-looking but false maps.
2. A stretched radial grid cannot be analyzed as `linspace(rmin,rmax,nr)`.
   This was a major source of misleading heights and derivatives in early
   tests.
3. Plot fixed-theta/fixed-phi slices in Mm with equal aspect when judging
   geometry. Degree-versus-Mm panels distort angles and widths.
4. Use common color limits across time. Independent auto-scaling can create a
   false impression of heating, cooling, or disappearance.
5. Dynamic AMR must be shown alongside suspicious structures. A feature that
   tracks a level transition deserves a numerical-control run.
6. A smaller CFL number delays a force imbalance but does not remove it.
   Track location, density, field strength, sound speed, and Alfven speed at
   the minimum-dt region.
7. Semi-relativistic correction limits magnetic characteristic speeds; it
   does not prevent real velocity generated by an inconsistent Lorentz force.
8. Optical-thin radiation is local. It removes internal energy from the cell
   through `n_e n_H Lambda(T)` and does not transport emitted radiation into
   neighboring cells as RMHD would.
9. Relaxation without radiation can be valid when the goal is a stable
   magnetic/atmospheric initial state. Radiation may be enabled after rope
   insertion, but its heating balance and turn-on transient must be tested.
10. High explicit resistivity accelerates magnetic diffusion and Joule
    conversion. Use the paper value for production unless a deliberately
    stronger relaxation experiment is being documented.
11. A momentum plateau is not by itself convergence. Check kinetic energy,
    magnetic energy, internal energy, spatial extrema, and whether trends are
    flat over several output intervals.
12. DAT does not contain boundary ghost values. Claims about the first active
    layer must distinguish stored data from reconstructed ghosts.
13. `simesh.mesh.rnode` currently reconstructs AMR block radii uniformly. For
    a geometrically stretched radial mesh, use AMRVAC VTU coordinates for
    point sampling such as RBSL balance-current estimates.

## Adding a New Diagnostic

New scripts should follow four rules:

1. Read paths and units from `analysis.json`.
2. Reuse `extract_axis_slice`, `derive_state`, and spherical geometry instead
   of copying constants or interpolation formulas.
3. Write both a figure and machine-readable JSON/CSV metrics.
4. Record whether values come from native leaf cells, repeated coarse cells,
   interpolation, or reconstructed ghost cells.

Case-specific science criteria belong in `examples/`; reusable I/O,
coordinate, and MHD formulas belong in `simulation_analysis/`.

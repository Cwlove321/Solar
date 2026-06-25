# Spherical RBSL AMRVAC Templates

This directory contains the baseline MPI-AMRVAC setup files for a spherical RBSL
simulation. The observational boundary data are generated on the Python side and
then copied into the expected AMRVAC input locations.

## Directory Layout

```text
sphere_rbsl/
├── pfss/            # PFSS background field and RBSL flux-rope insertion
├── mf/              # Magnetofrictional relaxation toward an NLFFF state
├── mhd_0beta/       # Zero-beta MHD evolution from the relaxed field
└── mhd_internal_e/  # MHD evolution with internal energy and filament loading
```

Each stage contains an `amrvac.par`, a `mod_usr.t`, a batch `job` script, and a
`data/` output directory. The stages are intended to run in sequence because
later stages restart from earlier-stage outputs.

## Python-Generated Inputs

The Python workflow currently writes the boundary products with
`workflows/simulations/spherical/prepare_rbsl_boundary.py`.

Required files:

```text
pfss/potential_boundary/pfss_boundary.dat
mf/nlfff_boundary/cutout.dat
pfss/rbsl/path.dat
```

`pfss_boundary.dat` is the global radial lower boundary. It is read by
`pfss/mod_usr.t`, which subtracts the RBSL radial contribution from the selected
patch and writes:

```text
pfss/potential_boundary/pfss_boundary1.dat
```

`cutout.dat` is the local vector magnetogram boundary. It is read by:

```text
mf/mod_usr.t
mhd_0beta/mod_usr.t
mhd_internal_e/mod_usr.t
```

The binary layout of `cutout.dat` is:

```text
int32 nx
int32 ny
float64 Br[nx * ny]
float64 Btheta[nx * ny]
float64 Bphi[nx * ny]
```

`path.dat` is the flux-rope/filament axis file. The Python-side generation
method will be added later. The AMRVAC side currently expects a binary stream
file with exactly 6000 points:

```text
float64 lon_deg
float64 lat_deg
float64 radius_rsun
```

`pfss/mod_usr.t` reads it as `./rbsl/path.dat` for RBSL insertion.
`mhd_internal_e/mod_usr.t` reads it as `../pfss/rbsl/path.dat` to load cool,
dense filament material around the same axis.

## Parameter Synchronization

After generating Python products, copy the values from `parameters.txt` or
`parameters.json` into the AMRVAC files.

Patch indices in `pfss/mod_usr.t`:

```fortran
lat0 = ...
lat1 = ...
lon0 = ...
lon1 = ...
```

Mesh and physical domain in every active `amrvac.par`:

```fortran
domain_nx1 = ...
domain_nx2 = ...
domain_nx3 = ...
xprobmin1 = ...
xprobmax1 = ...
xprobmin2 = ...
xprobmax2 = ...
xprobmin3 = ...
xprobmax3 = ...
```

For the current 2024-10-08, 23:00, `res=0.1 deg`, `lon=[13,61]`,
`lat=[8,32]` product, the Python output reports:

```text
lat0 = 580
lat1 = 819
lon0 = 1930
lon1 = 2409

domain_nx1 = 240
domain_nx2 = 240
domain_nx3 = 480
xprobmin1 = 1.001
xprobmax1 = 1.801
xprobmin2 = 0.16111111111111112
xprobmax2 = 0.22777777777777777
xprobmin3 = 0.5361111111111111
xprobmax3 = 0.6694444444444444
```

Some template files may still contain exploratory or older values. Treat the
Python-generated `parameters.json` as the source of truth for a given run.

## Run Order

1. Generate Python-side boundary files:

   ```bash
   PYTHONPATH=src python workflows/simulations/spherical/prepare_rbsl_boundary.py ...
   ```

2. Copy the generated files into the AMRVAC template:

   ```bash
   cp <python-output>/pfss_boundary.dat pfss/potential_boundary/pfss_boundary.dat
   cp <python-output>/cutout.dat mf/nlfff_boundary/cutout.dat
   cp <python-output>/path.dat pfss/rbsl/path.dat
   ```

3. Build and run `pfss/`.

   `pfss` reads `potential_boundary/pfss_boundary.dat` and `rbsl/path.dat`,
   writes `potential_boundary/pfss_boundary1.dat`, computes harmonic
   coefficients, and produces a restart in `pfss/data/`.

4. Build and run `mf/`.

   `mf` restarts from the PFSS output, reads `mf/nlfff_boundary/cutout.dat`,
   and relaxes the magnetic field with magnetofriction.

5. Build and run either `mhd_0beta/` or `mhd_internal_e/`.

   `mhd_0beta` is the zero-beta dynamic evolution. `mhd_internal_e` includes
   internal energy, gravity, and filament-density loading around `path.dat`.

## Notes

- The coordinate system is spherical in all stages.
- Magnetic-field normalization is defined in each `mod_usr.t`; the boundary
  files are written in Gauss by the Python workflow and normalized after read.
- `cutout.dat` is the standard local-boundary filename. Older names such as
  `04field0_cutout.dat` should not be used in new runs.
- The batch `job` files are LSF examples and may need queue, MPI count, and log
  directory edits on a different cluster.

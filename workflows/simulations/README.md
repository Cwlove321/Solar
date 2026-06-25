# Simulation Workflows

Reusable simulation workflows belong here.

Simulation workflows are split by grid geometry:

- [cartesian](cartesian/README.md): local-box simulations in `(x, y, z)`.
- [spherical](spherical/README.md): global or shell simulations in `(r, theta, phi)` or longitude/latitude grids.

## Boundary

- Put reusable simulation readers, grid helpers, and physical diagnostics in `src/solarutils`.
- Put generic command-line scripts in `workflows/simulations`.
- Put event-specific parameters and formal runs in `cases/<case-id>`.
- Keep large simulation snapshots and generated products outside Git; record their paths in the case
  `data_manifest.md`.

## Initial Targets

- Read existing MHD/binary/VTK simulation outputs.
- Convert simulation variables into common physical units.
- Make slices, projections, field-line diagnostics, and observation-comparison products.
- Reuse plotting and movie utilities from `src/solarutils/plotting` where possible.

## Shared Conventions

- Use SI or cgs units explicitly in variable names or metadata.
- Keep coordinate order explicit. Do not assume `(z, y, x)` or `(r, theta, phi)` silently.
- Keep geometry-specific interpolation and derivatives in the matching coordinate module.
- Put cross-geometry diagnostics, such as magnetic energy and divergence checks, behind a geometry
  argument or separate wrappers.

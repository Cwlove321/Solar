# workflows

Repeatable analysis workflows belong here. A workflow should combine reusable functions from
`solarutils` with project-specific configuration and command-line scripts.

## Cartesian magnetic-field projects

- `Cartesian_NLFFF/`: Cartesian potential-field and NLFFF extrapolation templates,
  including single-frame HMI/CEA preprocessing.
- `Cartesian_TMF/`: Data-driven Cartesian TMF templates and time-series
  preprocessing.
- `Cartesian_RBSL_Bipolar/`: Analytic bipolar-field model with an inserted
  RBSL flux rope, magnetofrictional relaxation, and zero-beta MHD evolution.
- `Spherical_RBSL/`: Portable spherical RBSL workflow covering boundary/path
  preprocessing, PFSS, magnetofriction, and MHD evolution.
- `Synthesize_Radiation/`: Synthetic-radiation workflows for Cartesian AMRVAC
  snapshots, including optically thin and optically thick AIA synthesis.
- `Helicity/`: Cartesian AMRVAC magnetic-energy and relative-helicity workflow,
  including uniform-grid preprocessing with the upstream `simesh` cloned into
  the local, Git-ignored `Packages/simesh` directory.

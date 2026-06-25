# Migration Notes

This file tracks how reusable code from the older working directories is promoted into
`src/solarutils`.

## First batch

Source directory: `/Users/xiaoxingyun/Desktop/pycode/solar_package/mylib`

| Old function or file | New location | Notes |
| --- | --- | --- |
| `mylib.mylib.curl` | `solarutils.fields.vector.curl` | Renamed `dim` to `axis_order`; validates shapes and dimensions. |
| `mylib.helicity.h_flux` | `solarutils.fields.helicity.h_flux` | Preserves original component order convention. |
| `mylib.helicity.e_flux` | `solarutils.fields.helicity.e_flux` | Preserves original component order convention. |
| `mylib.helicity.h_transpose` | `solarutils.fields.helicity.h_transpose` | Preserves original formula and adds clearer names. |
| `mylib.helicity.volume_integral` | `solarutils.fields.helicity.volume_integral` | Preserves boundary half-cell weighting. |
| `mylib.mylib.write_binary_data` | `solarutils.io.binary.write_magnetic_cube` | Uses a `MagneticCube` dataclass. |
| `mylib.mylib.read_binary_data` | `solarutils.io.binary.read_magnetic_cube` | Uses explicit `has_dx` option for legacy variants. |
| `mylib.mylib.init_boundary` | `solarutils.io.binary.init_boundary` | Preserves `(z, y, x)` convention. |
| `mylib.mylib.open_chase` | `solarutils.io.chase.read_chase` | Uses `NAXIS3` when available instead of hard-coded wavelength count. |
| `mylib.mylib.get_submap` | `solarutils.obs.maps.get_submap` | Keeps arcsec extent API. |
| `mylib.mylib.de_rotate` | `solarutils.obs.maps.derotate_to_reference` | Same core purpose, clearer name. |
| `mylib.mylib.aia_plot` | `solarutils.plotting.maps.plot_aia_map` | Removes hidden globals. |
| `mylib.mylib.hmi_plot` | `solarutils.plotting.maps.plot_hmi_los` | Removes hidden globals. |
| `mylib.mylib.vector_plot` | `solarutils.plotting.vectors.plot_signed_vectors` | Adds explicit `ax`, `bzmax`, and threshold arguments. |

## Next candidates

- `sharp/Calculate_sharpkeys_masked.py`: promote stable SHARP metric functions.
- `sharp/PIL_mask.py`: promote PIL mask utilities after adding tests on small arrays.
- `dem/`: package DEM inversion helpers with explicit relative imports.
- `dave/`: keep external DAVE4VM code separate; promote only thin wrappers and derived diagnostics.
- `tmf_prep/`: split reusable multigrid/boundary helpers from case-specific scripts.
- `rbsl_py/rbsl.py`: promote geometry and HMI disambiguation helpers after separating plotting.


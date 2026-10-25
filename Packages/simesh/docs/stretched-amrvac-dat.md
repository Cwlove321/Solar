# Stretched MPI-AMRVAC DAT coordinates

MPI-AMRVAC DAT snapshots do **not** store `stretch_dim`, `qstretch_baselevel`,
`nstretchedblocks_baselevel`, or per-block `rnode` coordinates. The header only
contains domain bounds and integer grid sizes. Consequently a stretched DAT
needs the original `amrvac.par` (or an equivalent configuration) in addition
to the snapshot.

Uniform DAT files keep their existing API:

```python
from simesh.amrvac import open_dataset
ds = open_dataset("data0000.dat")
```

For a stretched run, pass the par file:

```python
ds = open_dataset("data0000.dat", grid_config="amrvac.par")
```

or explicitly provide the coordinate settings:

```python
ds = open_dataset("data0000.dat", grid_config={
    "stretch_dim": ["uni", "none", "none"],
    "qstretch_baselevel": [1.041541581503301, None, None],
})
```

`stretch_type` is also accepted in place of `stretch_dim`; supported values are
`none`, `uni`, and `symm`. Symmetric stretching additionally requires
`nstretchedblocks_baselevel`. Header `xmin/xmax/domain_nx/block_nx` remain the
authoritative defaults; the supplied configuration only describes stretching.

Use `ds.uses_stretched_coordinates` to check the selected coordinate model.
The true local coordinates are `ds.block_node_coordinates` and
`ds.block_center_coordinates`, tuples ordered as `(r, theta, phi)` for
spherical data. Each entry has shape `(nleafs, block_nx[dimension] + 1)` for
nodes (or `block_nx[dimension]` for centers). `mesh.rnode[:, :2*ndim]` is also
corrected to the true block bounds; its historical single `dx` slot is a block
mean and must not be used for a stretched cell coordinate.

For `uni` stretching at AMR level 1, with `N=domain_nx` and `q` equal to
`qstretch_baselevel`, MPI-AMRVAC uses

`dx_first = (xmax-xmin) (1-q) / (1-q**N)` and
`x_i = xmin + dx_first (q**i-1)/(q-1)`.

At level `l>1`, `q_l=q_(l-1)**0.5` and
`dx_first,l=dx_first,l-1/(1+q_l)`. This makes each parent cell interval equal
to the sum of its two child intervals; simesh evaluates node indices at the
leaf level, rather than placing AMR blocks with a global uniform `dx`.

The zero-order `uniform_grid` sampler uses these exact nodes. For arbitrary
physical points, use `sample_points(points, interpolation="linear")`. It uses
the true leaf node and centre arrays, selects shared node faces with a
half-open ownership rule, and does not extrapolate a leaf's centre data into
its boundary half-cells. Same-level interfaces use a physical-coordinate
stencil from both sides. At a 2:1 coarse-fine interface the nearest physical
centre donors from all touching leaves are used and the affine reconstruction
is bounded by their value range; a rank-deficient stencil falls back to
nearest. This is intentionally variable-agnostic, so positive conserved
variables cannot become negative solely due to a block-face interpolation.

Ghost-cell exchange still assumes uniform cells. If a spherical/polar DAT is
opened without a configuration, simesh emits a strong warning because DAT
alone cannot tell a uniform run from a stretched one.

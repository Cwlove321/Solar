import numpy as np
import warnings
from .boundary import normalize_boundary_conditions
from .coordinates import AMRVACCoordinates, StretchConfig
from .datio import get_metadata, read_blocks_sequential
from .datio import update_header, write_datfile_from_sfc
from .dataset_base import DataSet
from .layouts import datau_to_udata
from simesh.utils.lib.amr.forest import AMRForest
from simesh.utils.lib.amr.mesh import AMRMesh

class AMRVACDataSet(DataSet):
    """
    AMRVAC specific implementation of DataSet.
    """
    def __init__(self, sfile: str, ghost_width: int = 0, boundary_conditions=None, grid_config=None):
        self.ghost_width = int(ghost_width)
        if self.ghost_width < 0:
            raise ValueError("ghost_width must be non-negative")
        self.boundary_conditions = boundary_conditions
        self.grid_config = grid_config
        super().__init__(sfile)

    def _build_mesh(
        self,
        nfields: int | None = None,
        field_indices: list[int] | None = None,
        boundary_table=None,
        normal_velocity_fields=None,
    ):
        if nfields is None:
            nfields = int(self.nw)
        if boundary_table is None or normal_velocity_fields is None:
            if field_indices is None:
                field_names = list(self.wnames[:nfields])
            else:
                field_names = [self.wnames[i] for i in field_indices]
            boundary_table, normal_velocity_fields = normalize_boundary_conditions(
                self.boundary_conditions,
                field_names,
                int(self.ndim),
            )

        return AMRMesh(
            self.ndim,
            self.block_nx,
            self.domain_nx,
            np.array(self.physical_domain[0], dtype=np.double),
            np.array(self.physical_domain[1], dtype=np.double),
            np.uint32(self.ghost_width),
            np.uint32(nfields),
            self.forest,
            boundary_table,
            normal_velocity_fields,
        )

    def load_metadata(self):
        """
        Load the metadata from the data file.
        """

        self.ng = self.ghost_width

        header, is_leaf, tree_info = get_metadata(self.sfile)
        self.metadata = header.copy()
        self.is_leaf = is_leaf.copy().astype(np.int32)
        self.tree_info = tree_info

        # Basic metadata
        self.ndim = np.uint32(header['ndim'])
        self.ndir = np.uint32(header['ndir'])
        self.nw = np.uint32(header['nw'])
        self.wnames = header['w_names']

        # AMR specific metadata
        self.nleafs = np.uint32(header['nleafs'])
        self.nparents = np.uint32(header['nparents'])
        self.levmax = np.uint32(header['levmax'])
        self.block_nx = header['block_nx'].astype(np.uint32)

        # Domain specific metadata
        self.domain_nx = header['domain_nx'].astype(np.uint32)
        self.physical_domain = np.array((header['xmin'], header['xmax']))
        self.periodic = header['periodic']
        self.geometry = header['geometry']
        self.coordinate_config = StretchConfig.from_user(self.grid_config, int(self.ndim))
        if self.coordinate_config.is_stretched and self.ghost_width:
            raise NotImplementedError(
                "Ghost-cell exchange currently assumes uniform cells. Open stretched DAT data with ghost_width=0; "
                "the exact physical-node zero-order sampler remains available."
            )
        self.coordinates = AMRVACCoordinates(
            header['xmin'], header['xmax'], header['domain_nx'], header['block_nx'],
            tree_info[0], tree_info[1], self.coordinate_config,
        )
        self.block_node_coordinates = self.coordinates.block_nodes()
        self.block_center_coordinates = self.coordinates.block_centers()
        self.uses_stretched_coordinates = self.coordinates.is_stretched
        self.metadata['coordinate_source'] = (
            'MPI-AMRVAC meshlist reconstruction' if self.uses_stretched_coordinates else 'DAT uniform metadata'
        )
        if self.uses_stretched_coordinates:
            # rnode still has its historical (min, max, one-dx) representation.
            # Its bounds are corrected below; consumers requiring cell coordinates must
            # use block_node_coordinates/block_center_coordinates.
            self.metadata['coordinate_warning'] = (
                'DAT has no stretch metadata; coordinates were reconstructed from grid_config.'
            )
        elif 'spherical' in self.geometry.lower() or 'polar' in self.geometry.lower():
            warnings.warn(
                'MPI-AMRVAC DAT does not record stretch_dim/qstretch_baselevel. '
                'This non-Cartesian DAT is being treated as uniform; pass grid_config=<amrvac.par or mapping> '
                'when the run used stretching.', RuntimeWarning, stacklevel=2,
            )

        # use nghostcells = 0
        ng3 = np.uint32(1)
        if int(self.ndim) == 3:
            ng3 = np.uint32(self.domain_nx[2]//self.block_nx[2])
        self.forest = AMRForest(self.ndim, 
                                np.uint32(self.domain_nx[0]//self.block_nx[0]), 
                                np.uint32(self.domain_nx[1]//self.block_nx[1]), 
                                ng3, 
                                self.is_leaf)
        self.mesh = self._build_mesh()
        if self.uses_stretched_coordinates:
            self._apply_true_block_bounds_to_mesh()

    def _apply_true_block_bounds_to_mesh(self):
        """Expose true block bounds through legacy ``mesh.rnode`` as well."""
        rnode = np.asarray(self.mesh.rnode)
        for dim, nodes in enumerate(self.block_node_coordinates):
            rnode[:, dim] = nodes[:, 0]
            rnode[:, int(self.ndim) + dim] = nodes[:, -1]
            # A single dx cannot represent a stretched block. Keep its mean solely
            # for legacy inspection; sampling is routed through true node arrays.
            rnode[:, 2 * int(self.ndim) + dim] = np.diff(nodes, axis=1).mean(axis=1)

    def print_metadata_impl(self):
        """
        furtherly print the metadata
        """
        print("====specific metadata=====")
        print("AMRVAC specific metadata:")
        print("Number of leaves: {self.nleafs}")
        print("Size of each block: {self.block_nx}")
        print("Maximum level: {self.levmax}")

    def load_data(self, field_indices: list[int] = None, boundary_conditions=None):
        """
        Load the amr 1d managed block data from the data file
        Can reload the data if the field_indices are different
        """
        if field_indices is None:
            field_indices = list(range(int(self.nw)))
        else:
            field_indices = [int(i) for i in field_indices]

        if boundary_conditions is not None:
            self.boundary_conditions = boundary_conditions

        field_names = [self.wnames[i] for i in field_indices]
        boundary_table, normal_velocity_fields = normalize_boundary_conditions(
            self.boundary_conditions,
            field_names,
            int(self.ndim),
        )
        self.boundary_condition_table = boundary_table
        self.normal_velocity_fields = normal_velocity_fields

        data = read_blocks_sequential(self.sfile, field_indices)
        if int(self.ndim) == 2:
            data = data[:, :, :, :, np.newaxis]
        if self.ghost_width > 0:
            self.mesh = self._build_mesh(
                len(field_indices),
                field_indices,
                boundary_table=boundary_table,
                normal_velocity_fields=normal_velocity_fields,
            )
            self.mesh.load_interior_data(np.ascontiguousarray(data, dtype=np.double))
            self.mesh.apply_ghost_cells()
            self.data = self.mesh.interior_view()
            if self.uses_stretched_coordinates:
                self._apply_true_block_bounds_to_mesh()
        else:
            self.data = data
        self.loaded_field_indices = field_indices

    def update_ghost_cells(self):
        if self.ghost_width <= 0:
            return
        if self.data is None:
            self.load_data(None)
            return
        self.mesh.apply_ghost_cells()
        self.data = self.mesh.interior_view()

    def exchange_ghost_cells(self):
        """
        Refresh ghost cells after mutating loaded interior block data.
        """
        self.update_ghost_cells()

    @property
    def has_ghost_cells(self):
        """
        Whether this dataset was opened with ghost-cell storage enabled.
        """
        return self.ghost_width > 0

    def blocks(self, include_ghosts: bool = False):
        """
        Return loaded block data in SFC layout.

        The default interior layout is ``(nleafs, nw, bx, by, bz)``. When
        ``include_ghosts`` is true and ``ghost_width > 0``, the returned view
        has shape ``(nleafs, nw, bx + 2g, by + 2g, bz + 2g)``.
        """
        if self.data is None:
            self.load_data(None)

        if not include_ghosts or self.ghost_width <= 0:
            return self.data

        padded = self.mesh.padded_view()
        return np.transpose(padded, (0, 4, 1, 2, 3))

    @property
    def loaded_field_indices(self):
        """
        Original file/header field indices currently loaded in self.data.
        """
        if not hasattr(self, "_loaded_field_indices"):
            self._loaded_field_indices = list(range(int(self.nw)))
        return self._loaded_field_indices

    @loaded_field_indices.setter
    def loaded_field_indices(self, value):
        self._loaded_field_indices = [int(i) for i in value]

    def _loaded_field_map(self):
        """
        Map original file/header field indices to column positions in self.data.
        """
        return {field_idx: i for i, field_idx in enumerate(self.loaded_field_indices)}
    
    def __getitem__(self, key):
        """
        Get data from the dataset on a uniform grid using NumPy-like slicing.

        The indexing uses a convention similar to ``numpy.mgrid`` where a
        complex-valued step encodes the total number of points in that
        direction.

        Examples
        --------
        - ``ds[::100j, ::100j, ::100j]``:
            build a uniform grid over the *entire* domain with a resolution of
            100 × 100 × 100 and return the full grid.

        - ``ds[50:100:200j, 50:100:200j, 50:100:200j]``:
            build a uniform grid over the *entire* domain with a resolution of
            200 × 200 × 200 and then return the sub-box ``50:100`` in each
            direction.

        Notes
        -----
        - Only 3D slicing with complex steps is currently supported.
        - The internal compute layout is ``(n_fields, nx, ny, nz)``.
        - This method converts that to the user-facing ``(nx, ny, nz, n_fields)``
          layout before applying the requested slices.
        """
        # Expect a 3D indexing key
        if not isinstance(key, tuple) or len(key) != 3:
            raise TypeError(
                "Indexing expects a tuple of three slices, e.g. "
                "[::100j, ::100j, ::100j] or [50:100:200j, 50:100:200j, 50:100:200j]"
            )

        slices = []
        nx = []

        for s in key:
            if not isinstance(s, slice) or not isinstance(s.step, complex):
                raise TypeError(
                    "Each index must be a slice with a complex step, e.g. ::100j or 50:100:200j"
                )

            # Total number of points along this axis
            n_axis = int(round(s.step.imag))
            if n_axis <= 0:
                raise ValueError("Complex step imag part must give a positive number of points")

            nx.append(n_axis)

            start = 0 if s.start is None else s.start
            stop = n_axis if s.stop is None else s.stop
            slices.append(slice(start, stop))

        # Build the full-domain uniform grid at the requested resolution.
        # This returns data with shape (nx, ny, nz, n_fields)
        full_grid = datau_to_udata(self.uniform_grid(nx))

        # Extract the requested sub-box
        return full_grid[slices[0], slices[1], slices[2], :]

    def uniform_grid(
        self,
        nx,
        xmin: list = None,
        xmax: list = None,
        field_indices: list[int] = None,
        interpolation: str = "zero",
    ):
        """
        Get the uniform grid data from the 1d AMR-managed data.

        Returns data in compute-oriented ``datau`` layout with shape
        ``(n_fields, nx, ny, nz)``. If you want user-facing ``udata`` layout
        ``(nx, ny, nz, n_fields)``, convert it with
        ``simesh.amrvac.layouts.datau_to_udata``.
        
        Parameters:
        -----------
        field_indices : list[int], optional
            Original file/header field indices to extract, corresponding to self.wnames.
            If None, uses all currently loaded fields.
        interpolation : {"zero", "linear"}, optional
            ``"zero"`` uses piecewise-constant sampling. ``"linear"`` uses
            trilinear interpolation from ghost-cell-padded mesh storage.
        """
        # Default to full domain if bounds are not specified
        if xmin is None:
            xmin = self.physical_domain[0]
        if xmax is None:
            xmax = self.physical_domain[1]

        nx = np.asarray(nx, dtype=np.uint32)
        if int(self.ndim) == 2:
            if nx.shape == (2,):
                nx = np.array([nx[0], nx[1], 1], dtype=np.uint32)
            elif nx.shape != (3,) or int(nx[2]) != 1:
                raise ValueError(f"2D resolution must have shape (nx, ny) or (nx, ny, 1), got {tuple(nx)}")
        elif nx.shape != (3,):
            raise ValueError(f"3D resolution must have three entries, got {tuple(nx)}")

        # Load all fields lazily if nothing is loaded yet.
        if self.data is None:
            self.load_data(None)

        loaded_field_map = self._loaded_field_map()

        if field_indices is not None:
            field_indices = [int(i) for i in field_indices]
            missing = [i for i in field_indices if i not in loaded_field_map]
            if missing:
                raise ValueError(
                    f"Requested field indices {missing} are not loaded. "
                    f"Currently loaded original field indices: {self.loaded_field_indices}"
                )

            loaded_columns = [loaded_field_map[i] for i in field_indices]
            data_to_use = self.data[:, loaded_columns, :, :, :]
            n_fields = len(loaded_columns)
        else:
            data_to_use = self.data
            n_fields = len(self.loaded_field_indices)
        
        uniform_grid = np.zeros((n_fields, int(nx[0]), int(nx[1]), int(nx[2])), dtype=np.double)

        if self.uses_stretched_coordinates:
            if interpolation.lower() not in ("zero", "nearest", "piecewise_constant"):
                raise NotImplementedError(
                    'Linear interpolation on stretched DAT grids is not yet available; '
                    'zero-order sampling uses the exact physical cell nodes.'
                )
            return self._uniform_grid_stretched(data_to_use, nx, np.asarray(xmin, dtype=float), np.asarray(xmax, dtype=float))

        interpolation = interpolation.lower()
        if interpolation in ("zero", "nearest", "piecewise_constant"):
            self.mesh.uniform_grid_zero_order(
                data_to_use,
                uniform_grid,
                nx,
                np.array(xmin, dtype=np.double),
                np.array(xmax, dtype=np.double),
            )
        elif interpolation in ("linear", "trilinear"):
            if self.ghost_width <= 0:
                raise ValueError("Linear interpolation requires opening the dataset with ghost_width > 0.")
            if field_indices is None:
                field_positions = np.arange(n_fields, dtype=np.uint32)
            else:
                field_positions = np.array(loaded_columns, dtype=np.uint32)
            self.mesh.uniform_grid_linear(
                uniform_grid,
                nx,
                np.array(xmin, dtype=np.double),
                np.array(xmax, dtype=np.double),
                field_positions,
            )
        else:
            raise ValueError(f"Unknown interpolation mode: {interpolation}")
        return uniform_grid

    def _uniform_grid_stretched(self, data, nx, xmin, xmax):
        """Sample leaf cells using their true physical node intervals."""
        nfields = data.shape[1]
        out = np.zeros((nfields, int(nx[0]), int(nx[1]), int(nx[2])), dtype=np.double)
        sample_centers = [xmin[d] + (np.arange(int(nx[d])) + 0.5) * (xmax[d] - xmin[d]) / int(nx[d]) for d in range(3)]
        nodes = self.block_node_coordinates
        active_dims = int(self.ndim)
        for leaf in range(data.shape[0]):
            output_indices = []
            local_indices = []
            valid = True
            for dim in range(active_dims):
                x = sample_centers[dim]
                leaf_nodes = nodes[dim][leaf]
                selected = np.flatnonzero((x >= leaf_nodes[0]) & ((x < leaf_nodes[-1]) | np.isclose(x, xmax[dim])))
                if not len(selected):
                    valid = False
                    break
                cell = np.searchsorted(leaf_nodes, x[selected], side='right') - 1
                cell = np.clip(cell, 0, len(leaf_nodes) - 2)
                output_indices.append(selected)
                local_indices.append(cell)
            if not valid:
                continue
            values = data[leaf]
            for dim, indices in enumerate(local_indices):
                values = np.take(values, indices, axis=dim + 1)
            target_slices = tuple(slice(int(indices[0]), int(indices[-1]) + 1) for indices in output_indices)
            if active_dims == 2:
                out[(slice(None),) + target_slices + (0,)] = values[..., 0]
            else:
                out[(slice(None),) + target_slices] = values
        return out

    def sample_points(self, points, field_indices: list[int] | None = None, interpolation: str = "linear"):
        """Sample loaded fields at physical points using true leaf coordinates.

        For stretched meshes, this is the coordinate-safe replacement for
        treating a block as uniformly spaced.  Values are cell centred, so
        linear interpolation is performed on ``block_center_coordinates``.

        A leaf node interval is used only to choose an owning leaf.  It is not
        a licence to extrapolate the owner's cell-centred data into its two
        half-cells at a block face.  Within the cell-centre box we use the
        usual tensor-product interpolation.  In a face half-cell we build a
        local, physical-coordinate affine stencil from all leaf blocks that
        touch the owner.  The result is limited to the donor range, making the
        coarse-fine fallback bounded even when its stencil is not a tensor
        product.  Points outside the mesh return ``nan``.
        """
        points = np.asarray(points, dtype=float)
        if points.ndim != 2 or points.shape[1] != int(self.ndim):
            raise ValueError(f"points must have shape (n, {int(self.ndim)})")
        interpolation = interpolation.lower()
        if interpolation not in {"linear", "trilinear", "nearest"}:
            raise ValueError("sample_points supports interpolation='linear' or 'nearest'")
        if self.data is None:
            self.load_data(None)

        loaded_map = self._loaded_field_map()
        if field_indices is None:
            columns = list(range(len(self.loaded_field_indices)))
        else:
            missing = [int(i) for i in field_indices if int(i) not in loaded_map]
            if missing:
                raise ValueError(f"Requested field indices {missing} are not loaded")
            columns = [loaded_map[int(i)] for i in field_indices]
        result = np.full((len(points), len(columns)), np.nan, dtype=float)
        nodes, centers = self.block_node_coordinates, self.block_center_coordinates
        ndim = int(self.ndim)

        # Use half-open node boxes.  Therefore a shared face belongs to its
        # positive-side leaf, independently of SFC/block traversal order.  The
        # physical outer upper boundary is the sole closed upper face.
        owner = np.full(len(points), -1, dtype=int)
        global_high = np.array([np.max(nodes[dim][:, -1]) for dim in range(ndim)])
        for leaf in range(int(self.nleafs)):
            inside = np.ones(len(points), dtype=bool)
            for dim in range(ndim):
                low, high = nodes[dim][leaf, 0], nodes[dim][leaf, -1]
                upper = (points[:, dim] < high) | ((high == global_high[dim]) & (points[:, dim] <= high))
                inside &= (points[:, dim] >= low) & upper
            owner[(owner < 0) & inside] = leaf

        for leaf in np.unique(owner[owner >= 0]):
            target = np.flatnonzero(owner == leaf)
            if not len(target):
                continue
            if interpolation == "nearest":
                cell_index = [np.abs(centers[dim][leaf, None, :] - points[target, dim, None]).argmin(axis=1) for dim in range(ndim)]
                values = self.data[leaf, columns]
                if ndim == 2:
                    result[target] = values[:, cell_index[0], cell_index[1], 0].T
                else:
                    result[target] = values[:, cell_index[0], cell_index[1], cell_index[2]].T
                continue

            # Never extrapolate an individual leaf's cell-centre array.  Keep
            # the fast, exact multilinear path for genuinely interior points.
            centre_inside = np.ones(len(target), dtype=bool)
            for dim in range(ndim):
                centre_inside &= ((points[target, dim] >= centers[dim][leaf, 0]) &
                                  (points[target, dim] <= centers[dim][leaf, -1]))
            interior = target[centre_inside]
            boundary = target[~centre_inside]

            if len(interior):
                self._sample_leaf_centres_linear(leaf, points, interior, columns, result)
            if len(boundary):
                self._sample_leaf_boundary_linear(leaf, points, boundary, columns, result)
        return result

    def _sample_leaf_centres_linear(self, leaf, points, target, columns, result):
        """Fast multilinear interpolation strictly inside one leaf centre box."""
        centers = self.block_center_coordinates
        ndim = int(self.ndim)
        lower_index = []
        weights = []
        for dim in range(ndim):
            axis = centers[dim][leaf]
            idx = np.searchsorted(axis, points[target, dim], side="right") - 1
            idx = np.clip(idx, 0, len(axis) - 2)
            lower_index.append(idx)
            weights.append((points[target, dim] - axis[idx]) / (axis[idx + 1] - axis[idx]))

        values = self.data[leaf, columns]
        if ndim == 2:
            i, j = lower_index
            wx, wy = weights
            sampled = (
                values[:, i, j, 0] * ((1 - wx) * (1 - wy))
                + values[:, i + 1, j, 0] * (wx * (1 - wy))
                + values[:, i, j + 1, 0] * ((1 - wx) * wy)
                + values[:, i + 1, j + 1, 0] * (wx * wy)
            )
        else:
            i, j, k = lower_index
            wx, wy, wz = weights
            sampled = np.zeros((len(columns), len(target)), dtype=float)
            for dx in (0, 1):
                for dy in (0, 1):
                    for dz in (0, 1):
                        weight = (wx if dx else 1 - wx) * (wy if dy else 1 - wy) * (wz if dz else 1 - wz)
                        sampled += values[:, i + dx, j + dy, k + dz] * weight
        result[target] = sampled.T

    def _touching_leaf_indices(self, leaf):
        """Return leaf blocks sharing a face, edge, or corner with ``leaf``."""
        cache = getattr(self, "_sampling_touching_leaf_cache", None)
        if cache is None:
            cache = self._sampling_touching_leaf_cache = {}
        if leaf in cache:
            return cache[leaf]
        nodes, ndim = self.block_node_coordinates, int(self.ndim)
        low = np.array([nodes[d][leaf, 0] for d in range(ndim)])
        high = np.array([nodes[d][leaf, -1] for d in range(ndim)])
        all_low = np.column_stack([nodes[d][:, 0] for d in range(ndim)])
        all_high = np.column_stack([nodes[d][:, -1] for d in range(ndim)])
        scale = max(1.0, float(np.max(np.abs(np.r_[low, high]))))
        tol = 64 * np.finfo(float).eps * scale
        # For two closed intervals their separation is zero when they overlap
        # or touch.  AMR leaves have disjoint interiors, so this selects only
        # the owner and geometrical neighbours, including 2:1 interfaces.
        separation = np.maximum(np.maximum(all_low - high, low - all_high), 0.0)
        cache[leaf] = np.flatnonzero(np.all(separation <= tol, axis=1))
        return cache[leaf]

    def _sample_leaf_boundary_linear(self, leaf, points, target, columns, result):
        """Bounded local affine reconstruction across a leaf boundary.

        The donor search is geometric, not SFC-order based.  At a same-level
        face its nearest centres contain the two opposing cell rows and an
        affine field is recovered exactly.  At coarse-fine faces the physically
        nearest centres (hence the locally finest available donor set) are
        used; rank-deficient stencils deliberately fall back to nearest.
        """
        ndim = int(self.ndim)
        centres = self.block_center_coordinates
        neighbour_leaves = self._touching_leaf_indices(int(leaf))
        donor_coordinates = []
        donor_values = []
        for neighbour in neighbour_leaves:
            axes = [centres[dim][neighbour] for dim in range(ndim)]
            mesh = np.meshgrid(*axes, indexing="ij")
            donor_coordinates.append(np.column_stack([axis.ravel() for axis in mesh]))
            values = self.data[neighbour, columns]
            if ndim == 2:
                donor_values.append(values[:, :, :, 0].reshape(len(columns), -1).T)
            else:
                donor_values.append(values.reshape(len(columns), -1).T)
        donor_coordinates = np.concatenate(donor_coordinates, axis=0)
        donor_values = np.concatenate(donor_values, axis=0)
        # Canonicalise donor order before resolving equal distances.  This is
        # observable on a shared face, where several cell centres can be
        # equidistant, and prevents SFC ordering from changing a result.
        canonical_order = np.lexsort(tuple(donor_coordinates[:, dim] for dim in range(ndim - 1, -1, -1)))
        donor_coordinates = donor_coordinates[canonical_order]
        donor_values = donor_values[canonical_order]

        # Normalising by the owner's physical cell widths makes donor choice
        # robust for stretched r/theta/phi coordinates.
        local_scale = np.array([
            np.median(np.diff(centres[dim][leaf])) for dim in range(ndim)
        ])
        local_scale = np.maximum(local_scale, np.finfo(float).eps)
        nstencil = min(len(donor_coordinates), max(12, 4 * (ndim + 1)))
        try:
            # SciPy is an optional dependency of simesh, hence the NumPy
            # fallback below.  Batched KD-tree queries are important for a
            # plotting plane that crosses thousands of AMR leaf faces.
            from scipy.spatial import cKDTree
            selected = cKDTree(donor_coordinates / local_scale).query(
                points[target] / local_scale, k=nstencil
            )[1]
            selected = np.asarray(selected).reshape(len(target), nstencil)
        except ImportError:
            selected = np.empty((len(target), nstencil), dtype=int)
            for row, point_index in enumerate(target):
                distance2 = np.sum(((donor_coordinates - points[point_index]) / local_scale) ** 2, axis=1)
                selected[row] = np.argsort(distance2, kind="stable")[:nstencil]

        for row, point_index in enumerate(target):
            # A tree is allowed to return equally distant points in any order;
            # restore canonical donor order before the least-squares solve.
            selected_row = np.sort(selected[row])
            x = donor_coordinates[selected_row]
            y = donor_values[selected_row]
            # Centre coordinates relative to the query point avoid an
            # ill-conditioned intercept for spherical radii near one.
            design = np.column_stack((np.ones(len(selected_row)), x - points[point_index]))
            if np.linalg.matrix_rank(design) < ndim + 1:
                distance2 = np.sum(((x - points[point_index]) / local_scale) ** 2, axis=1)
                result[point_index] = y[np.argmin(distance2)]
                continue
            coefficients, _, _, _ = np.linalg.lstsq(design, y, rcond=None)
            reconstructed = coefficients[0]
            # This is a slope limiter rather than a variable-specific fix: no
            # positive cell-centred scalar can become negative in a boundary
            # half-cell, and discontinuous AMR fields stay bounded.
            result[point_index] = np.clip(reconstructed, y.min(axis=0), y.max(axis=0))

    def uniform_full(self, field_indices: list[int] = None):
        """
        Return the full-domain uniform grid in compute-oriented ``datau``
        layout with shape ``(n_fields, nx, ny, nz)`` for datasets without
        refinement.
        """
        if int(self.levmax) != 1:
            raise ValueError("uniform_full() is only available when levmax == 1.")

        if self.data is None:
            self.load_data(None)

        loaded_field_map = self._loaded_field_map()

        if field_indices is not None:
            field_indices = [int(i) for i in field_indices]
            missing = [i for i in field_indices if i not in loaded_field_map]
            if missing:
                raise ValueError(
                    f"Requested field indices {missing} are not loaded. "
                    f"Currently loaded original field indices: {self.loaded_field_indices}"
                )

            loaded_columns = [loaded_field_map[i] for i in field_indices]
            data_to_use = self.data[:, loaded_columns, :, :, :]
            n_fields = len(loaded_columns)
        else:
            data_to_use = self.data
            n_fields = len(self.loaded_field_indices)

        if int(self.ndim) == 2:
            uniform_shape = (n_fields, int(self.domain_nx[0]), int(self.domain_nx[1]), 1)
        else:
            uniform_shape = (n_fields, int(self.domain_nx[0]), int(self.domain_nx[1]), int(self.domain_nx[2]))
        uniform_grid = np.zeros(uniform_shape, dtype=np.double)
        self.mesh.uniform_full_level1(data_to_use, uniform_grid)

        expected_shape = uniform_shape
        assert uniform_grid.shape == expected_shape, \
            f"uniform_full result shape mismatch: {uniform_grid.shape} != {expected_shape}"

        return uniform_grid

    def write_datfile(self, sfile: str, overwrite: bool = False):

        if self.data is None:
            self.load_data(None)

        updated_header = update_header(
            self.metadata,
            nw=len(self.loaded_field_indices),
            w_names=[self.wnames[i] for i in self.loaded_field_indices],
        )
        data0 = np.asarray(self.data)
        return write_datfile_from_sfc(sfile, data0, updated_header, self.is_leaf, self.tree_info, overwrite=overwrite)

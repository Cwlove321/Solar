import numpy as np

from simesh.amrvac.amrvac_dataset import AMRVACDataSet


def _dataset_from_boxes(boxes, cells_per_box, field):
    """Small in-memory AMR dataset with cell-centred values for sampler tests."""
    ndim = 3
    ds = object.__new__(AMRVACDataSet)
    ds.ndim = ndim
    ds.nleafs = len(boxes)
    ds.loaded_field_indices = [0]
    ds.data = np.empty((len(boxes), 1, cells_per_box, cells_per_box, cells_per_box))
    nodes = [np.empty((len(boxes), cells_per_box + 1)) for _ in range(ndim)]
    centres = [np.empty((len(boxes), cells_per_box)) for _ in range(ndim)]
    for leaf, (low, high) in enumerate(boxes):
        axes = []
        for dim in range(ndim):
            nodes[dim][leaf] = np.linspace(low[dim], high[dim], cells_per_box + 1)
            centres[dim][leaf] = 0.5 * (nodes[dim][leaf, :-1] + nodes[dim][leaf, 1:])
            axes.append(centres[dim][leaf])
        x, y, z = np.meshgrid(*axes, indexing="ij")
        ds.data[leaf, 0] = field(x, y, z)
    ds.block_node_coordinates = tuple(nodes)
    ds.block_center_coordinates = tuple(centres)
    return ds


def _two_level_one_boxes():
    return [
        (np.array([0.0, 0.0, 0.0]), np.array([0.5, 1.0, 1.0])),
        (np.array([0.5, 0.0, 0.0]), np.array([1.0, 1.0, 1.0])),
    ]


def test_same_level_face_linear_field_is_continuous_and_exact():
    field = lambda x, y, z: 2.0 * x - 3.0 * y + 0.5 * z + 7.0
    ds = _dataset_from_boxes(_two_level_one_boxes(), 4, field)
    points = np.array([[x, 0.41, 0.62] for x in [0.43, 0.49, 0.5, 0.51, 0.57]])
    assert np.allclose(ds.sample_points(points)[:, 0], field(*points.T), atol=1e-12)


def test_shared_face_result_does_not_depend_on_leaf_order():
    field = lambda x, y, z: x + 2 * y - z
    original = _dataset_from_boxes(_two_level_one_boxes(), 4, field)
    reversed_ds = _dataset_from_boxes(list(reversed(_two_level_one_boxes())), 4, field)
    points = np.array([[0.5, 0.2, 0.7], [0.49, 0.2, 0.7], [0.51, 0.2, 0.7]])
    assert np.array_equal(original.sample_points(points), reversed_ds.sample_points(points))


def test_positive_steep_field_stays_finite_and_positive_at_face():
    field = lambda x, y, z: np.exp(12.0 * x) * (1.0 + 0.1 * y)
    ds = _dataset_from_boxes(_two_level_one_boxes(), 4, field)
    points = np.array([[x, 0.45, 0.55] for x in np.linspace(0.44, 0.56, 31)])
    sampled = ds.sample_points(points)[:, 0]
    assert np.all(np.isfinite(sampled))
    assert np.all(sampled > 0.0)


def test_two_to_one_coarse_fine_interface_is_bounded_and_exact_for_linear_field():
    # The left leaf has dx=0.125.  The two right leaves have dx=0.0625,
    # representing the fine side of a 2:1 interface.
    boxes = [
        (np.array([0.0, 0.0, 0.0]), np.array([0.5, 1.0, 1.0])),
        (np.array([0.5, 0.0, 0.0]), np.array([0.75, 1.0, 1.0])),
        (np.array([0.75, 0.0, 0.0]), np.array([1.0, 1.0, 1.0])),
    ]
    field = lambda x, y, z: 3.0 * x + y - 2.0 * z
    ds = _dataset_from_boxes(boxes, 4, field)
    points = np.array([[x, 0.37, 0.58] for x in np.linspace(0.43, 0.60, 20)])
    sampled = ds.sample_points(points)[:, 0]
    assert np.all(np.isfinite(sampled))
    assert np.allclose(sampled, field(*points.T), atol=1e-12)


def test_single_block_internal_linear_matches_previous_multilinear_result():
    field = lambda x, y, z: 4.0 * x - 2.0 * y + z
    ds = _dataset_from_boxes([(np.zeros(3), np.ones(3))], 4, field)
    points = np.array([[0.31, 0.46, 0.62], [0.69, 0.53, 0.38]])
    assert np.allclose(ds.sample_points(points)[:, 0], field(*points.T), atol=1e-12)

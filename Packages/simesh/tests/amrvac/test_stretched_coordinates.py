import numpy as np
import pytest

from simesh.amrvac.coordinates import AMRVACCoordinates, StretchConfig, grid_config_from_par
from simesh.amrvac import open_dataset, write_datfile_from_uniform


def test_unidirectional_nodes_follow_amrvac_and_refine_parent_intervals():
    config = StretchConfig.from_user(
        {"stretch_dim": ["uni"], "qstretch_baselevel": [2.0]}, ndim=1
    )
    # Two level-one blocks; the two level-two blocks are the children of the
    # first level-one block. DAT tree ids are one-based.
    coordinates = AMRVACCoordinates(
        xmin=[1.0], xmax=[16.0], domain_nx=[4], block_nx=[2],
        block_levels=[1, 2, 2], block_ids=[[1], [1], [2]], config=config,
    )
    nodes = coordinates.block_nodes()[0]
    centers = coordinates.block_centers()[0]

    assert np.allclose(nodes[0], [1.0, 2.0, 4.0])
    assert np.allclose(nodes[1], [1.0, 1.0 + 1.0 / (1.0 + np.sqrt(2.0)), 2.0])
    assert np.allclose(nodes[2], [2.0, 2.0 + 2.0 / (1.0 + np.sqrt(2.0)), 4.0])
    assert np.allclose(np.diff(nodes[0]), [np.diff(nodes[1]).sum(), np.diff(nodes[2]).sum()])
    assert np.allclose(centers[0], 0.5 * (nodes[0, :-1] + nodes[0, 1:]))


def test_symmetric_nodes_are_mirrored_and_cover_domain():
    config = StretchConfig.from_user(
        {"stretch_type": ["symm"], "qstretch_baselevel": [2.0], "nstretchedblocks_baselevel": [2]}, ndim=1
    )
    coordinates = AMRVACCoordinates(
        xmin=[0.0], xmax=[10.0], domain_nx=[8], block_nx=[2],
        block_levels=[1], block_ids=[[1]], config=config,
    )
    all_nodes = coordinates._level_nodes(1, 0)
    assert all_nodes[0] == 0.0
    assert all_nodes[-1] == 10.0
    assert np.allclose(np.diff(all_nodes), np.diff(all_nodes)[::-1])


def test_stretched_config_requires_q_and_par_parser(tmp_path):
    with pytest.raises(ValueError, match="qstretch_baselevel"):
        StretchConfig.from_user({"stretch_dim": ["uni"]}, ndim=1)
    par = tmp_path / "amrvac.par"
    par.write_text("&meshlist\n stretch_dim(1)='uni'\n qstretch_baselevel(1)=1.041d0\n/\n")
    parsed = grid_config_from_par(par, ndim=1)
    assert parsed["stretch_type"] == ["uni"]
    assert parsed["qstretch_baselevel"] == [1.041]


def test_dataset_exposes_true_nodes_and_corrected_rnode_bounds(tmp_path):
    dat = tmp_path / "uniform-topology.dat"
    write_datfile_from_uniform(
        dat, np.zeros((4, 4, 4, 1)), ["w"],
        xmin=np.array([1.0, 0.0, 0.0]), xmax=np.array([16.0, 1.0, 1.0]),
        block_nx=np.array([2, 2, 2]), overwrite=True,
    )
    ds = open_dataset(dat, grid_config={
        "stretch_dim": ["uni", "none", "none"],
        "qstretch_baselevel": [2.0, None, None],
    })
    assert ds.uses_stretched_coordinates
    radial_nodes = ds.block_node_coordinates[0]
    assert np.allclose(radial_nodes[0], [1.0, 2.0, 4.0])
    assert np.allclose(np.asarray(ds.mesh.rnode)[0, [0, 3]], [1.0, 4.0])


def test_sample_points_uses_true_cell_centers(tmp_path):
    dat = tmp_path / "linear.dat"
    centers = (np.arange(4) + 0.5) / 4.0
    x, y, z = np.meshgrid(centers, centers, centers, indexing="ij")
    data = (2.0 * x + 3.0 * y - z)[..., None]
    write_datfile_from_uniform(
        dat, data, ["w"], xmin=np.zeros(3), xmax=np.ones(3),
        block_nx=np.array([2, 2, 2]), overwrite=True,
    )
    ds = open_dataset(dat)
    ds.load_data()
    points = np.array([[0.3, 0.3, 0.3], [0.7, 0.7, 0.7]])
    assert np.allclose(ds.sample_points(points)[:, 0], 2 * points[:, 0] + 3 * points[:, 1] - points[:, 2])

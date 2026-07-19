import numpy as np

from workflows.Cartesian_NLFFF.Preprocess.prepare_nlfff_boundary import (
    parse_args as parse_nlfff_args,
)
from workflows.Cartesian_NLFFF.Preprocess.nlfff_core.cartesian import (
    multigrid_reduce_components,
)
from workflows.Cartesian_TMF.Preprocess.driven_boundary import reduce_components
from workflows.Cartesian_TMF.Preprocess.prepare_tmf_boundaries import (
    parse_args as parse_tmf_args,
)


def test_tmf_reducer_matches_shared_cartesian_block_average():
    bx = np.arange(24.0).reshape(4, 6)
    by = bx + 100.0
    bz = bx - 50.0

    expected = multigrid_reduce_components(bx, by, bz, level=2)
    actual = reduce_components(bx, by, bz, level=2, preprocess=False)

    for actual_component, expected_component in zip(actual, expected, strict=True):
        np.testing.assert_allclose(actual_component, expected_component)


def test_nlfff_and_tmf_preprocess_defaults():
    nlfff = parse_nlfff_args(
        [
            "--output-dir", "out",
            "--x0", "0", "--y0", "0", "--nx", "4", "--ny", "4",
        ]
    )
    assert nlfff.preprocess is True
    assert nlfff.preprocess_mu3 == 0.1
    assert nlfff.preprocess_mu4 == 0.1
    assert nlfff.preprocess_max_iter == 5000
    assert nlfff.preprocess_tol == 1.0e-4

    tmf = parse_tmf_args(
        [
            "--input-dir", "in", "--output-dir", "out",
            "--x1", "0", "--x2", "4", "--y1", "0", "--y2", "4",
        ]
    )
    assert tmf.prep == 0
    assert tmf.nlfff_prep == 0

"""Check the grid policy before any potential/NLFFF boundary output is written."""

import importlib
import json

import numpy as np
import pytest
from astropy.io import fits


@pytest.fixture(params=["Cartesian_NLFFF", "Cartesian_TMF"])
def workflow(request):
    return importlib.import_module(f"workflows.{request.param}.Preprocess.prepare_nlfff_boundary")


def test_grid_factorization_and_invalid_dimensions(workflow):
    report = workflow.validate_physical_grid(360, 216)
    assert report["factorization"] == {"nx": {"2": 3, "3": 2, "5": 1}, "ny": {"2": 3, "3": 3}}
    with pytest.raises(ValueError, match=r"nx=370.*remaining factor 37.*ny=220.*remaining factor 11"):
        workflow.validate_physical_grid(370, 220)
    for bad in (0, -2, 2.5, True):
        with pytest.raises(ValueError, match="positive integer"):
            workflow.validate_physical_grid(bad, 216)


def test_cli_rejects_bad_binned_grid_before_loading_or_overwriting(workflow, tmp_path):
    out = tmp_path / "existing"
    out.mkdir()
    sentinel = out / "potential_boundary.dat"
    sentinel.write_bytes(b"preserve existing input")
    with pytest.raises(ValueError, match="remaining factor 37"):
        workflow.main(["--output-dir", str(out), "--x0", "10", "--y0", "10",
                       "--nx", "1480", "--ny", "880", "--level", "3"])
    assert sentinel.read_bytes() == b"preserve existing input"
    assert list(out.iterdir()) == [sentinel]


def test_block_validation_uses_physical_grid(workflow):
    base = ["--output-dir", "unused", "--x0", "30", "--y0", "18",
            "--nx", "1440", "--ny", "864", "--level", "3"]
    args = workflow.parse_args(base + ["--block-nx1", "12", "--block-nx2", "12", "--block-nx3", "12"])
    assert workflow.validate_physical_window(args, 4)["domain_shape_xyz"] == [360, 216, 216]
    args.block_nx2 = 10
    with pytest.raises(ValueError, match="divide domain_nx2=216"):
        workflow.validate_physical_window(args, 4)


def test_shared_writer_checks_grid_but_accepts_non235_file_ghost_shape(workflow, tmp_path):
    bx = np.ones((8, 10))
    result = workflow.write_nlfff_from_cartesian(tmp_path / "good", bx, bx, bx, level=1, dx_cm=1e8)
    assert result["grid_validation"]["physical_shape_xy"] == [10, 8]
    header = np.fromfile(tmp_path / "good/nlfff_boundary.dat", dtype=np.int32, count=2)
    assert header.tolist() == [14, 12]  # 14 contains 7, but these are file ghosts.
    bad = np.ones((8, 14))
    with pytest.raises(ValueError, match="remaining factor 7"):
        workflow.write_nlfff_from_cartesian(tmp_path / "bad", bad, bad, bad, level=1, dx_cm=1e8)
    assert not (tmp_path / "bad").exists()


def test_cli_writes_validation_metadata_and_consistent_binary_headers(workflow, tmp_path):
    cea = tmp_path / "cea"
    cea.mkdir()
    header = fits.Header({"CDELT1": 0.03, "CDELT2": 0.03})
    for name in ("Br", "Bt", "Bp"):
        fits.writeto(cea / f"{name}.fits", np.ones((16, 20)), header)
    out = tmp_path / "output"
    assert workflow.main(["--cea-dir", str(cea), "--output-dir", str(out),
                          "--x0", "4", "--y0", "4", "--nx", "12", "--ny", "8", "--level", "2",
                          "--block-nx1", "2", "--block-nx2", "2", "--block-nx3", "2",
                          "--no-preprocess", "--no-quicklook"]) == 0
    metadata = json.loads((out / "parameters.json").read_text())
    assert metadata["grid_validation"]["passed"] is True
    assert metadata["grid_validation"]["physical_shape_xy"] == [6, 4]
    assert np.fromfile(out / "potential_boundary.dat", np.int32, count=2).tolist() == [6, 4]
    assert np.fromfile(out / "nlfff_boundary.dat", np.int32, count=2).tolist() == [10, 8]

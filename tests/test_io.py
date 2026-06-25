import numpy as np
from astropy.io import fits

from solarutils.io.binary import MagneticCube, init_boundary, read_magnetic_cube, write_magnetic_cube
from solarutils.io.chase import normalize_chase_header


def test_magnetic_cube_roundtrip(tmp_path):
    cube = MagneticCube(
        b1=np.arange(24, dtype=float).reshape(2, 3, 4),
        b2=np.arange(24, 48, dtype=float).reshape(2, 3, 4),
        b3=np.arange(48, 72, dtype=float).reshape(2, 3, 4),
        dx=0.5,
    )
    path = tmp_path / "cube.binary"

    write_magnetic_cube(path, cube)
    loaded = read_magnetic_cube(path)

    assert loaded.dx == 0.5
    assert np.array_equal(loaded.b1, cube.b1)
    assert np.array_equal(loaded.b2, cube.b2)
    assert np.array_equal(loaded.b3, cube.b3)


def test_init_boundary_returns_signed_flux(tmp_path):
    b1 = np.ones((2, 3, 4))
    b2 = np.ones((2, 3, 4)) * 2
    b3 = np.ones((2, 3, 4)) * 3

    total = init_boundary(tmp_path / "bd.binary", b1, b2, b3)

    assert total == 0


def test_normalize_chase_header_adds_sunpy_map_keywords():
    header = fits.Header()
    header["DATE_OBS"] = "2024-10-08T23:00:00"
    header["RSUN_REF"] = 696000000.0
    header["RSUN_OBS"] = 960.0

    normalized = normalize_chase_header(header)

    assert normalized["CTYPE1"] == "HPLN-TAN"
    assert normalized["CTYPE2"] == "HPLT-TAN"
    assert normalized["DATE-OBS"] == "2024-10-08T23:00:00"
    assert normalized["DSUN_OBS"] > 0
    assert normalized["HGLN_OBS"] == 0.0
    assert normalized["HGLT_OBS"] == 0.0

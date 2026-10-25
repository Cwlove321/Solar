from types import SimpleNamespace

import numpy as np

from solarutils.obs import maps


def test_aia_prep_runs_calibration_steps_in_order(monkeypatch):
    calls = []
    input_map = SimpleNamespace(wavelength="171 Angstrom")
    pointing_table = object()
    psf = object()

    def update_pointing(smap, *, pointing_table=None):
        calls.append(("update_pointing", smap, pointing_table))
        return SimpleNamespace(wavelength=smap.wavelength, step="pointing")

    def respike(smap, *, spikes=None):
        calls.append(("respike", smap.step, spikes))
        return SimpleNamespace(wavelength=smap.wavelength, step="respike")

    def deconvolve(smap, *, psf=None, iterations=25):
        calls.append(("deconvolve", smap.step, psf, iterations))
        return SimpleNamespace(wavelength=smap.wavelength, step="deconvolve")

    def register(smap):
        calls.append(("register", smap.step))
        return SimpleNamespace(wavelength=smap.wavelength, step="register")

    fake_calibrate = SimpleNamespace(update_pointing=update_pointing, respike=respike, register=register)
    fake_psf = SimpleNamespace(deconvolve=deconvolve)
    monkeypatch.setattr(maps, "_aiapy_imports", lambda: (fake_calibrate, fake_psf))

    result = maps.aia_prep(
        input_map,
        pointing_table=pointing_table,
        psf=psf,
        deconvolve_psf=True,
        respike_spikes=("positions", "values"),
        deconvolve_iterations=5,
    )

    assert result.step == "register"
    assert [call[0] for call in calls] == ["update_pointing", "respike", "deconvolve", "register"]
    assert calls[0][2] is pointing_table
    assert calls[2][2] is psf
    assert calls[2][3] == 5


def test_aia_prep_can_skip_expensive_steps(monkeypatch):
    calls = []
    input_map = SimpleNamespace(wavelength="171 Angstrom")

    def update_pointing(smap, *, pointing_table=None):
        calls.append("update_pointing")
        return smap

    fake_calibrate = SimpleNamespace(update_pointing=update_pointing)
    fake_psf = SimpleNamespace()
    monkeypatch.setattr(maps, "_aiapy_imports", lambda: (fake_calibrate, fake_psf))

    result = maps.aia_prep(input_map, deconvolve_psf=False, register_map=False)

    assert result is input_map
    assert calls == ["update_pointing"]


def test_aia_psf_path_uses_zero_padded_wavelength(tmp_path):
    assert maps.aia_psf_path(94, tmp_path) == tmp_path / "aia_0094_psf.npz"


def test_load_aia_psf_reads_cached_array(tmp_path):
    psf = np.ones((2, 2), dtype=np.float32)
    path = maps.aia_psf_path(171, tmp_path)
    np.savez_compressed(path, psf=psf, wavelength_angstrom=171)

    loaded = maps.load_aia_psf(171, tmp_path)

    assert loaded.dtype == np.float32
    assert np.array_equal(loaded, psf)


def test_crop_psf_to_shape_returns_normalized_center_crop():
    psf = np.arange(64, dtype=np.float32).reshape(8, 8)

    cropped = maps.crop_psf_to_shape(psf, (4, 4))

    assert cropped.shape == (4, 4)
    assert np.isclose(cropped.sum(), 1.0)
    assert np.allclose(cropped, psf[2:6, 2:6] / psf[2:6, 2:6].sum())


def test_crop_psf_to_shape_rejects_odd_dimensions():
    psf = np.ones((8, 8), dtype=np.float32)

    try:
        maps.crop_psf_to_shape(psf, (4, 5))
    except ValueError as exc:
        assert "even dimensions" in str(exc)
    else:
        raise AssertionError("odd target dimensions should fail")


def test_pad_extent_clips_to_map_footprint():
    smap = SimpleNamespace(
        bottom_left_coord=SimpleNamespace(
            Tx=SimpleNamespace(to_value=lambda unit: -100.0),
            Ty=SimpleNamespace(to_value=lambda unit: -200.0),
        ),
        top_right_coord=SimpleNamespace(
            Tx=SimpleNamespace(to_value=lambda unit: 100.0),
            Ty=SimpleNamespace(to_value=lambda unit: 200.0),
        ),
    )

    padded = maps.pad_extent((-50.0, 50.0, -50.0, 50.0), 100.0, smap)

    assert padded == (-100.0, 100.0, -150.0, 150.0)


def test_sharpen_map_preserves_metadata_and_matches_unsharp_mask(monkeypatch):
    data = np.arange(9, dtype=np.float32).reshape(3, 3)
    meta = {"instrument": "AIA"}
    input_map = SimpleNamespace(data=data, meta=meta)

    def fake_filter(values, sigma=1.0):
        assert sigma == 2.0
        return np.ones_like(values)

    class FakeMap:
        def __init__(self, data, meta):
            self.data = data
            self.meta = meta

    monkeypatch.setattr(maps.ndimage, "gaussian_filter", fake_filter)
    monkeypatch.setattr(maps, "Map", FakeMap)

    sharpened = maps.sharpen_map(input_map, sigma=2.0, weight=1.5)

    expected = data.astype(np.float64) + (data.astype(np.float64) - 1.0) * 1.5
    assert np.array_equal(sharpened.data, expected)
    assert sharpened.meta == meta
    assert sharpened.meta is not meta


def test_multiscale_sharpen_map_combines_unsharp_masks(monkeypatch):
    data = np.arange(9, dtype=np.float32).reshape(3, 3)
    meta = {"instrument": "AIA"}
    input_map = SimpleNamespace(data=data, meta=meta)
    blurred_by_sigma = {
        1.0: np.ones_like(data) * 1.0,
        2.0: np.ones_like(data) * 2.0,
    }
    calls = []

    def fake_filter(values, sigma=1.0):
        calls.append(sigma)
        return blurred_by_sigma[sigma]

    class FakeMap:
        def __init__(self, data, meta):
            self.data = data
            self.meta = meta

    monkeypatch.setattr(maps.ndimage, "gaussian_filter", fake_filter)
    monkeypatch.setattr(maps, "Map", FakeMap)

    sharpened = maps.multiscale_sharpen_map(input_map, sigmas=(1.0, 2.0), weights=(1.5, 0.25))

    expected = data.astype(np.float64)
    expected += (data.astype(np.float64) - 1.0) * 1.5
    expected += (data.astype(np.float64) - 2.0) * 0.25
    assert calls == [1.0, 2.0]
    assert np.array_equal(sharpened.data, expected)
    assert sharpened.meta == meta
    assert sharpened.meta is not meta


def test_multiscale_sharpen_map_accepts_scalar_weight(monkeypatch):
    data = np.arange(4, dtype=np.float32).reshape(2, 2)
    input_map = SimpleNamespace(data=data, meta={})

    class FakeMap:
        def __init__(self, data, meta):
            self.data = data
            self.meta = meta

    monkeypatch.setattr(maps.ndimage, "gaussian_filter", lambda values, sigma=1.0: np.zeros_like(values))
    monkeypatch.setattr(maps, "Map", FakeMap)

    sharpened = maps.multiscale_sharpen_map(input_map, sigmas=(1.0, 2.0), weights=0.5)

    assert np.array_equal(sharpened.data, data.astype(np.float64) * 2.0)


def test_multiscale_sharpen_map_rejects_mismatched_scales_and_weights():
    input_map = SimpleNamespace(data=np.ones((2, 2)), meta={})

    try:
        maps.multiscale_sharpen_map(input_map, sigmas=(1.0, 2.0), weights=(1.0,))
    except ValueError as exc:
        assert "same length" in str(exc)
    else:
        raise AssertionError("mismatched sigmas and weights should fail")


def test_derotated_reference_meta_uses_reference_geometry_and_source_time():
    reference_map = SimpleNamespace(
        meta={
            "crpix1": 10.0,
            "crpix2": 20.0,
            "date-obs": "2024-10-08T23:00:00.000",
            "t_obs": "2024-10-08T23:00:01.000",
            "exptime": 2.0,
            "bunit": "DN",
            "blank": -2147483648,
        }
    )
    source_map = SimpleNamespace(
        meta={
            "date-obs": "2024-10-08T23:10:00.000",
            "t_obs": "2024-10-08T23:10:01.000",
            "exptime": 2.9,
            "bunit": "DN",
        },
        date=SimpleNamespace(isot="2024-10-08T23:10:00.000"),
    )

    meta = maps._derotated_reference_meta(
        reference_map,
        source_map,
        exposure_time=2.0,
        exposure_normalized=True,
    )

    assert meta["crpix1"] == 10.0
    assert meta["crpix2"] == 20.0
    assert meta["date-obs"] == "2024-10-08T23:10:00.000"
    assert meta["t_obs"] == "2024-10-08T23:10:01.000"
    assert meta["exptime"] == 2.0
    assert meta["bunit"] == "DN / s"
    assert "blank" not in meta
    assert reference_map.meta["date-obs"] == "2024-10-08T23:00:00.000"


def test_map_exposure_value_returns_quantity_value():
    exposure = SimpleNamespace(value=2.5)
    smap = SimpleNamespace(exposure_time=exposure)

    assert maps._map_exposure_value(smap) == 2.5


def test_map_exposure_value_returns_default_for_missing_exposure():
    smap = SimpleNamespace(exposure_time=None)

    assert maps._map_exposure_value(smap, default=None) is None


def test_derotated_reference_meta_can_omit_exposure_keywords():
    reference_map = SimpleNamespace(meta={"date-obs": "reference", "exptime": 2.0, "bunit": "Gauss"})
    source_map = SimpleNamespace(meta={"date-obs": "source", "bunit": "Gauss"}, date=SimpleNamespace(isot="source"))

    meta = maps._derotated_reference_meta(reference_map, source_map, exposure_time=None, exposure_normalized=False)

    assert meta["date-obs"] == "source"
    assert "exptime" not in meta
    assert meta["bunit"] == "Gauss"


def test_derotate_to_reference_rejects_unknown_method():
    try:
        maps.derotate_to_reference("unused.fits", object(), method="bad")
    except ValueError as exc:
        assert "method must be" in str(exc)
    else:
        raise AssertionError("unknown derotation method should fail")


def test_derotate_to_reference_pixel_method_uses_pixel_helper(monkeypatch):
    source_map = SimpleNamespace()
    reference_map = SimpleNamespace()
    calls = []

    class FakeMap:
        def __init__(self, path):
            self.path = path

    def fake_pixel(smap, reference, *, order=1, fill_value=None):
        calls.append((smap, reference, order, fill_value))
        return source_map

    monkeypatch.setattr(maps, "Map", FakeMap)
    monkeypatch.setattr(maps, "derotate_to_reference_pixel", fake_pixel)

    result = maps.derotate_to_reference("aia_test.fits", reference_map, order=3, fill_value=-1.0)

    assert result is source_map
    assert len(calls) == 1
    assert calls[0][1:] == (reference_map, 3, -1.0)
    assert isinstance(calls[0][0], FakeMap)


def test_align_map_sequence_by_content_removes_residual_pointing_jump(monkeypatch):
    class FakeMap:
        def __init__(self, data, meta):
            self.data = data
            self.meta = meta

    rng = np.random.default_rng(17)
    reference_data = rng.normal(size=(128, 128))
    jumped_data = np.roll(reference_data, shift=(4, 1), axis=(0, 1))
    monkeypatch.setattr(maps, "Map", FakeMap)

    corrected, shifts = maps.align_map_sequence_by_content(
        [FakeMap(jumped_data, {"id": 0}), FakeMap(reference_data, {"id": 1})],
        reference_index=1,
        upsample_factor=50,
        region=(slice(8, -8), slice(8, -8)),
    )

    residual, _, _ = __import__("skimage.registration", fromlist=["phase_cross_correlation"]).phase_cross_correlation(
        corrected[0].data[8:-8, 8:-8], corrected[1].data[8:-8, 8:-8], upsample_factor=50
    )
    assert np.allclose(shifts[0], (-4.0, -1.0), atol=0.1)
    assert np.allclose(shifts[1], (0.0, 0.0))
    assert np.linalg.norm(residual) < 0.1
    assert corrected[1].meta == {"id": 1}


def test_align_map_sequence_by_content_rejects_mismatched_shapes():
    maps_in = [
        SimpleNamespace(data=np.ones((2, 2))),
        SimpleNamespace(data=np.ones((3, 2))),
    ]

    try:
        maps.align_map_sequence_by_content(maps_in)
    except ValueError as exc:
        assert "identical data shapes" in str(exc)
    else:
        raise AssertionError("mismatched maps should be rejected")


def test_align_map_sequence_by_content_rejects_empty_region():
    maps_in = [SimpleNamespace(data=np.ones((2, 2))), SimpleNamespace(data=np.ones((2, 2)))]

    try:
        maps.align_map_sequence_by_content(maps_in, region=(slice(0, 0), slice(None)))
    except ValueError as exc:
        assert "at least one pixel" in str(exc)
    else:
        raise AssertionError("empty registration region should be rejected")

import numpy as np

from solarutils.sim import CartesianGrid, SphericalGrid
from solarutils.sim.cartesian import (
    boundary_metadata,
    crop_components,
    crop_components_with_padding,
    format_guo_amrvac_parameters,
    local_cartesian_from_heliographic,
    multigrid_reduce_components,
    preprocess_cartesian_field,
    write_nlfff_boundary,
    write_potential_boundary,
    _periodic_laplace,
)
from solarutils.sim.cea import (
    cea_lon_lat_from_xy,
    cea_xy_centers,
    gaussian_block_reduce,
    heliographic_to_local_cartesian_boundary,
    hmi_native_vector_components,
    make_cea_patch_grid,
    native_to_heliographic_components,
)
from solarutils.sim.spherical import (
    angular_distance,
    flux_rope_axis_from_points,
    great_circle_points,
    hmi_disambiguate_azimuth,
    hmi_vector_components,
    match_template_location,
    patch_bounds_from_lon_lat,
    pfss_amrvac_parameters,
    pfss_mod_usr_parameters,
    resample_great_circle_path,
    spherical_cap_flux,
    spherical_flux_balance,
    spherical_preprocess_metrics,
    spherical_to_cartesian,
    read_path_dat,
    preprocess_spherical_field,
    write_path_dat,
    write_cutout,
    write_pfss_boundary,
)


def test_cartesian_grid_shape_uses_zyx_data_order():
    grid = CartesianGrid(x=np.arange(4), y=np.arange(3), z=np.arange(2), unit="Mm")

    assert grid.shape == (2, 3, 4)
    assert grid.spacing == (1.0, 1.0, 1.0)


def test_local_cartesian_from_heliographic_uses_cea_boundary_convention():
    br = np.array([[3.0]])
    bt = np.array([[-2.0]])
    bp = np.array([[1.0]])

    bx, by, bz = local_cartesian_from_heliographic(br, bt, bp)

    assert np.allclose(bx, [[1.0]])
    assert np.allclose(by, [[2.0]])
    assert np.allclose(bz, [[3.0]])


def test_crop_and_multigrid_reduce_components_block_average():
    data = np.arange(24, dtype=float).reshape(4, 6)
    bx, by, bz = crop_components(data, data + 100, data + 200, x0=1, y0=0, nx=4, ny=4)

    assert bx.shape == (4, 4)
    bx2, by2, bz2 = multigrid_reduce_components(bx, by, bz, level=2)

    assert bx2.shape == (2, 2)
    assert np.allclose(bx2, [[4.5, 6.5], [16.5, 18.5]])
    assert np.allclose(by2, bx2 + 100)
    assert np.allclose(bz2, bx2 + 200)


def test_crop_components_with_padding_fills_out_of_frame_ghost_pixels():
    data = np.arange(20, dtype=float).reshape(4, 5)

    bx, by, bz, padding = crop_components_with_padding(
        data,
        data + 100,
        data + 200,
        x0=1,
        y0=0,
        nx=3,
        ny=4,
        pad_x=1,
        pad_y=2,
        fill_value=0.0,
    )

    assert bx.shape == (8, 5)
    assert padding["bottom"] == 2
    assert padding["top"] == 2
    assert np.allclose(bx[:2], 0.0)
    assert np.allclose(bx[-2:], 0.0)
    assert np.array_equal(bx[2:6], data[:, 0:5])
    assert np.array_equal(by[2:6], data[:, 0:5] + 100)
    assert np.array_equal(bz[2:6], data[:, 0:5] + 200)


def test_preprocess_cartesian_field_returns_finite_same_shape():
    rng = np.random.default_rng(42)
    bz = rng.normal(size=(8, 8))
    bx = rng.normal(scale=0.2, size=bz.shape)
    by = rng.normal(scale=0.2, size=bz.shape)

    bx_p, by_p, bz_p, metrics = preprocess_cartesian_field(
        bx,
        by,
        bz,
        max_iter=5,
        report_interval=100,
    )

    assert bx_p.shape == bx.shape
    assert by_p.shape == by.shape
    assert bz_p.shape == bz.shape
    assert np.isfinite(bx_p).all()
    assert np.isfinite(by_p).all()
    assert np.isfinite(bz_p).all()
    assert metrics["iteration"] == 5


def test_preprocess_smoothing_term_matches_shift_laplace_form():
    rng = np.random.default_rng(123)
    field = rng.normal(size=(7, 9))
    lap = _periodic_laplace(field)

    optimized = 2.0 * _periodic_laplace(lap)
    idl_form = 2.0 * (
        -4.0 * lap
        + _periodic_laplace(np.roll(field, 1, axis=1))
        + _periodic_laplace(np.roll(field, -1, axis=1))
        + _periodic_laplace(np.roll(field, 1, axis=0))
        + _periodic_laplace(np.roll(field, -1, axis=0))
    )

    assert np.allclose(optimized, idl_form)


def test_guo_cartesian_boundary_binary_headers_and_parameters(tmp_path):
    bx = np.arange(48, dtype=float).reshape(6, 8)
    by = bx + 100
    bz = bx + 200
    meta = boundary_metadata(nx=8, ny=6, xc_cm=10.0e9, yc_cm=-5.0e9, dx_cm=2.0e8, dy_cm=3.0e8)

    potential = write_potential_boundary(tmp_path / "potential_boundary.dat", bz, meta)
    nlfff = write_nlfff_boundary(tmp_path / "nlfff_boundary.dat", bx, by, bz, meta)

    with potential.open("rb") as handle:
        dims = np.fromfile(handle, dtype=np.int32, count=2)
        header = np.fromfile(handle, dtype=np.float64, count=4)
        data = np.fromfile(handle, dtype=np.float64)
    assert np.array_equal(dims, [4, 2])
    assert np.allclose(header, [10.0e9, -5.0e9, 2.0e8, 3.0e8])
    assert np.array_equal(data, bz[2:-2, 2:-2].ravel(order="C"))

    with nlfff.open("rb") as handle:
        dims = np.fromfile(handle, dtype=np.int32, count=2)
        header = np.fromfile(handle, dtype=np.float64, count=4)
        data = np.fromfile(handle, dtype=np.float64)
    assert np.array_equal(dims, [8, 6])
    assert np.allclose(header, [10.0e9, -5.0e9, 2.0e8, 3.0e8])
    assert np.array_equal(data[: bx.size], bx.ravel(order="C"))

    text = format_guo_amrvac_parameters(meta)
    assert "domain_nx1=4" in text
    assert "domain_nx2=2" in text
    assert "xprobmin3=0.1d0" in text


def test_spherical_grid_shape_uses_r_theta_phi_order():
    grid = SphericalGrid(r=np.arange(2), theta=np.arange(3), phi=np.arange(4))

    assert grid.shape == (2, 3, 4)


def test_spherical_grid_latitude_from_colatitude():
    grid = SphericalGrid(r=np.ones(1), theta=np.array([0.0, np.pi / 2.0, np.pi]), phi=np.ones(1))

    assert np.allclose(grid.latitude, [np.pi / 2.0, 0.0, -np.pi / 2.0])


def test_cea_grid_centers_and_equatorial_lon_lat_mapping():
    x, y = cea_xy_centers(3, 3, 0.5)

    assert np.allclose(x, [-0.5, 0.0, 0.5])
    assert np.allclose(y, [-0.5, 0.0, 0.5])

    lon, lat = cea_lon_lat_from_xy(x, y, 10.0, 0.0)

    assert np.allclose(lon, [9.5, 10.0, 10.5])
    assert np.allclose(lat, np.degrees(np.arcsin(np.radians(y))))


def test_make_cea_patch_grid_shape_uses_y_x_order():
    grid = make_cea_patch_grid(10.0, 20.0, width_degree=1.0, height_degree=0.5, resolution_degree=0.25)

    assert grid.shape == (2, 4)
    assert grid.lon.shape == (2, 4)
    assert grid.lat.shape == (2, 4)


def test_hmi_native_vector_components_matches_generated_sharp_comparison():
    bxi, beta, bzeta = hmi_native_vector_components(
        field=np.array([2.0]),
        inclination=np.array([90.0]),
        azimuth=np.array([90.0]),
    )

    assert np.allclose(bxi, [2.0])
    assert np.allclose(beta, [0.0], atol=1e-12)
    assert np.allclose(bzeta, [0.0], atol=1e-12)


def test_gaussian_block_reduce_preserves_constant_field():
    data = np.full((6, 9), 5.0)

    reduced = gaussian_block_reduce(data, 2, 3, sigma_y=0.5, sigma_x=0.5)

    assert reduced.shape == (3, 3)
    assert np.allclose(reduced, 5.0)


def test_gaussian_block_reduce_handles_nan_weights():
    data = np.full((4, 4), 2.0)
    data[0, 0] = np.nan

    reduced = gaussian_block_reduce(data, 2, 2, sigma_y=0.0, sigma_x=0.0)

    assert reduced.shape == (2, 2)
    assert np.isfinite(reduced).all()
    assert np.allclose(reduced, 2.0)


def test_native_to_heliographic_components_at_disk_center():
    br, bt, bp = native_to_heliographic_components(
        bxi=np.array([[1.0]]),
        beta=np.array([[2.0]]),
        bzeta=np.array([[3.0]]),
        lon_degree=np.array([[0.0]]),
        lat_degree=np.array([[0.0]]),
        disk_lon_degree=0.0,
        disk_lat_degree=0.0,
        p_angle_degree=0.0,
    )

    assert np.allclose(br, [[3.0]])
    assert np.allclose(bt, [[-2.0]])
    assert np.allclose(bp, [[1.0]])


def test_heliographic_to_local_cartesian_boundary_uses_northward_y():
    bx, by, bz = heliographic_to_local_cartesian_boundary(
        br=np.array([[3.0]]),
        bt=np.array([[-2.0]]),
        bp=np.array([[1.0]]),
    )

    assert np.allclose(bx, [[1.0]])
    assert np.allclose(by, [[2.0]])
    assert np.allclose(bz, [[3.0]])


def test_spherical_to_cartesian_equator_x_axis():
    x, y, z = spherical_to_cartesian(2.0, np.pi / 2.0, 0.0)

    assert np.isclose(x, 2.0)
    assert np.isclose(y, 0.0)
    assert np.isclose(z, 0.0)


def test_angular_distance_and_great_circle_points():
    assert np.isclose(angular_distance(0.0, 0.0, 90.0, 0.0), 90.0)

    points = great_circle_points(0.0, 0.0, 90.0, 0.0, 3)

    assert np.allclose(points[:, 0], [0.0, 45.0, 90.0])
    assert np.allclose(points[:, 1], [0.0, 0.0, 0.0], atol=1e-12)


def test_resample_great_circle_path_returns_equal_distance_nodes():
    path = np.array([[0.0, 0.0], [90.0, 0.0], [180.0, 0.0]])

    points, total = resample_great_circle_path(path, 5)

    assert np.isclose(total, 180.0)
    assert np.allclose(points[:, 0], [0.0, 45.0, 90.0, 135.0, 180.0])
    assert np.allclose(points[:, 1], 0.0, atol=1e-12)


def test_spherical_cap_flux_uses_standard_angular_distance():
    lon = np.array([-1.0, 0.0, 1.0])
    lat = np.array([-1.0, 0.0, 1.0])
    br = np.full((3, 3), 2.0)

    result = spherical_cap_flux(lon, lat, br, 0.0, 0.0, radius=10.0, solar_radius=696.0)

    assert result["pixel_count"] == 1
    assert np.isclose(result["mean_br"], 2.0)
    assert np.isclose(result["flux_mx"], 2.0 * np.pi * (10.0e8) ** 2)


def test_flux_rope_axis_from_points_and_path_dat_roundtrip(tmp_path):
    points = np.array([[30.0, 20.0], [35.0, 21.0], [40.0, 20.5]])

    result = flux_rope_axis_from_points(points, height_mm=20.0, npoints=11, solar_radius_mm=696.0)
    axis = result["path"]

    assert axis.shape == (11, 3)
    assert np.allclose(axis[0, :2], points[0])
    assert np.allclose(axis[-1, :2], points[-1])
    assert np.isclose(axis[0, 2], 1.0)
    assert np.isclose(axis[-1, 2], 1.0)
    assert axis[:, 2].max() > 1.0
    assert result["depth_mm"] > 0.0
    assert result["major_radius_mm"] > 20.0

    path = write_path_dat(tmp_path / "path.dat", axis)
    raw = read_path_dat(path)
    restored = read_path_dat(path, subtract_lon_180=True)

    assert raw.shape == (11, 3)
    assert np.allclose(raw[:, 0], axis[:, 0] + 180.0)
    assert np.allclose(restored, axis)


def test_preprocess_spherical_field_reduces_force_torque_metric():
    rng = np.random.default_rng(1234)
    theta = np.linspace(0.4, 1.2, 12)
    phi = np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False)
    br = rng.normal(size=(theta.size, phi.size))
    bt = rng.normal(scale=0.4, size=br.shape)
    bp = rng.normal(scale=0.4, size=br.shape)

    before = spherical_preprocess_metrics(br, bt, bp, theta, phi)
    brn, btn, bpn = preprocess_spherical_field(
        br,
        bt,
        bp,
        theta,
        phi,
        mu3=0.1,
        mu4=1.0,
        max_iter=60,
        dt=1000,
    )
    after = spherical_preprocess_metrics(brn, btn, bpn, theta, phi)

    assert brn.shape == br.shape
    assert np.isfinite(brn).all()
    assert np.isfinite(btn).all()
    assert np.isfinite(bpn).all()
    assert after["l12_force_torque"] < before["l12_force_torque"]


def test_hmi_disambiguate_azimuth_applies_selected_bit():
    azimuth = np.array([[10.0, 20.0], [30.0, 40.0]])
    disambig = np.array([[0, 2], [4, 6]])

    corrected = hmi_disambiguate_azimuth(azimuth, disambig, method=1)

    assert np.allclose(corrected, [[10.0, 200.0], [30.0, 220.0]])


def test_hmi_vector_components_matches_notebook_convention():
    bx, by, bz = hmi_vector_components(
        field=np.array([2.0]),
        inclination=np.array([90.0]),
        azimuth=np.array([90.0]),
    )

    assert np.allclose(bx, [2.0])
    assert np.allclose(by, [0.0], atol=1e-12)
    assert np.allclose(bz, [0.0], atol=1e-12)


def test_spherical_flux_balance_uses_cos_latitude_weight():
    br = np.array([[2.0, -1.0], [1.0, -1.0]])
    flux = spherical_flux_balance(br, np.array([0.0, 60.0]))

    assert np.isclose(flux["positive"], 2.5)
    assert np.isclose(flux["negative"], -1.5)
    assert np.isclose(flux["signed"], 1.0)
    assert np.isclose(flux["unsigned"], 4.0)
    assert np.isclose(flux["balance"], 0.25)


def test_pfss_parameter_helpers_match_notebook_indexing():
    x0, y0, nx, ny = patch_bounds_from_lon_lat((-7.0, 65.0), (0.0, 40.0), 0.05)

    assert (x0, y0, nx, ny) == (3460, 1800, 1440, 800)
    assert pfss_mod_usr_parameters((3600, 7200), x0, y0, nx, ny) == {
        "lat0": 1000,
        "lat1": 1799,
        "lon0": 3460,
        "lon1": 4899,
    }

    amrvac = pfss_amrvac_parameters((3600, 7200), x0, y0, nx, ny, nr=400, length_r=0.8)
    assert amrvac["domain_nx1"] == 800
    assert amrvac["domain_nx2"] == 800
    assert amrvac["domain_nx3"] == 1440
    assert np.isclose(amrvac["xprobmin1"], 1.001)
    assert np.isclose(amrvac["xprobmax1"], 1.801)


def test_match_template_location_finds_inserted_patch():
    big = np.zeros((8, 10), dtype=float)
    template = np.array([[1.0, 2.0, 3.0], [0.0, 4.0, 5.0]])
    big[4:6, 3:6] = template

    x, y, score = match_template_location(big, template)

    assert (x, y) == (3, 4)
    assert score > 0.99


def test_write_pfss_boundary_and_cutout_binary_headers(tmp_path):
    pfss = np.arange(6, dtype=float).reshape(2, 3)
    pfss_path = write_pfss_boundary(tmp_path / "pfss_boundary.dat", pfss, [-179.5, -178.5, -177.5], [-89.5, -88.5])
    with pfss_path.open("rb") as handle:
        dims = np.fromfile(handle, dtype=np.int32, count=2)
        theta = np.fromfile(handle, dtype=np.float64, count=2)
        phi = np.fromfile(handle, dtype=np.float64, count=3)
        data = np.fromfile(handle, dtype=np.float64, count=6)

    assert np.array_equal(dims, [3, 2])
    assert np.allclose(theta, np.pi / 2 - np.radians([-89.5, -88.5]))
    assert np.allclose(phi, np.pi + np.radians([-179.5, -178.5, -177.5]))
    assert np.array_equal(data, pfss.ravel())

    cutout_path = write_cutout(tmp_path / "cutout.dat", pfss, pfss + 10, pfss + 20)
    with cutout_path.open("rb") as handle:
        dims = np.fromfile(handle, dtype=np.int32, count=2)
        br = np.fromfile(handle, dtype=np.float64, count=6)

    assert np.array_equal(dims, [2, 3])
    assert np.array_equal(br, np.rot90(pfss, k=-1).ravel())

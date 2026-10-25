"""Regression coverage for CEA reference time, pixel geometry and saved AIA."""
from pathlib import Path
import astropy.units as u
import numpy as np
import pytest
from sunpy.map import Map

from Packages.sharp.cea_pipeline import build_target, Coordinates, write_components
from workflows.observations.project_aia_to_sharp import project_aia_product, parse_args


@pytest.fixture
def hmi():
    # DATE-OBS intentionally differs from T_OBS, as in real HMI FITS.
    return Map(np.zeros((8,8)), dict(
        instrume='HMI_COMBINED', telescop='SDO/HMI', bunit='G',
        ctype1='HPLN-TAN', ctype2='HPLT-TAN', cunit1='arcsec', cunit2='arcsec',
        crpix1=2040.42334, crpix2=2050.6958, crval1=0., crval2=0.,
        cdelt1=.504039764, cdelt2=.504039764, crota2=180.013489,
        dsun_obs=149428849142.56348, rsun_ref=696000000., rsun_obs=960.7302856,
        crlt_obs=6.27475405, crln_obs=103.91481,
        **{'date-obs':'2024-10-08T23:58:37.400',
           't_obs':'2024.10.08_23:59:59.436_TAI','t_rec':'2024.10.09_00:00:00_TAI'}))


def test_written_magnetic_reference_time_and_history_survive(tmp_path,hmi):
    target=build_target(hmi,center=(30.,25.),width=.21,height=.15)
    write_components(tmp_path,target,{'Br':np.ones(target.data.shape)})
    saved=Map(tmp_path/'Br.fits')
    assert abs((saved.reference_date-hmi.reference_date).to_value(u.s))<1e-6
    assert abs((saved.reference_date-saved.date).to_value(u.s))>40
    assert saved.wcs.wcs.compare(target.wcs.wcs,tolerance=1e-10)


def test_official_geometry_matches_wcs_with_shifted_crpix(hmi):
    if not Path(__file__).parents[1].joinpath('Packages/sharp/build/libcustom_coordinates.dylib').exists():
        pytest.skip('Build official numerical core to test compiled geometry.')
    target=build_target(hmi,center=(30.,25.),width=.51,height=.33)
    target=target.submap((2,1)*u.pix,top_right=(11,7)*u.pix)
    _,_,lon,lat=Coordinates(hmi,target)(0,target.data.shape[0])
    yy,xx=np.indices(target.data.shape)
    world=target.pixel_to_world(xx*u.pix,yy*u.pix)
    assert np.max(np.abs(lon*180/np.pi-world.lon.deg))<1e-5
    assert np.max(np.abs(lat*180/np.pi-world.lat.deg))<1e-5


def test_aia_odd_grid_mask_units_time_and_fits_roundtrip(tmp_path,hmi):
    target=build_target(hmi,center=(30.,25.),width=.21,height=.15)
    target.data[0,0]=np.nan
    header=target.wcs.to_header()
    header.update(INSTRUME='AIA',TELESCOP='SDO/AIA',WAVELNTH=193,WAVEUNIT='angstrom',BUNIT='DN',EXPTIME=2.)
    header['DATE-OBS']=target.reference_date.utc.isot;header['T_OBS']=target.reference_date.utc.isot
    aia=Map(np.full(target.data.shape,8.),header)
    output=project_aia_product(aia,target)
    output.save(tmp_path/'aia.fits',overwrite=True)
    saved=Map(tmp_path/'aia.fits')
    assert saved.data.shape==(5,7)
    assert np.isnan(saved.data[0,0])
    assert np.allclose(saved.data[np.isfinite(target.data)],4.)
    assert saved.unit.is_equivalent(u.DN/u.s)
    assert abs((saved.reference_date-target.reference_date).to_value(u.s))<1e-6
    assert saved.wcs.wcs.compare(output.wcs.wcs,tolerance=1e-10)
    # Reprojecting an already normalized product must not divide by exposure again.
    again=project_aia_product(output,target)
    assert np.nanmedian(again.data)==pytest.approx(4.)


def test_hpc_grid_rejects_off_disk_bounds(hmi):
    with pytest.raises(ValueError,match='visible solar'):
        build_target(hmi,bounds=[1200,1300,100,200])


def test_aia_cli_full_target_and_optional_bin2():
    args=parse_args(['--sharp','Br.fits','--aia','aia.fits','--out','out'])
    assert args.x0 is None and not args.bin2 and not args.raw_dn

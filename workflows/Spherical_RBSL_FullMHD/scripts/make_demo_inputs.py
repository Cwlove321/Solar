#!/usr/bin/env python3
"""Generate small synthetic configuration fixtures, not observational data."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import struct

import numpy as np
from astropy.io import fits


def make_demo(output: Path) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite fixture directory: {output}')
    boundary = output / 'boundary'
    axis = output / 'axis'
    boundary.mkdir(parents=True)
    axis.mkdir()
    nphi, ntheta, resolution = 360, 180, 1.0
    nlon, nlat, x0, y0 = 24, 24, 168, 90
    lon_range, lat_range = [-12., 12.], [0., 24.]
    lon = np.arange(nlon) + .5
    lat = np.arange(nlat) + .5
    xx, yy = np.meshgrid(lon, lat)
    patch = 100.0 * (np.exp(-((xx - 7.)**2 + (yy - 12.)**2)/12.)
                     - np.exp(-((xx - 17.)**2 + (yy - 12.)**2)/12.))
    global_br = np.zeros((ntheta, nphi), dtype='<f8')
    global_br[y0:y0+nlat, x0:x0+nlon] = patch
    theta = np.radians(180 - (np.arange(ntheta) + .5) * resolution)
    phi = np.radians((np.arange(nphi) + .5) * resolution)
    with (boundary / 'pfss_boundary.dat').open('wb') as stream:
        stream.write(struct.pack('<ii', nphi, ntheta))
        for values in (theta, phi, global_br):
            np.asarray(values, dtype='<f8').tofile(stream)
    header = fits.Header({'BUNIT': 'G', 'NGHOST': 0, 'LONMIN': lon_range[0],
                          'LONMAX': lon_range[1], 'LATMIN': lat_range[0], 'LATMAX': lat_range[1],
                          'CONTENT': 'SYNTHETIC CONFIGURATION TEST ONLY'})
    fits.PrimaryHDU(patch[::-1, :].T.astype(np.float32), header=header).writeto(boundary / 'Brll0000.fits')
    metadata = {'schema_version': 1, 'synthetic': True,
                'coordinate_convention': 'Stonyhurst longitude + 180 degrees in simulation phi',
                'grid': {'resolution_degree': resolution, 'effective_nghost': 0,
                         'base_nghost': 0, 'refine_level': 1,
                         'nlon': nlon, 'nlat': nlat, 'x0': x0, 'y0': y0,
                         'lon_range': lon_range, 'lat_range': lat_range}}
    (boundary / 'parameters.json').write_text(json.dumps(metadata, indent=2) + '\n')
    s = np.linspace(0., 1., 32)
    path = np.column_stack((171 + 18*s, 12 + 2*np.sin(np.pi*s), 1 + .06*np.sin(np.pi*s)))
    np.asarray(path, dtype='<f8').tofile(axis / 'path.dat')
    (axis / 'path_summary.json').write_text(json.dumps({'npoints': len(path), 'synthetic': True,
        'coordinate_convention': 'longitude+180 deg, latitude deg, radius/Rsun'}, indent=2) + '\n')
    settings = {'pfss_stage': '01_pfss', 'fullmhd_source': 'insert',
                'thermal_profile': 'default', 'amr_mode': 'fixed', 'atmosphere': {},
                'grid': {'domain_nx1': 12, 'xprobmin1': 1., 'xprobmax1': 1.5,
                         'refine_max_level': 2, 'qstretch_baselevel1': 1.04},
                'rbsl': {'bflux_mx': 1.e20, 'radius_mm': 10., 'positive_helicity': False, 'flip_axis': False},
                'global_pfss': {'lmax': 90, 'r0': 1., 'rs': 2.5},
                'filament': {'method': 'none'},
                'stage_times': {'02_relax': {'time_max_code': .001, 'save_interval_code': .001},
                                '03_insert': {'time_max_code': 0., 'save_interval_code': .001},
                                '04_fullmhd': {'time_max_code': .001, 'save_interval_code': .001}}}
    settings_path = output / 'settings.json'
    settings_path.write_text(json.dumps(settings, indent=2) + '\n')
    (output / 'README.txt').write_text('Synthetic arrays for configuration and pipeline tests only. These are not solar observations or a physically validated MHD benchmark.\n')
    return {'boundary': str(boundary), 'axis': str(axis / 'path.dat'), 'settings': str(settings_path)}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    print(json.dumps(make_demo(args.output), indent=2))

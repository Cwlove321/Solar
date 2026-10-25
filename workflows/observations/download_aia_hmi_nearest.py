#!/usr/bin/env python3
"""Download nearest good AIA/HMI FITS and plot a common-time solar-coordinate pair.

Run from solar_code using /opt/anaconda3/bin/python. JSOC email is passed via
--email. Raw exported FITS are preserved; only the preview is reprojected.
"""
from __future__ import annotations
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import socket
import requests
import time

os.environ.setdefault('SUNPY_CONFIGDIR', '/private/tmp/solarcode_sunpy')
os.environ.setdefault('MPLCONFIGDIR', '/private/tmp/solarcode_mpl')
import astropy.units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.time import Time
from astropy.visualization import ImageNormalize, AsinhStretch
import drms
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
from sunpy.coordinates import HeliographicStonyhurst, propagate_with_solar_surface
from sunpy.map import Map, make_fitswcs_header
from sunpy.time import parse_time


def jsoc_time(value):
    if str(value).endswith('_TAI'):
        return Time(str(value).removesuffix('_TAI').replace('.', '-', 2).replace('_', 'T'), scale='tai')
    return parse_time(value)


def save_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + '\n')


def download_ranges(url, destination, workers=8):
    """Fetch verified byte ranges concurrently on high-latency archive links."""
    def session():
        s = requests.Session()
        s.trust_env = False
        return s
    with session() as client:
        head = client.head(url, timeout=60)
        head.raise_for_status()
        total = int(head.headers['Content-Length'])
    width = (total + workers - 1)//workers
    def fetch(i):
        lo, hi = i*width, min((i+1)*width, total)-1
        path = destination.with_name(destination.name + f'.chunk{i}')
        expected = hi-lo+1
        if path.exists() and path.stat().st_size == expected:
            return path
        for attempt in range(3):
            try:
                with session() as client:
                    with client.get(url, headers={'Range': f'bytes={lo}-{hi}'}, stream=True, timeout=(30,120)) as response:
                        response.raise_for_status()
                        if response.status_code != 206 or response.headers.get('Content-Range') != f'bytes {lo}-{hi}/{total}':
                            raise ValueError('Server returned an unexpected byte range')
                        with path.open('wb') as f:
                            for chunk in response.iter_content(262144):
                                f.write(chunk)
                if path.stat().st_size != expected:
                    raise ValueError('Incomplete byte range')
                print(f'Download part {i+1}/{workers}: {destination.name}', flush=True)
                return path
            except (requests.RequestException, ValueError):
                if attempt == 2:
                    raise
    with ThreadPoolExecutor(max_workers=workers) as pool:
        parts = list(pool.map(fetch, range(workers)))
    with destination.open('wb') as f:
        for path in parts:
            f.write(path.read_bytes())
    if destination.stat().st_size != total:
        raise ValueError('Assembled file size mismatch')
    for path in parts:
        path.unlink()


def acquire(kind, series, segment, out, args):
    out.mkdir(parents=True, exist_ok=True)
    state_path = out / f'nearest_{args.tag}.json'
    if state_path.exists():
        state = json.loads(state_path.read_text())
        if Path(state.get('file', '')).is_file():
            with fits.open(state['file']) as hdul:
                hdul.verify('exception')
            actual = Map(state['file']).reference_date.utc
            state['query_t_obs_utc'] = state.get('query_t_obs_utc', state['t_obs_utc'])
            state['t_obs_utc'] = actual.isot
            state['offset_seconds'] = float((actual-Time(args.time, scale='utc')).to_value(u.s))
            save_json(state_path, state)
            print(f'Reusing {state["file"]}', flush=True)
            return state
    else:
        state = {}
    target = Time(args.time, scale='utc')
    if 'url' in state:
        url = state['url']
    else:
        c = drms.Client(email=args.email)
        target = Time(args.time, scale='utc')
        start = (target - 90*u.s).utc.isot + 'Z'
        ds = f'{series}[{start}/180s]' + (f'[{args.wavelength}]' if kind == 'aia' else '')
        if 'request_id' not in state:
            candidates = c.query(ds, key='T_REC,T_OBS,DATE__OBS,QUALITY', rec_index=True)
            candidates = candidates[candidates.QUALITY == 0].copy()
            if candidates.empty:
                raise RuntimeError(f'No good records in {ds}')
            candidates['offset_s'] = [(jsoc_time(t)-target).to_value(u.s) for t in candidates.T_OBS]
            row = candidates.loc[candidates.offset_s.abs().idxmin()]
            record = row.name + '{' + segment + '}'
            print(f'{kind}: selected {record}; T_OBS={jsoc_time(row.T_OBS).utc.isot} UTC', flush=True)
            request = c.export(record, method='url', protocol='fits')
            state = dict(kind=kind, series=series, record=record,
                         query=ds, selection='nearest QUALITY=0 by T_OBS exposure midpoint',
                         reference_time_utc=target.utc.isot,
                         t_obs_utc=jsoc_time(row.T_OBS).utc.isot,
                         date_obs_utc=row.DATE__OBS, t_rec=row.T_REC,
                         offset_seconds=float(row.offset_s), quality=int(row.QUALITY),
                         request_id=request.id)
            save_json(state_path, state)
        else:
            request = c.export_from_id(state['request_id'])
        deadline = time.monotonic() + 1200
        while not request.wait(timeout=30, sleep=5):
            if time.monotonic() > deadline:
                raise TimeoutError(f'Export pending: {request.id}; rerun to resume')
            print(f'{kind}: export pending {request.id}', flush=True)
        urls = request.urls
        if len(urls) != 1:
            raise RuntimeError(f'Expected one file, got {len(urls)}')
        url = str(urls.iloc[0].url).replace('http://', 'https://', 1)
        url = url.replace('https://jsoc.stanford.edu/', 'https://jsoc1.stanford.edu/')
        state['url'] = url
        save_json(state_path, state)
    dest = out / Path(url.split('?')[0]).name
    partial = dest.with_suffix(dest.suffix + '.part')
    print(f'{kind}: downloading {url}', flush=True)
    download_ranges(url, partial)
    with fits.open(partial) as hdul:
        hdul.verify('exception')
        data_hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
        if data_hdu.data.shape != (4096, 4096):
            raise ValueError(f'Unexpected full disk shape: {data_hdu.data.shape}')
    m = Map(partial)
    if kind == 'aia' and int(m.wavelength.to_value(u.AA)) != args.wavelength:
        raise ValueError('Unexpected AIA wavelength')
    actual = m.reference_date.utc
    if abs((actual-target).to_value(u.s)) > 90:
        raise ValueError('FITS observation time outside search window')
    partial.rename(dest)
    state.update(file=str(dest), url=url, bytes=dest.stat().st_size,
                 sha256=hashlib.sha256(dest.read_bytes()).hexdigest(),
                 fits_reference_time_utc=actual.isot, shape=list(m.data.shape),
                 query_t_obs_utc=state['t_obs_utc'], t_obs_utc=actual.isot,
                 offset_seconds=float((actual-target).to_value(u.s)))
    save_json(state_path, state)
    return state


def plot_pair(a, h, args):
    am, hm = Map(a['file']), Map(h['file'])
    target = Time(args.time, scale='utc')
    xmin, xmax, ymin, ymax = args.extent
    nx, ny = int(round((xmax-xmin)/args.scale)), int(round((ymax-ymin)/args.scale))
    if nx <= 0 or ny <= 0:
        raise ValueError('Invalid extent or pixel scale')
    observer = am.observer_coordinate.transform_to(HeliographicStonyhurst(obstime=target))
    center = SkyCoord((xmin+xmax)/2*u.arcsec, (ymin+ymax)/2*u.arcsec,
                      frame='helioprojective', observer=observer, obstime=target, rsun=am.rsun_meters)
    header = make_fitswcs_header((ny, nx), center,
              scale=[(xmax-xmin)/nx, (ymax-ymin)/ny]*u.arcsec/u.pix,
              projection_code='TAN')
    # Use identical solar radius in both coordinate transforms.
    hm.meta['rsun_ref'] = am.rsun_meters.to_value(u.m)
    with propagate_with_solar_surface():
        ap = am.reproject_to(header)
        hp = hm.reproject_to(header)
    intensity = ap.data/am.exposure_time.to_value(u.s)
    b_los = hp.data
    valid_a, valid_h = float(np.isfinite(intensity).mean()), float(np.isfinite(b_los).mean())
    if min(valid_a, valid_h) < 0.99:
        raise ValueError(f'Insufficient coverage: {valid_a}, {valid_h}')
    vmax = float(args.aia_vmax) if args.aia_vmax is not None else float(np.nanpercentile(intensity, 99.7))
    colors_q_d = ['#98D9E4', '#49ABC1', '#F4DF86', '#ED9B62', '#C64146']
    cmap = (LinearSegmentedColormap.from_list('white_cyan_yellow_red_D', colors_q_d, N=256)
            if args.aia_colormap == 'q-d' else am.plot_settings['cmap'])
    if args.aia_vmin >= vmax:
        raise ValueError('AIA display minimum must be less than maximum')
    norm = ImageNormalize(vmin=args.aia_vmin, vmax=vmax, stretch=AsinhStretch(0.03), clip=True)
    plt.rcParams.update({'font.size': 11, 'font.family':'DejaVu Sans', 'axes.titlesize':12})
    fig, axes = plt.subplots(1, 2, figsize=(12.7, 6.6), sharex=True, sharey=True)
    titles = [f'AIA {args.wavelength} Å | {a["t_obs_utc"].split("T")[1]} UTC',
              f'HMI LOS (45 s) | {h["t_obs_utc"].split("T")[1]} UTC']
    im_a = axes[0].imshow(intensity, origin='lower', extent=args.extent,
                          cmap=cmap, norm=norm, interpolation='nearest')
    im_h = axes[1].imshow(b_los, origin='lower', extent=args.extent,
                          cmap='gray', vmin=-args.hmi_limit, vmax=args.hmi_limit, interpolation='nearest')
    overlay_alpha = None
    if args.overlay_threshold is not None:
        overlay_alpha = np.where(np.isfinite(intensity) & (intensity >= args.overlay_threshold), args.overlay_alpha, 0.0)
        assert np.all(overlay_alpha[intensity < args.overlay_threshold] == 0)
        assert np.all(overlay_alpha[intensity >= args.overlay_threshold] == args.overlay_alpha)
        axes[1].imshow(intensity, origin='lower', extent=args.extent,
                       cmap=cmap, norm=norm,
                       alpha=overlay_alpha, interpolation='nearest')
        titles[1] = f'HMI LOS + AIA {args.wavelength} Å (≥{args.overlay_threshold:g} DN/s; α={args.overlay_alpha:g})'
    for ax, title in zip(axes, titles):
        ax.set(xlim=(xmin,xmax), ylim=(ymin,ymax), title=title, xlabel='Solar X [arcsec]', aspect='equal')
        ax.set_xticks(np.arange(np.ceil(xmin/100)*100, xmax+1, 100))
        ax.set_yticks(np.arange(np.ceil(ymin/100)*100, ymax+1, 100))
        ax.tick_params(direction='out')
    axes[0].set_ylabel('Solar Y [arcsec]')
    for ax, im, label in [(axes[0],im_a,f'AIA {args.wavelength} Å intensity [DN s⁻¹]'), (axes[1],im_h,'LOS magnetic field [G]')]:
        bar=fig.colorbar(im, ax=ax, orientation='horizontal', fraction=0.047, pad=0.12, extend='both' if ax is axes[1] else 'max')
        if ax is axes[0]:
            ticks = args.aia_vmin + (vmax-args.aia_vmin) * np.array([0, 0.025, 0.125, 0.5, 1.0])
            bar.set_ticks(ticks)
            bar.set_ticklabels([f'{value:.0f}' for value in ticks])
        bar.set_label(label)
    fig.suptitle(f'SDO observations · {target.utc.isot[:10]}', y=0.98, fontsize=15)
    fig.text(0.5,0.92, f'Common reference: {target.utc.isot.split("T")[1]} UTC   |   X [{xmin:g}, {xmax:g}]″   Y [{ymin:g}, {ymax:g}]″', ha='center', fontsize=10)
    fig.subplots_adjust(left=0.07,right=0.98,bottom=0.12,top=0.86,wspace=0.12)
    args.figure_dir.mkdir(parents=True,exist_ok=True)
    stem=f'aia{args.wavelength}_hmi45_{args.tag}_x{xmin:g}_{xmax:g}_y{ymin:g}_{ymax:g}'
    if overlay_alpha is not None:
        stem += f'_overlay_ge{args.overlay_threshold:g}_alpha{args.overlay_alpha:g}'
    if args.aia_colormap != 'aia':
        stem += f'_cmap_{args.aia_colormap}_vmin{args.aia_vmin:g}'
    if args.aia_vmax is not None:
        stem += f'_vmax{args.aia_vmax:g}'
    outputs=[]
    for suffix in ('png','pdf'):
        path=args.figure_dir/f'{stem}.{suffix}'
        fig.savefig(path,dpi=220,bbox_inches='tight',facecolor='white')
        outputs.append(str(path))
    plt.close(fig)
    if overlay_alpha is not None:
        fig, ax = plt.subplots(figsize=(8.8, 8.8))
        base = ax.imshow(b_los, origin='lower', extent=args.extent,
                         cmap='gray', vmin=-args.hmi_limit, vmax=args.hmi_limit,
                         interpolation='nearest')
        ax.imshow(intensity, origin='lower', extent=args.extent,
                  cmap=cmap, norm=norm,
                  alpha=overlay_alpha, interpolation='nearest')
        ax.set(xlim=(xmin,xmax), ylim=(ymin,ymax), aspect='equal',
               xlabel='Solar X [arcsec]', ylabel='Solar Y [arcsec]')
        ax.set_xticks(np.arange(np.ceil(xmin/100)*100, xmax+1, 100))
        ax.set_yticks(np.arange(np.ceil(ymin/100)*100, ymax+1, 100))
        ax.set_title(f'HMI LOS + AIA {args.wavelength} Å | {target.utc.isot[:10]}', fontsize=15, pad=42)
        ax.text(0.5, 1.035, f'Common reference: {target.utc.isot.split("T")[1]} UTC\nAIA ≥ {args.overlay_threshold:g} DN/s: α = {args.overlay_alpha:g}; below threshold: transparent',
                transform=ax.transAxes, ha='center', va='bottom', fontsize=10)
        cb=fig.colorbar(base, ax=ax, orientation='horizontal', fraction=0.04, pad=0.10, extend='both')
        cb.set_label('HMI LOS magnetic field [G]')
        fig.subplots_adjust(left=0.11, right=0.97, top=0.85, bottom=0.12)
        for suffix in ('png','pdf'):
            path=args.figure_dir/f'{stem}_single.{suffix}'
            fig.savefig(path, dpi=220, bbox_inches='tight', facecolor='white')
            outputs.append(str(path))
        plt.close(fig)
    manifest=dict(reference_time_utc=target.utc.isot, extent_arcsec=args.extent,
       grid_shape=[ny,nx], pixel_scale_arcsec=args.scale, aia=a, hmi=h,
       processing='Raw FITS preserved. Bilinear WCS reprojection with solar surface differential rotation to target UTC; common AIA observer and solar radius. AIA exposure normalized; no PSF deconvolution, degradation correction, or external pointing refinement. TAN-grid display extent approximates helioprojective arcsec to subpixel precision.',
       aia_display=dict(stretch='asinh',a=0.03,vmin=args.aia_vmin,vmax=vmax,unit='DN/s',
                        colormap=args.aia_colormap, colors=colors_q_d if args.aia_colormap == 'q-d' else None),
       hmi_display=dict(vmin=-args.hmi_limit,vmax=args.hmi_limit,unit='G'),
       finite_fraction=dict(aia=valid_a,hmi=valid_h), outputs=outputs,
       source_documentation='https://docs.sunpy.org/en/stable/generated/api/sunpy.coordinates.propagate_with_solar_surface.html')
    if overlay_alpha is not None:
        manifest['aia_overlay'] = dict(threshold_dn_per_s=args.overlay_threshold,
             alpha_at_or_above_threshold=args.overlay_alpha, alpha_below_threshold=0,
             visible_pixels=int(np.count_nonzero(overlay_alpha)),
             alpha_values=np.unique(overlay_alpha).tolist(),
             rule='Threshold applied to exposure-normalized linear intensity before color normalization')
    save_json(args.figure_dir/f'{stem}.json',manifest)
    print(json.dumps(manifest,indent=2),flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--time', required=True, help='Target observation time in UTC')
    p.add_argument('--email',required=True)
    p.add_argument('--wavelength',type=int,default=131)
    p.add_argument('--aia-dir',type=Path,required=True)
    p.add_argument('--hmi-dir',type=Path,required=True)
    p.add_argument('--figure-dir',type=Path,required=True)
    p.add_argument('--extent', nargs=4, type=float, required=True, metavar=('XMIN', 'XMAX', 'YMIN', 'YMAX'), help='Field of view in arcsec')
    p.add_argument('--scale',type=float,default=0.5)
    p.add_argument('--hmi-limit',type=float,default=800)
    p.add_argument('--aia-colormap', choices=['aia', 'q-d'], default='aia')
    p.add_argument('--aia-vmin', type=float, default=0.0)
    p.add_argument('--aia-vmax', type=float, default=None, help='Fixed AIA color scale maximum in DN/s')
    p.add_argument('--overlay-threshold', type=float, default=None, help='AIA DN/s threshold for overlay on HMI')
    p.add_argument('--overlay-alpha', type=float, default=0.5)
    args=p.parse_args()
    if not 0 <= args.overlay_alpha <= 1:
        p.error('--overlay-alpha must be between 0 and 1')
    args.tag=Time(args.time,scale='utc').strftime('%Y%m%d_%H%M%S_UTC')
    socket.setdefaulttimeout(60)
    with ThreadPoolExecutor(max_workers=2) as pool:
        fa=pool.submit(acquire,'aia','aia.lev1_euv_12s','image',args.aia_dir,args)
        fh=pool.submit(acquire,'hmi','hmi.M_45s','magnetogram',args.hmi_dir,args)
        a,h=fa.result(),fh.result()
    plot_pair(a,h,args)

if __name__=='__main__':
    main()

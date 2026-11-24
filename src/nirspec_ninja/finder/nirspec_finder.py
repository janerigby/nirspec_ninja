#!/usr/bin/env python3
"""NIRSpec overview: WISE W1 with four MSA quadrants, IFU, and fixed slits.

Requires astropy, numpy, scipy, matplotlib, pysiaf, astroquery, requests.
IRSA retrieval uses astroquery's documented SIA2 interface. Network retrieval
must be verified on a machine with internet access. No JWST pipeline run needed.
"""
import argparse
import json
import os
from pathlib import Path

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parent / '.mplconfig'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as path_effects
import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.visualization import ImageNormalize, AsinhStretch, PercentileInterval
from astropy.wcs import WCS
import pysiaf
from pysiaf.utils import rotations

NAMES = ([f'NRS_FULL_MSA{i}' for i in range(1, 5)] + ['NRS_FULL_IFU'] +
         [f'NRS_{s}_SLIT' for s in ['S200A1', 'S200A2', 'S400A1', 'S1600A1', 'S200B1']])


def geometry(path):
    with fits.open(path) as hdus:
        primary = hdus[0].header.copy()
        if primary.get('INSTRUME', '').upper() != 'NIRSPEC':
            raise ValueError('Input must identify INSTRUME=NIRSPEC.')
        # Prefer SCI attitude metadata, falling back to the primary header.
        header = primary.copy()
        if 'SCI' in hdus:
            header.update(hdus['SCI'].header)
    ref_keys = ['V2_REF', 'V3_REF', 'RA_REF', 'DEC_REF', 'ROLL_REF']
    if all(k in header for k in ref_keys):
        attitude = rotations.attitude(*[float(header[k]) for k in ref_keys])
        source = {k: float(header[k]) for k in ref_keys}
    elif all(k in header for k in ['RA_V1', 'DEC_V1', 'PA_V3']):
        attitude = rotations.attitude(0, 0, header['RA_V1'], header['DEC_V1'], header['PA_V3'])
        source = {k: float(header[k]) for k in ['RA_V1', 'DEC_V1', 'PA_V3']}
    else:
        raise ValueError('Missing attitude metadata: need V2_REF/V3_REF/RA_REF/DEC_REF/ROLL_REF '
                         'or RA_V1/DEC_V1/PA_V3. Supply a constituent exposure with these headers.')
    siaf = pysiaf.Siaf('NIRSpec')
    footprints = {}
    for name in NAMES:
        # Stored SIAF ideal vertices define the aperture boundaries.
        v2, v3 = siaf[name].closed_polygon_points('tel', rederive=False)
        ra, dec = rotations.pointing(attitude, v2, v3)
        footprints[name] = SkyCoord(ra, dec, unit='deg')
    all_corners = SkyCoord(np.concatenate([v.ra.deg for v in footprints.values()]),
                          np.concatenate([v.dec.deg for v in footprints.values()]), unit='deg')
    anchor = footprints['NRS_FULL_IFU'][0]
    offsets = all_corners.transform_to(anchor.skyoffset_frame())
    center = SkyCoord(lon=(offsets.lon.min() + offsets.lon.max()) / 2,
                      lat=(offsets.lat.min() + offsets.lat.max()) / 2,
                      frame=anchor.skyoffset_frame()).icrs
    size = max(np.ptp(offsets.lon.arcmin), np.ptp(offsets.lat.arcmin)) * 1.25
    target = None
    if 'TARG_RA' in header and 'TARG_DEC' in header:
        target = SkyCoord(header['TARG_RA'], header['TARG_DEC'], unit='deg')
    return footprints, center, size, target, source, header


def image_hdu(path):
    with fits.open(path) as hdus:
        for hdu in hdus:
            if hdu.data is not None and hdu.data.ndim == 2:
                wcs = WCS(hdu.header)
                if wcs.has_celestial:
                    return np.array(hdu.data, dtype=float), wcs.celestial
    raise ValueError(f'No 2D celestial-WCS image found in {path}')


def wise_image(center, size, cache):
    """Discover W1 intensity FITS images via IRSA SIA; crop a covering tile.

    Do not assume a collection identifier or construct undocumented cutout URLs.
    Save the discovery table for troubleshooting service-dependent responses.
    """
    from astroquery.ipac.irsa import Irsa
    import requests
    cache.mkdir(parents=True, exist_ok=True)
    key = f'w1_{center.ra.deg:.6f}_{center.dec.deg:.6f}_{size:.3f}'
    cutout_path = cache / f'{key}.fits'
    if cutout_path.exists():
        return image_hdu(cutout_path), str(cutout_path)
    print('Querying IRSA SIA for WISE W1 images...', flush=True)
    results = Irsa.query_sia(pos=(center.ra.deg, center.dec.deg, size / 120),
                             band=3.4e-6, data_type='image')
    table_path = cache / f'{key}_sia.ecsv'
    results.write(table_path, format='ascii.ecsv', overwrite=True)
    candidates = []
    for row in results:
        def value(k):
            return str(row[k]) if k in results.colnames else ''
        description = ' '.join(value(k) for k in ['obs_collection', 'instrument_name',
                                                   'obs_title', 'obs_id', 'access_url']).lower()
        url, fmt = value('access_url'), value('access_format').lower()
        if 'wise' not in description or 'fits' not in fmt or not url.startswith(('https://', 'http://')):
            continue
        if any(word in description for word in ['uncert', 'coverage', '-unc-', '-cov-', '-msk-']):
            continue
        candidates.append((0 if 'allwise' in description else 1, url))
    if not candidates:
        raise RuntimeError(f'IRSA returned no directly accessible WISE W1 FITS image. '
                           f'Inspect {table_path}, or supply --image a local W1 FITS image.')
    failures = []
    for _, url in sorted(set(candidates)):
        tile_path = cache / f'{key}_tile.fits'
        try:
            print(f'Downloading {url}', flush=True)
            with requests.get(url, stream=True, timeout=(30, 180)) as response:
                response.raise_for_status()
                with tile_path.open('wb') as output:
                    for chunk in response.iter_content(1024 * 1024):
                        output.write(chunk)
            data, wcs = image_hdu(tile_path)
            cut = Cutout2D(data, center, size * u.arcmin, wcs=wcs, mode='strict')
            fits.PrimaryHDU(cut.data, cut.wcs.to_header()).writeto(cutout_path, overwrite=True)
            tile_path.unlink()
            (cache / f'{key}.json').write_text(json.dumps({'access_url': url}, indent=2))
            return (cut.data, cut.wcs), url
        except Exception as error:
            failures.append(f'{url}: {error}')
        finally:
            tile_path.unlink(missing_ok=True)
    raise RuntimeError('No returned WISE tile covered the chart. Use --image with a mosaic.\n' +
                       '\n'.join(failures))



def orient_msa(data, wcs, center, size, footprints, height_arcmin=None):
    """Align average MSA edges with plot axes, choosing IFU-left orientation.

    SIAF quadrants have slight geometric distortion; align their mean edge angle.
    """
    from scipy.ndimage import map_coordinates
    from astropy.wcs.utils import proj_plane_pixel_scales
    vertices = footprints['NRS_FULL_IFU'][:-1].transform_to(center.skyoffset_frame())
    east, north = np.mean(vertices.lon.deg), np.mean(vertices.lat.deg)
    length = np.hypot(east, north)
    if length == 0:
        raise ValueError('Cannot set orientation: IFU coincides with chart center.')
    # Pixel +x points away from IFU; +y preserves the usual sky handedness.
    preferred_x = np.array([-east / length, -north / length])
    edge_angles = []
    for i in range(1, 5):
        corners = footprints[f'NRS_FULL_MSA{i}'].transform_to(center.skyoffset_frame())
        de, dn = np.diff(corners.lon.deg), np.diff(corners.lat.deg)
        edge_angles.extend(np.arctan2(dn, de))
    base_angle = np.angle(np.mean(np.exp(4j * np.array(edge_angles)))) / 4
    candidates = [np.array([np.cos(base_angle + i * np.pi / 2),
                            np.sin(base_angle + i * np.pi / 2)]) for i in range(4)]
    x_east, x_north = max(candidates, key=lambda vector: np.dot(vector, preferred_x))
    scale = float(np.min(np.abs(proj_plane_pixel_scales(wcs))))
    count = int(np.ceil(size / 60 / scale))
    rows = count if height_arcmin is None else int(np.ceil(height_arcmin / 60 / scale))
    output_wcs = WCS(naxis=2)
    output_wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    output_wcs.wcs.crval = [center.ra.deg, center.dec.deg]
    output_wcs.wcs.crpix = [(count + 1) / 2, (rows + 1) / 2]
    output_wcs.wcs.cd = scale * np.array([[x_east, x_north], [x_north, -x_east]])
    yy, xx = np.indices((rows, count), dtype=float)
    sky = output_wcs.pixel_to_world(xx, yy)
    sx, sy = wcs.world_to_pixel(sky)
    output = map_coordinates(data, [sy, sx], order=1, mode='constant', cval=np.nan)
    return output, output_wcs


def write_regions(path, footprints):
    lines = ['# Region file format: DS9 version 4.1', 'global width=2', 'fk5']
    for name, coordinates in footprints.items():
        color = 'cyan' if 'MSA' in name else ('magenta' if 'IFU' in name else 'yellow')
        vertices = ','.join(f'{ra:.9f},{dec:.9f}' for ra, dec in
                            zip(coordinates.ra.deg[:-1], coordinates.dec.deg[:-1]))
        lines.append(f'polygon({vertices}) # color={color} text={{{name}}}')
    path.write_text('\n'.join(lines) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--output', type=Path, help='PNG path; default beside input')
    parser.add_argument('--image', type=Path, help='Use any user-supplied 2D celestial-WCS FITS background instead of WISE')
    parser.add_argument('--cache', type=Path, default=Path(__file__).resolve().parent / 'finder_cache')
    parser.add_argument('--size', type=float, help='Square chart width in arcminutes')
    parser.add_argument('--FSzoom', action='store_true', help='Zoom to IFU and adjacent fixed slits, about 100 x 60 arcsec; --size overrides width')
    parser.add_argument('--regions', action='store_true', help='Also write DS9 aperture polygons')
    args = parser.parse_args()
    footprints, center, automatic_size, target, attitude, header = geometry(args.input)
    if args.FSzoom:
        zoom_names = ['NRS_FULL_IFU', 'NRS_S200A1_SLIT', 'NRS_S200A2_SLIT',
                      'NRS_S400A1_SLIT', 'NRS_S1600A1_SLIT']
        corners = SkyCoord(np.concatenate([footprints[name].ra.deg[:-1] for name in zoom_names]),
                           np.concatenate([footprints[name].dec.deg[:-1] for name in zoom_names]), unit='deg')
        frame = footprints['NRS_S200A1_SLIT'][0].skyoffset_frame()
        offsets = corners.transform_to(frame)
        center = SkyCoord(lon=(offsets.lon.min() + offsets.lon.max()) / 2,
                          lat=(offsets.lat.min() + offsets.lat.max()) / 2, frame=frame).icrs
        automatic_size = 100 / 60
    size = args.size if args.size is not None else automatic_size
    height = size * .6 if args.FSzoom else size
    if not np.isfinite(size) or size <= 0:
        parser.error('--size must be positive and finite')
    output = args.output or args.input.with_name(args.input.stem + '_finder.png')
    output.parent.mkdir(parents=True, exist_ok=True)
    # Regions can be useful even when image retrieval fails.
    if args.regions:
        write_regions(output.with_suffix('.reg'), footprints)
    if args.image:
        data, wcs = image_hdu(args.image)
        cut = Cutout2D(data, center, (height, size) * u.arcmin, wcs=wcs, mode='partial', fill_value=np.nan)
        (data, wcs), provenance = (cut.data, cut.wcs), str(args.image)
    else:
        (data, wcs), provenance = wise_image(center, size, args.cache)
    data, wcs = orient_msa(data, wcs, center, size, footprints,
                           height_arcmin=height if args.FSzoom else None)
    finite = data[np.isfinite(data)]
    if not len(finite):
        raise ValueError('Image has no finite pixels.')
    vmin, vmax = PercentileInterval(99).get_limits(finite)
    fig = plt.figure(figsize=(12, 7.2) if args.FSzoom else (9, 9))
    ax = fig.add_subplot(111, projection=wcs)
    ax.imshow(data, origin='lower', cmap='gray', norm=ImageNormalize(vmin=vmin, vmax=vmax,
                                                                  stretch=AsinhStretch()))
    # Separate the crowded central apertures into a label column. Preserve
    # their vertical order so their leader lines do not cross one another.
    central_names = ['NRS_FULL_IFU', 'NRS_S200A1_SLIT', 'NRS_S200A2_SLIT',
                     'NRS_S400A1_SLIT', 'NRS_S1600A1_SLIT']
    label_order = sorted(central_names, key=lambda name:
                         np.mean(wcs.world_to_pixel(footprints[name])[1][:-1]), reverse=True)
    label_positions = {name: (0.065, y) for name, y in
                       zip(label_order, np.linspace(0.61, 0.39, len(label_order)))}
    for name, coordinates in footprints.items():
        color = 'cyan' if 'MSA' in name else ('magenta' if 'IFU' in name else 'yellow')
        x, y = wcs.world_to_pixel(coordinates)
        # Narrow fixed slits are only a few display pixels wide in FSzoom;
        # thick outlines hide the background inside their true boundaries.
        ax.plot(x, y, color=color, linewidth=0.4 if args.FSzoom else 1.3)
        label = name.replace('NRS_FULL_', '').replace('NRS_', '').replace('_SLIT', '')
        if args.FSzoom:
            if name in central_names:
                ax.annotate(label, (np.mean(x[:-1]), np.max(y)), xytext=(0, 7),
                            textcoords='offset points', ha='center', color=color, fontsize=8,
                            bbox=dict(facecolor='black', edgecolor='none', alpha=.5, pad=1))
        elif name in label_positions:
            ax.annotate(label, (np.mean(x[:-1]), np.mean(y[:-1])),
                        xytext=label_positions[name], textcoords='axes fraction',
                        ha='left', va='center', color=color, fontsize=8,
                        bbox=dict(facecolor='black', edgecolor='none', alpha=0.65, pad=2),
                        arrowprops=dict(arrowstyle='-', color=color, linewidth=0.7,
                                        shrinkA=4, shrinkB=3))
        else:
            ax.annotate(label, (np.mean(x[:-1]), np.mean(y[:-1])), xytext=(4, 5),
                        textcoords='offset points', color=color, fontsize=8)
    ax.set_xlim(-0.5, data.shape[1] - 0.5)
    ax.set_ylim(-0.5, data.shape[0] - 0.5)
    ax.coords[0].set_axislabel('Right ascension (ICRS)')
    ax.coords[1].set_axislabel('Declination (ICRS)')
    ax.coords.grid(color='white', alpha=0.25)
    # Conventional two-arm N/E compass, derived from the displayed celestial WCS.
    # Axes coordinates keep the marker a fixed size regardless of image dimensions.
    cx, cy = wcs.world_to_pixel(center)
    origin = np.array([0.12, 0.88])
    ax.plot(origin[0], origin[1], 'o', color='black', markersize=3,
            transform=ax.transAxes, zorder=10)
    for label, angle in [('N', 0 * u.deg), ('E', 90 * u.deg)]:
        sky_tip = center.directional_offset_by(angle, 10 * u.arcsec)
        tx, ty = wcs.world_to_pixel(sky_tip)
        direction = np.array([(tx - cx) / data.shape[1], (ty - cy) / data.shape[0]])
        direction /= np.linalg.norm(direction)
        arm_length = 0.08 if label == 'N' else 0.05
        tip = origin + arm_length * direction
        ax.plot([origin[0], tip[0]], [origin[1], tip[1]], color='black', linewidth=2,
                transform=ax.transAxes, zorder=10,
                path_effects=[path_effects.Stroke(linewidth=4, foreground='white'),
                              path_effects.Normal()])
        text_position = origin + (arm_length + 0.018) * direction
        ax.text(*text_position, label, color='black', fontsize=11, fontweight='bold',
                ha='center', va='center', transform=ax.transAxes, zorder=10,
                bbox=dict(facecolor='white', edgecolor='none', alpha=0.8, pad=1))
    fid_aperture = pysiaf.Siaf('NIRSpec')['NRS_S200A1_SLIT']
    if 'V2_REF' in attitude:
        matrix = rotations.attitude(*[attitude[k] for k in
                                      ['V2_REF', 'V3_REF', 'RA_REF', 'DEC_REF', 'ROLL_REF']])
    else:
        matrix = rotations.attitude(0, 0, attitude['RA_V1'], attitude['DEC_V1'], attitude['PA_V3'])
    fid_ra, fid_dec = rotations.pointing(matrix, fid_aperture.V2Ref, fid_aperture.V3Ref)
    fiducial = SkyCoord(fid_ra, fid_dec, unit='deg')
    fid_text = fiducial.to_string('hmsdms', sep=':', precision=2)
    background_label = args.image.name if args.image else 'WISE W1 (3.4 μm)'
    ax.set_title(f'NIRSpec overview • {background_label}\n{args.input.name}\n'
                 f'S200A1 center (ICRS): {fid_text}', fontsize=10, pad=12)
    fig.text(0.5, 0.025, 'MSA: cyan   IFU: magenta   Fixed slits: yellow',
             ha='center', fontsize=9)
    fig.text(0.5, 0.005, f'Field of view: {size:.2f} × {height:.2f} arcmin',
             ha='center', fontsize=9)
    fig.savefig(output, dpi=180, bbox_inches='tight')
    plt.close(fig)
    output.with_suffix('.json').write_text(json.dumps({
        'input': str(args.input.resolve()), 'image': provenance, 'attitude': attitude,
        'pysiaf_version': pysiaf.__version__, 'siaf_prd': str(pysiaf.JWST_PRD_VERSION),
        'orientation': 'Mean MSA edges aligned with plot axes; IFU left',
        'fiducial': {'aperture': 'S200A1', 'ra_deg': float(fid_ra), 'dec_deg': float(fid_dec)},
        'size_arcmin': size, 'height_arcmin': height, 'FSzoom': args.FSzoom, 'exposure_type': header.get('EXP_TYPE'),
        'note': 'Single header attitude; combined-product constituent pointings are not reconstructed.'
    }, indent=2))
    print(f'Saved {output}')


if __name__ == '__main__':
    main()

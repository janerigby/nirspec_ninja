"""Contemporaneous NIRSpec fixed-slit sky backgrounds for calibrated IFU data.

No dependency on the LEGGOS driver. Public functions accept paths/arrays;
science decisions and output provenance are explicit. Surface brightness is
F_nu in MJy/sr. Background statistical errors are approximate and are NOT
added to science ERR: one background is shared by many science pixels.
"""
from dataclasses import dataclass
from pathlib import Path
import csv
import json
import re
import warnings

import numpy as np
from astropy import units as u
from astropy.io import fits
from jwst import datamodels
from jwst.lib.wcs_utils import get_wavelengths

USABLE_SLITS = ('S200A1', 'S200A2', 'S400A1')
TRACE_THRESHOLD = 2.5


def exposure_key(filename):
    """Match detector pairs of the exact JWST exposure, never by target/date."""
    name = Path(filename).name
    match = re.fullmatch(r'(jw\d+_\d+_\d+)_nrs[12]_(?:cal|s2d)\.fits', name, re.I)
    if not match:
        raise ValueError(f'Cannot identify exact JWST exposure from {name}')
    return match.group(1).lower()


def exposure_interval(filename):
    header = fits.getheader(filename, 0)
    start, end = float(header['EXPSTART']), float(header['EXPEND'])
    if not np.isfinite(start + end) or end <= start:
        raise ValueError(f'Invalid exposure time interval: {filename}')
    return start, end


def contemporaneous_s2d(ifu_file, s2d_files):
    """Require matching exposure IDs, wavelength configuration and time overlap."""
    key = exposure_key(ifu_file)
    sci_header = fits.getheader(ifu_file, 0)
    start, end = exposure_interval(ifu_file)
    matches = []
    for filename in s2d_files:
        if exposure_key(filename) != key:
            continue
        header = fits.getheader(filename, 0)
        if any(header[k] != sci_header[k] for k in ('GRATING', 'FILTER')):
            raise ValueError(f'Configuration mismatch for contemporaneous FS file {filename}')
        bg_start, bg_end = exposure_interval(filename)
        if min(end, bg_end) <= max(start, bg_start):
            raise ValueError(f'FS and IFU exposure intervals do not overlap: {filename}')
        matches.append(str(Path(filename).resolve()))
    if not matches:
        raise FileNotFoundError(f'No contemporaneous FS S2D for {ifu_file}; regenerate FS Spec2')
    return matches


@dataclass
class SlitPixels:
    name: str
    brightness: np.ndarray
    wavelength: np.ndarray
    error: np.ndarray
    dq: np.ndarray | None
    calibrated: bool = True


def read_slit_pixels(filename):
    """Read S2D pixels; DQ is optional because resampling defaults not to propagate it.

    Mask nonpositive/nonfinite WHT when supplied. Keep absent DQ as None, not
    invented flags; detection and row selection use SCI/ERR validity as well.
    """
    with datamodels.open(filename) as model:
        slits = model.slits if hasattr(model, 'slits') else [model]
        for slit in slits:
            name = str(slit.name).upper()
            unit = slit.meta.bunit_data
            try:
                if not unit:
                    raise ValueError('Missing brightness unit')
                scale = u.Unit(unit).to(u.MJy / u.sr)
                calibrated = True
            except (ValueError, u.UnitConversionError):
                if name in USABLE_SLITS:
                    raise ValueError(f'Expected MJy/sr brightness in {filename}, {name}')
                # Ignored, potentially uncalibrated slits can still be screened
                # in their native units: detection significance is invariant
                # to common SCI/ERR scaling. Never use these values for sky.
                scale, calibrated = 1.0, False
            if name in USABLE_SLITS and slit.source_type != 'EXTENDED':
                raise ValueError(f'Expected EXTENDED calibration in {filename}, {name}')
            wavelength = get_wavelengths(slit, 'NRS_FIXEDSLIT')
            if wavelength is None:
                raise ValueError(f'No wavelength array for {filename}, {slit.name}')
            brightness = np.array(slit.data, dtype=float) * scale
            if slit.err is None:
                raise ValueError(f'Missing required ERR array in {filename}, {name}')
            error = np.array(slit.err, dtype=float) * scale
            wavelength = np.array(wavelength, dtype=float)
            dq = None if slit.dq is None else np.array(slit.dq, dtype=np.uint32)
            if brightness.ndim != 2 or error.shape != brightness.shape or wavelength.shape != brightness.shape:
                raise ValueError(f'Incompatible S2D SCI/ERR/wavelength shapes in {filename}, {name}')
            if dq is not None and dq.shape != brightness.shape:
                raise ValueError(f'Incompatible S2D DQ shape in {filename}, {name}')
            weight = getattr(slit, 'wht', None)
            if weight is not None:
                weight = np.asarray(weight)
                if weight.shape != brightness.shape:
                    raise ValueError(f'Incompatible S2D WHT shape in {filename}, {name}')
                unusable = ~np.isfinite(weight) | (weight <= 0)
                brightness[unusable] = np.nan
                error[unusable] = np.nan
            yield SlitPixels(name, brightness, wavelength, error, dq, calibrated)



def _finite_scalar(value):
    if value is None:
        return None
    number = float(value)
    return number if np.isfinite(number) else None


def save_trace_plot(slit, result, filename):
    """Save a detector diagnostic with a file canvas; no GUI/backend changes."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from astropy.visualization import ImageNormalize, ZScaleInterval
    figure = Figure(figsize=(13, 5))
    FigureCanvasAgg(figure)
    image, profile = figure.subplots(1, 2, sharey=True, width_ratios=[3, 1.3])
    from matplotlib import colormaps
    cmap = colormaps['gray_r'].copy()
    cmap.set_bad('cyan')
    finite = np.isfinite(slit.brightness)
    occupied_rows = np.flatnonzero(np.any(finite, axis=1))
    occupied_columns = np.flatnonzero(np.any(finite, axis=0))
    if not occupied_rows.size or not occupied_columns.size:
        raise ValueError('Cannot plot a slit with no finite image pixels')
    first_row, last_row = occupied_rows[[0, -1]]
    first_column, last_column = occupied_columns[[0, -1]]
    normalization = ImageNormalize(slit.brightness, interval=ZScaleInterval(contrast=.5))
    colors = cmap(normalization(np.ma.masked_invalid(slit.brightness)))
    # Do not render entirely empty rows/columns, including interior gaps.
    empty = ~np.any(finite, axis=1)[:, None] | ~np.any(finite, axis=0)[None, :]
    colors[empty, 3] = 0
    image.imshow(colors[first_row:last_row + 1, first_column:last_column + 1],
                 origin='lower', aspect='auto',
                 extent=(first_column - .5, last_column + .5, first_row - .5, last_row + .5))
    image.axhline(result['center'], color='red', linestyle='--')
    half_width = result['fwhm'] / 2
    for boundary in (result['center'] - half_width, result['center'] + half_width):
        image.axhline(boundary, color='red', linewidth=.8)
    excluded = result.get('excluded_rows', np.zeros(slit.brightness.shape[0], bool))
    for index, row in enumerate(np.flatnonzero(excluded)):
        image.axhspan(row - .5, row + .5, color='orange', alpha=.18)
        profile.axhspan(row - .5, row + .5, color='orange', alpha=.18,
                        label='Excluded rows' if index == 0 else None)
    image.set_ylim(first_row - .5, last_row + .5)
    image.set_xlim(first_column - .5, last_column + .5)
    label = 'TRACE DETECTED' if result['significance'] >= result.get('threshold', TRACE_THRESHOLD) else 'NO TRACE DETECTED'
    image.set(title=f"{slit.name}: {label} ({result['significance']:.2f} sigma)",
              xlabel='Dispersion pixel', ylabel='Spatial row')
    rows = np.arange(len(result['profile']))
    profile.errorbar(result['profile'], rows, xerr=result['profile_error'], fmt='.', label='Row median')
    profile.plot(np.where(excluded, np.nan, result['background']), rows, '--', label='Fitted background')
    profile.plot(np.where(excluded, np.nan, result['model']), rows, label='Background + trace')
    profile.axhline(result['center'], color='red', linestyle='--')
    for index, boundary in enumerate((result['center'] - half_width, result['center'] + half_width)):
        profile.axhline(boundary, color='red', linewidth=.8,
                        label=f"FWHM = {result['fwhm']:.1f} rows" if index == 0 else None)
    profile.set(xlabel='Surface brightness [MJy/sr]', title='Spatial profile')
    profile.legend(fontsize=8, loc='upper center', bbox_to_anchor=(.5, -.18))
    figure.tight_layout()
    figure.savefig(filename, dpi=150)
    figure.clear()


def scan_fixed_slit_traces(s2d_files, output_dir, target, threshold=TRACE_THRESHOLD,
                          make_plots=True, screening_scope='exposure'):
    """Screen S2Ds with source_detection; configuration scope applies L3 decisions.

    Reports include every slit. Only the three approved slit names can be
    accepted for sky measurement. Failed fits are rejected, never treated as
    evidence of an empty slit. No emission-line screening is added.
    """
    from nirspec_ninja.fixed_slit import source_detection as detector
    if screening_scope not in ('exposure', 'configuration'):
        raise ValueError('screening_scope must be exposure or configuration')
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    records = []
    for filename in sorted(s2d_files):
        filename = str(Path(filename).resolve())
        header = fits.getheader(filename, 0)
        start, end = exposure_interval(filename)
        key = exposure_key(filename) if screening_scope == 'exposure' else 'combined L3'
        for slit in read_slit_pixels(filename):
            if slit.name.startswith('S200B'):
                continue
            good, excluded = spatial_pixel_mask(slit, clip_nan=5)
            fit_data = np.where(good & ~excluded[:, None], slit.brightness, np.nan)
            result = detector.detect_source_s2d(
                fit_data, error_2d=slit.error, dq_2d=slit.dq,
                allowed_source_center_rows=~excluded,
                slit_length_pixels=detector.NIRSPEC_FIXED_SLIT_LENGTH_ARCSEC[slit.name]
                / detector.NIRSPEC_SPATIAL_PIXEL_SCALE_ARCSEC,
            )
            sigma = _finite_scalar(result.get('significance'))
            ok = result['status'] == 'ok' and sigma is not None
            trace = ok and sigma >= threshold
            accepted = slit.name in USABLE_SLITS and ok and not trace
            reason = ('ignored slit' if slit.name not in USABLE_SLITS else
                      'failed detection' if not ok else
                      'continuum trace' if trace else 'accepted')
            plot = None
            if make_plots and ok and slit.calibrated:
                result['excluded_rows'] = excluded
                result['threshold'] = threshold
                plot_path = output_dir / (Path(filename).stem + '_' + slit.name + '_trace.png')
                save_trace_plot(slit, result, str(plot_path))
                plot = str(plot_path.resolve())
            records.append(dict(
                target=target, exposure=key, file=filename, detector=header['DETECTOR'],
                screening_scope=screening_scope,
                grating=header['GRATING'], filter=header['FILTER'],
                start_mjd=start, end_mjd=end, slit=slit.name,
                status=result['status'], dq_available=slit.dq is not None, corrected_sigma=sigma,
                raw_sigma=_finite_scalar(result.get('raw_significance')),
                center=_finite_scalar(result.get('center')),
                fwhm=_finite_scalar(result.get('fwhm')),
                excluded_rows=np.flatnonzero(excluded).tolist(),
                block_significances=[_finite_scalar(value) for value in result.get('block_significances', [])],
                blocks_above_3sigma=result.get('blocks_above_3sigma'),
                accepted=accepted, calibrated=slit.calibrated, decision=reason, reason=result.get('reason', ''), plot=plot,
            ))
    (output_dir / 'trace_report.json').write_text(json.dumps(records, indent=2, allow_nan=False) + '\n')
    columns = ['target', 'screening_scope', 'exposure', 'detector', 'grating', 'filter', 'slit',
               'corrected_sigma', 'raw_sigma', 'status', 'dq_available', 'accepted', 'decision',
               'start_mjd', 'end_mjd', 'file', 'reason', 'plot']
    with (output_dir / 'trace_report.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        writer.writerows(records)
    lines = [f'Target: {target}', f'Screening scope: {screening_scope}',
             f'Continuum rejection: corrected significance >= {threshold:g} sigma',
             'Failures excluded; no emission-line screening.',
             'EXPOSURE | DETECTOR | GRATING/FILTER | SLIT | SIGMA | DECISION']
    for record in records:
        sigma = record['corrected_sigma']
        lines.append(f"{record['exposure']} | {record['detector']} | "
                     f"{record['grating']}/{record['filter']} | {record['slit']} | "
                     f"{sigma:.2f}" if sigma is not None else
                     f"{record['exposure']} | {record['detector']} | "
                     f"{record['grating']}/{record['filter']} | {record['slit']} | NA")
        lines[-1] += ' | ' + record['decision']
    (output_dir / 'trace_report.txt').write_text('\n'.join(lines) + '\n')
    return records


def bin_edges(centers):
    centers = np.asarray(centers, dtype=float)
    if centers.ndim != 1 or len(centers) < 2 or not np.all(np.isfinite(centers)):
        raise ValueError('At least two finite spectral centers are required')
    if not np.all(np.diff(centers) > 0):
        raise ValueError('Wavelength centers must be strictly increasing')
    edges = np.r_[centers[0] - (centers[1] - centers[0]) / 2,
                  (centers[:-1] + centers[1:]) / 2,
                  centers[-1] + (centers[-1] - centers[-2]) / 2]
    if edges[0] <= 0:
        raise ValueError('Wavelength bin edges must be positive')
    return edges


def resample_surface_brightness(wavelength, brightness, output_wavelength):
    """Bin-overlap resampling of F_nu, integrating in frequency, not wavelength.

    Output bins must be completely covered by valid native bins. Missing
    pixels/gaps are not bridged or extrapolated. No bin-width multiplication
    is left in the resulting brightness density. Input values represent
    bin-averaged F_nu; wavelength edges are inferred from center midpoints.
    """
    wave = np.asarray(wavelength, dtype=float)
    values = np.asarray(brightness, dtype=float)
    if wave.shape != values.shape or wave.ndim != 1:
        raise ValueError('Native wavelengths and brightness must be matching 1D arrays')
    output_edges = 299792458.0e6 / bin_edges(output_wavelength)
    out_lo, out_hi = output_edges[1:], output_edges[:-1]
    total_integral = np.zeros(len(out_lo))
    total_coverage = np.zeros(len(out_lo))
    finite = np.isfinite(wave) & (wave > 0)
    starts = np.flatnonzero(finite & ~np.r_[False, finite[:-1]])
    ends = np.flatnonzero(finite & ~np.r_[finite[1:], False]) + 1
    for start, end in zip(starts, ends):
        if end - start < 2:
            continue
        native_wave = wave[start:end]
        native_values = values[start:end]
        if np.all(np.diff(native_wave) < 0):
            native_wave, native_values = native_wave[::-1], native_values[::-1]
        edges = (299792458.0e6 / bin_edges(native_wave))[::-1]
        native_values = native_values[::-1]
        good = np.isfinite(native_values)
        widths = np.diff(edges)
        integral = np.r_[0.0, np.cumsum(np.where(good, native_values, 0.0) * widths)]
        coverage = np.r_[0.0, np.cumsum(good * widths)]
        total_integral += np.interp(out_hi, edges, integral) - np.interp(out_lo, edges, integral)
        total_coverage += np.interp(out_hi, edges, coverage) - np.interp(out_lo, edges, coverage)
    output_widths = out_hi - out_lo
    complete = np.isclose(total_coverage, output_widths, rtol=1e-8, atol=0)
    result = np.full(len(out_lo), np.nan)
    result[complete] = total_integral[complete] / output_widths[complete]
    return result


def spatial_pixel_mask(slit, clip_nan=5):
    """Shared notebook row exclusions for trace fitting and sky estimation."""
    good = (np.isfinite(slit.brightness) & np.isfinite(slit.wavelength)
            & np.isfinite(slit.error) & (slit.error > 0))
    if slit.dq is not None:
        good &= (slit.dq & datamodels.dqflags.pixel['DO_NOT_USE']) == 0
    mostly_invalid = 2 * np.sum(~good, axis=1) >= good.shape[1]
    anchors = np.unique(np.r_[np.flatnonzero(mostly_invalid), 0, good.shape[0] - 1])
    excluded = np.zeros(good.shape[0], dtype=bool)
    for row in anchors:
        excluded[max(0, row - clip_nan):min(len(excluded), row + clip_nan + 1)] = True
    return good, excluded


def usable_spatial_rows(slit, clip_nan=5):
    """Notebook edge/mostly-invalid-row exclusion plus SCI/ERR/DQ masking."""
    good, excluded = spatial_pixel_mask(slit, clip_nan)
    values = np.where(good, slit.brightness, np.nan)
    return [(slit.wavelength[row], values[row]) for row in np.flatnonzero(~excluded)]


def _nanmedian(values, axis=0):
    with warnings.catch_warnings():
        warnings.filterwarnings('ignore', message='All-NaN slice encountered', category=RuntimeWarning)
        return np.nanmedian(values, axis=axis)


def estimate_background(s2d_files, trace_records, output_wavelength, clip_nan=5):
    """Median per slit, then equal-weight median across accepted slit spectra.

    Uncertainty: 1.253*1.4826*MAD/sqrt(N_rows) within each slit, then the
    largest contributing slit uncertainty. This is an approximate measure;
    S2D row correlations and systematic calibration errors are not included.
    """
    screening_enabled = trace_records is not None
    lookup = {(r['file'], r['slit']): r for r in (trace_records or [])
              if r.get('screening_scope', 'exposure') == 'exposure'}
    configuration_lookup = {}
    for record in (trace_records or []):
        if record.get('screening_scope') == 'configuration':
            key = (record['grating'].upper(), record['filter'].upper(), record['slit'])
            configuration_lookup.setdefault(key, []).append(record)
    rows_by_slit = {name: [] for name in USABLE_SLITS}
    used = []
    for filename in s2d_files:
        filename = str(Path(filename).resolve())
        header = fits.getheader(filename, 0)
        for slit in read_slit_pixels(filename):
            if slit.name not in USABLE_SLITS:
                continue
            key = (header['GRATING'].upper(), header['FILTER'].upper(), slit.name)
            decisions = configuration_lookup.get(key)
            if not screening_enabled:
                decisions = [dict(file=None, accepted=True)]
            if decisions is None:
                record = lookup.get((filename, slit.name))
                decisions = [] if record is None else [record]
            if not decisions:
                raise ValueError(f'Missing trace decision for {filename}, {slit.name}')
            if not all(record['accepted'] for record in decisions):
                continue
            rows = [resample_surface_brightness(wave, values, output_wavelength)
                    for wave, values in usable_spatial_rows(slit, clip_nan)]
            rows = [row for row in rows if np.any(np.isfinite(row))]
            if rows:
                rows_by_slit[slit.name].extend(rows)
                used.append(dict(file=filename, slit=slit.name, rows=len(rows),
                                 screening_enabled=screening_enabled,
                                 screening_files=[r['file'] for r in decisions if r['file'] is not None]))
    length = len(output_wavelength)
    slit_medians, slit_errors, counts = [], [], []
    for name in USABLE_SLITS:
        rows = rows_by_slit[name]
        if not rows:
            continue
        rows = np.asarray(rows)
        count = np.sum(np.isfinite(rows), axis=0)
        median = _nanmedian(rows)
        mad = _nanmedian(np.abs(rows - median))
        error = np.full(length, np.nan)
        # One sample can supply a brightness, but cannot estimate scatter.
        valid = count >= 2
        error[valid] = 1.253 * 1.4826 * mad[valid] / np.sqrt(count[valid])
        slit_medians.append(median)
        slit_errors.append(error)
        counts.append(count)
    if not slit_medians:
        return dict(wavelength=np.asarray(output_wavelength), brightness=np.full(length, np.nan),
                    uncertainty=np.full(length, np.nan), n_slits=np.zeros(length, int),
                    n_rows=np.zeros(length, int), used=used)
    medians = np.asarray(slit_medians)
    errors = np.asarray(slit_errors)
    n_slits = np.sum(np.isfinite(medians), axis=0)
    uncertainty = np.full(length, np.nan)
    for index in range(length):
        contributors = np.isfinite(medians[:, index])
        if np.any(contributors) and np.all(np.isfinite(errors[contributors, index])):
            uncertainty[index] = np.max(errors[contributors, index])
    return dict(wavelength=np.asarray(output_wavelength), brightness=_nanmedian(medians),
                uncertainty=uncertainty, n_slits=n_slits,
                n_rows=np.sum(counts, axis=0), used=used)


def cube_wavelengths(filename):
    """Read the saved cube WCS, including nonlinear PRISM wavelength axes."""
    with datamodels.open(filename) as cube:
        nz, ny, nx = cube.data.shape
        wcs = cube.meta.wcs
        if wcs is None:
            raise ValueError(f'Cube has no saved WCS: {filename}')
        wave = np.asarray(wcs(np.full(nz, (nx - 1) / 2),
                              np.full(nz, (ny - 1) / 2), np.arange(nz))[2], dtype=float)
        unit = wcs.output_frame.unit[2]
        wave *= u.Unit(unit).to(u.um)
    bin_edges(wave)  # Validate increasing finite centers without changing them.
    return wave


def save_background_spectrum(background, stem):
    """Write the exact cube-grid CSV and master-background FITS, including SB_ERROR."""
    from jwst.master_background.create_master_bkg import create_background
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    with stem.with_suffix('.csv').open('w', newline='') as handle:
        handle.write('# wavelength: micron; brightness and uncertainty: MJy/sr\n')
        handle.write('# uncertainty: max slit 1.253*1.4826*MAD/sqrt(N_rows); correlations/systematics excluded\n')
        writer = csv.writer(handle)
        writer.writerow(['wave', 'surf_bright', 'sb_error', 'n_slits', 'n_rows'])
        writer.writerows(zip(background['wavelength'], background['brightness'],
                             background['uncertainty'], background['n_slits'], background['n_rows']))
    model = create_background(background['wavelength'], background['brightness'])
    if model is None:
        raise RuntimeError('Cannot create JWST master background model')
    try:
        model.spec[0].spec_table['SB_ERROR'] = background['uncertainty']
        model.spec[0].spec_table['DQ'] = np.where(
            np.isfinite(background['brightness']), 0, datamodels.dqflags.pixel['DO_NOT_USE'])
        model.save(stem.with_suffix('.fits'))
    finally:
        model.close()
    return str(stem.with_suffix('.fits').resolve())


def subtract_ifu_background(ifu_file, background_file, output_file):
    """Run MasterBackgroundStep, then explicitly flag uncovered IFU pixels.

    The pipeline's IFU path otherwise silently applies zero beyond background
    coverage. Flag pixels for interpolation segments touching a missing bin,
    as well as outside the provided range. Preserve science ERR unchanged.
    """
    from jwst.master_background import MasterBackgroundStep
    from jwst.assign_wcs import nirspec
    from gwcs.wcstools import grid_from_bounding_box
    with datamodels.open(background_file) as background:
        table = background.spec[0].spec_table
        wave = np.array(table['WAVELENGTH'], dtype=float)
        valid = np.isfinite(table['SURF_BRIGHT'])
    step = MasterBackgroundStep(user_background=str(background_file),
                                save_results=False, save_background=False)
    result = step.run(str(ifu_file))
    try:
        if result.meta.cal_step.master_background != 'COMPLETE':
            raise RuntimeError(f'Master background was not applied to {ifu_file}')
        flagged = 0
        for wcs in nirspec.nrs_ifu_wcs(result):
            x, y = grid_from_bounding_box(wcs.bounding_box)
            wl = np.asarray(wcs(x, y)[2])
            science = np.isfinite(wl)
            coverage = np.interp(wl, wave, valid.astype(float), left=0.0, right=0.0)
            bad = science & (coverage < 1.0 - 1e-12)
            xx, yy = x[bad].astype(int), y[bad].astype(int)
            result.data[yy, xx] = np.nan
            result.dq[yy, xx] |= datamodels.dqflags.pixel['DO_NOT_USE']
            flagged += int(bad.sum())
        result.save(str(output_file))
        return flagged
    finally:
        result.close()


def write_background_report(records, output_dir):
    """Update a short target report and machine-readable per-exposure provenance."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    path = output_dir / 'background_report.json'
    previous = json.loads(path.read_text()) if path.exists() else []
    merged = {(r['cube'], r['ifu_file']): r for r in previous}
    merged.update({(r['cube'], r['ifu_file']): r for r in records})
    merged = sorted(merged.values(), key=lambda r: (r['cube'], r['exposure'], r['ifu_file']))
    path.write_text(json.dumps(merged, indent=2, allow_nan=False) + '\n')
    lines = ['Contemporaneous FS backgrounds; each IFU exposure uses only its own FS exposure.',
             'Slits: S200A1/S200A2/S400A1; source screening status is recorded per background.',
             'Uncertainty is approximate; retained separately, not propagated into cube ERR.',
             'TARGET | CUBE | EXPOSURE | MJD START-END | USED SLITS | MISSING BINS | SCREENING | STATUS']
    for record in merged:
        names = ','.join(sorted({r['slit'] for r in record['used']})) or 'NONE'
        lines.append(f"{record['target']} | {Path(record['cube']).name} | {record['exposure']} | "
                     f"{record['start_mjd']:.6f}-{record['end_mjd']:.6f} | {names} | "
                     f"{record['missing_bins']}/{record['total_bins']} | "
                     f"{'enabled' if record.get('trace_screening_enabled', True) else 'DISABLED'} | {record['status']}")
    (output_dir / 'background_report.txt').write_text('\n'.join(lines) + '\n')


def clone_spec3_parameters(pipeline):
    """Snapshot step settings plus reference overrides omitted by get_pars()."""
    from copy import deepcopy
    parameters = deepcopy(pipeline.get_pars())
    for name in pipeline.step_defs:
        step = getattr(pipeline, name)
        for reference_type in step.reference_file_types:
            override = step.get_ref_override(reference_type)
            if override is not None:
                parameters['steps'][name]['override_' + reference_type] = override
    return parameters


def subset_cube_offsets(offset_file, science_files, output_file):
    """Preserve offsets when a multi-grating association is split per cube."""
    import asdf
    names = [Path(filename).name for filename in science_files]
    with asdf.open(offset_file) as original:
        source_names = list(original['filename'])
        indices = []
        for name in names:
            matches = [i for i, candidate in enumerate(source_names) if Path(candidate).name == name]
            if len(matches) != 1:
                raise ValueError(f'Expected one cube offset for {name}; found {len(matches)}')
            indices.append(matches[0])
        tree = dict(filename=names, raoffset=[float(original['raoffset'][i]) for i in indices],
                    decoffset=[float(original['decoffset'][i]) for i in indices],
                    units=original['units'])
    with asdf.AsdfFile(tree) as output:
        output.write_to(output_file)
    return str(Path(output_file).resolve())


def build_subtracted_cube(reference_cube, association_file, pipeline_parameters,
                          s2d_files, trace_records, background_dir, target,
                          gap_policy='flag'):
    """Subtract exposure-matched backgrounds before combining the IFU images.

    reference_cube is the saved unsubtracted cube. Its actual WCS wavelength
    array is used for every background contributing to this cube. Corrected
    calibrations are temporary, with unchanged basenames to preserve cube
    offset-file lookup. Input science files and the unsubtracted cube remain
    intact. A provenance report plus background FITS/CSV spectra are retained.
    """
    from copy import deepcopy
    import tempfile
    from jwst.pipeline import Spec3Pipeline
    if gap_policy not in ('flag', 'withhold'):
        raise ValueError(f'Unknown background gap policy: {gap_policy}')
    reference_cube = Path(reference_cube).resolve()
    wavelength = cube_wavelengths(reference_cube)
    header = fits.getheader(reference_cube, 0)
    association = json.loads(Path(association_file).read_text())
    if len(association['products']) != 1:
        raise ValueError('Expected one association product per IFU cube build')
    original_product = association['products'][0]
    members = []
    for member in original_product['members']:
        if member['exptype'].lower() != 'science':
            raise ValueError('FS subtraction expects a science-only IFU association')
        filename = Path(member['expname'])
        if not filename.is_absolute():
            filename = Path(association_file).resolve().parent / filename
        member_header = fits.getheader(filename, 0)
        if all(member_header[k] == header[k] for k in ('GRATING', 'FILTER')):
            members.append(str(filename.resolve()))
    if not members:
        raise ValueError(f'No IFU association members match {reference_cube.name}')
    product = original_product['name'].replace('NOBG', 'FSBGSUB')
    if product == original_product['name']:
        product += '_FSBGSUB'
    cube_dir = reference_cube.parent
    background_dir = Path(background_dir)
    spectra_dir = background_dir / reference_cube.stem.removesuffix('_s3d')
    records, estimates = [], {}
    for filename in members:
        matches = contemporaneous_s2d(filename, s2d_files)
        key = exposure_key(filename)
        if key not in estimates:
            estimates[key] = estimate_background(matches, trace_records, wavelength)
        estimate = estimates[key]
        bg_file = save_background_spectrum(estimate, spectra_dir / (key + '_background'))
        start, end = exposure_interval(filename)
        missing = int(np.sum(~np.isfinite(estimate['brightness'])))
        records.append(dict(
            target=target, cube=str(reference_cube), ifu_file=filename, exposure=key,
            start_mjd=start, end_mjd=end, fs_files=matches, used=estimate['used'],
            background_file=bg_file, missing_bins=missing, total_bins=len(wavelength),
            min_rows=int(estimate['n_rows'].min()), max_rows=int(estimate['n_rows'].max()),
            uncertainty_missing_bins=int(np.sum(~np.isfinite(estimate['uncertainty']))),
            status='prepared', gap_policy=gap_policy,
            trace_screening_enabled=trace_records is not None,
        ))
    if (gap_policy == 'withhold' and any(r['missing_bins'] for r in records)) or all(
            r['missing_bins'] == r['total_bins'] for r in records):
        for record in records:
            record['status'] = 'cube withheld: insufficient background coverage'
        write_background_report(records, background_dir)
        print(f'FS background cube withheld; see {background_dir / "background_report.txt"}')
        return None
    write_background_report(records, background_dir)
    try:
        with tempfile.TemporaryDirectory(prefix='fsbg_', dir=cube_dir) as tmp:
            corrected_members = []
            for record in records:
                destination = Path(tmp) / Path(record['ifu_file']).name
                record['flagged_ifu_pixels'] = subtract_ifu_background(
                    record['ifu_file'], record['background_file'], destination)
                corrected_members.append({'expname': str(destination), 'exptype': 'science'})
            corrected_association = deepcopy(association)
            corrected_association['products'] = [{'name': product, 'members': corrected_members}]
            association_path = Path(tmp) / 'fs_background_asn.json'
            association_path.write_text(json.dumps(corrected_association, indent=2) + '\n')
            parameters = deepcopy(pipeline_parameters)
            parameters.update(output_dir=str(cube_dir), output_file=product, save_results=True)
            parameters['steps']['master_background']['skip'] = True  # already applied per exposure
            cube_parameters = parameters['steps'].get('cube_build', {})
            offset_file = cube_parameters.get('offset_file')
            if offset_file:
                cube_parameters['offset_file'] = subset_cube_offsets(
                    offset_file, members, Path(tmp) / 'offsets.asdf')
            corrected_pipeline = Spec3Pipeline(**parameters)
            corrected_pipeline.run(str(association_path))
    except Exception as error:
        for record in records:
            record['status'] = f'FAILED: {type(error).__name__}: {error}'
        write_background_report(records, background_dir)
        raise
    candidates = []
    for filename in cube_dir.glob(product + '*_s3d.fits'):
        candidate_header = fits.getheader(filename, 0)
        if all(candidate_header[k] == header[k] for k in ('GRATING', 'FILTER')):
            candidates.append(filename)
    if len(candidates) != 1:
        raise RuntimeError(f'Expected one background-subtracted cube for {reference_cube}; got {candidates}')
    result = candidates[0]
    with datamodels.open(result) as corrected_cube, datamodels.open(reference_cube) as original_cube:
        statuses = getattr(corrected_cube.meta.cal_step, '_instance', {})
        failed = [name for name, status in statuses.items() if status == 'FAILED']
        if failed or corrected_cube.data.shape != original_cube.data.shape:
            reason = f'failed steps {failed}' if failed else 'cube shape mismatch'
            for record in records:
                record['status'] = 'FAILED: ' + reason
            write_background_report(records, background_dir)
            raise RuntimeError(f'Background-subtracted cube validation failed: {reason}')
    if not np.array_equal(cube_wavelengths(result), wavelength):
        for record in records:
            record['status'] = 'FAILED: cube wavelength grid mismatch'
        write_background_report(records, background_dir)
        raise RuntimeError(f'Background-subtracted cube grid differs from {reference_cube.name}')
    for record in records:
        record['status'] = 'subtracted; gaps flagged' if record['missing_bins'] else 'subtracted'
        record['subtracted_cube'] = str(result.resolve())
    write_background_report(records, background_dir)
    return str(result.resolve())

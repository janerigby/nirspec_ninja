"""Compare exposure-level RATE products using the retained background mask.

Only NumPy and Astropy are needed; this module does not import JWST.
Mask SCI values of 1 mean background, 0 mean excluded/illuminated.
"""
from pathlib import Path
import csv
import fcntl
import numpy as np
from astropy.io import fits

LEGACY_COLUMNS = ('Rate_filename', 'RMS_before', 'RMS_after', 'MAD_sigma_before',
           'MAD_sigma_after', 'RMS_reduction_percent', 'N_pixels', 'Units')

COLUMNS = LEGACY_COLUMNS + ('Fit_method', 'MAD_reduction_percent',
    'EFFEXPTM', 'READPATT', 'NGROUPS', 'NINTS', 'NFRAMES', 'GROUPGAP',
    'TGROUP', 'TFRAME', 'DETECTOR', 'SUBARRAY')


HEADINGS = dict(zip(COLUMNS, ('File', 'RMS_pre', 'RMS_post', 'MAD_pre', 'MAD_post',
    'RMS_%', 'Npix', 'Units', 'Method', 'MAD_%', 'Exp_s', 'Readout', 'Ngroups',
    'Nints', 'Nframes', 'Gap', 'Group_s', 'Frame_s', 'Det', 'Subarray')))


def measure_noise(before_file, after_file, mask_file, fit_method=None):
    """Measure both 2D SCI images on identical finite, usable background pixels.

    RMS is sqrt(mean((x - median(x))**2)); MAD sigma is
    1.4826022185 * median(abs(x - median(x))). No sigma clipping is done:
    rejection based on the cleaned image could exaggerate improvement.
    DQ DO_NOT_USE (bit 1) and NON_SCIENCE (bit 512) are excluded in either
    image. Other flags, e.g. JUMP_DET, need not invalidate a fitted slope.
    Values retain SCI BUNIT; inconsistent units/shapes fail explicitly.
    """
    with fits.open(before_file, memmap=False) as pre, \
         fits.open(after_file, memmap=False) as post, \
         fits.open(mask_file, memmap=False) as mask:
        header = dict(pre[0].header)
        header.update(dict(pre['SCI'].header))
        header.update(dict(post[0].header))
        header.update(dict(post['SCI'].header))
        metadata = {key: header.get(key, 'unknown') for key in COLUMNS[10:]}
        a = np.asarray(pre['SCI'].data, dtype=float)
        b = np.asarray(post['SCI'].data, dtype=float)
        m = np.asarray(mask['SCI'].data)
        if a.ndim != 2 or a.shape != b.shape or a.shape != m.shape:
            raise ValueError('RATE images and background mask must have identical 2D shapes')
        units = pre['SCI'].header.get('BUNIT', '').strip()
        if units != post['SCI'].header.get('BUNIT', '').strip():
            raise ValueError('Before/after SCI units differ')
        valid = (m == 1) & np.isfinite(a) & np.isfinite(b)
        for hdus in (pre, post):
            if 'DQ' in hdus:
                dq = hdus['DQ'].data
                if dq.shape != a.shape:
                    raise ValueError('DQ and SCI shapes differ')
                valid &= (dq & (1 | 512)) == 0
        count = int(valid.sum())
        if count < 2:
            raise ValueError('Fewer than two usable background pixels')
        def scatter(data):
            residual = data[valid] - np.median(data[valid])
            return float(np.sqrt(np.mean(residual**2))), float(1.4826022185 * np.median(np.abs(residual)))
        rms_before, mad_before = scatter(a)
        rms_after, mad_after = scatter(b)
    row = dict(zip(LEGACY_COLUMNS, (Path(after_file).name, rms_before, rms_after,
                mad_before, mad_after,
                100 * (1 - rms_after / rms_before) if rms_before > 0 else float('nan'),
                count, units or 'unspecified')))
    method = str(fit_method).lower() if fit_method is not None else 'unknown'
    if method not in ('median', 'fft', 'unknown'):
        raise ValueError('fit_method must be median or fft')
    row.update(Fit_method='FFT' if method == 'fft' else method,
        MAD_reduction_percent=100 * (1 - mad_after / mad_before) if mad_before > 0 else float('nan'))
    row.update(metadata)
    return row


def _format_report_row(row):
    item = dict(row)
    for key in ('RMS_before', 'RMS_after', 'MAD_sigma_before', 'MAD_sigma_after'):
        item[key] = format(float(item[key]), '.4g')
    for key in ('RMS_reduction_percent', 'MAD_reduction_percent'):
        item[key] = format(float(item[key]), '.2f')
    for key in ('EFFEXPTM', 'TGROUP', 'TFRAME'):
        if item[key] != 'unknown':
            item[key] = format(float(item[key]), '.6g')
    return item


def _write_aligned_table(stream, rows):
    """Pad columns for people while retaining tab-separated machine parsing."""
    items = [_format_report_row(row) for row in rows]
    widths = {key: max(len(HEADINGS[key]), *(len(str(item[key])) for item in items))
              for key in COLUMNS}
    stream.write('  '.join(HEADINGS[key].ljust(widths[key]) for key in COLUMNS) + '\n')
    text_columns = {'Rate_filename', 'Units', 'Fit_method', 'READPATT', 'DETECTOR', 'SUBARRAY'}
    for item in items:
        stream.write('  '.join(str(item[key]).ljust(widths[key]) if key in text_columns
                     else str(item[key]).rjust(widths[key]) for key in COLUMNS) + '\n')


def write_noise_report(before_file, after_file, mask_file, report_file=None, fit_method=None):
    """Upsert one exposure in a tab-separated directory report (process-safe).

    Positive reduction means improved RMS; negative means worse. Re-running
    a filename replaces its row. The lock is separate from the atomically
    replaced report so parallel LEGGOS workers cannot lose each other's rows.
    """
    row = measure_noise(before_file, after_file, mask_file, fit_method=fit_method)
    report = Path(report_file) if report_file else Path(after_file).parent / 'flicker_noise_rms.txt'
    with report.with_suffix(report.suffix + '.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        rows = {}
        if report.exists():
            with report.open(newline='') as stream:
                lines = [line for line in stream if not line.startswith('#')]
                if lines and '\t' not in lines[0]:
                    import re
                    lines = ['\t'.join(re.split(r' {2,}', line.rstrip())) + '\n' for line in lines]
                reader = csv.DictReader(lines, delimiter='\t')
                reverse_headings = {value: key for key, value in HEADINGS.items()}
                reader.fieldnames = [reverse_headings.get(name.strip(), name.strip()) for name in (reader.fieldnames or [])]
                if tuple(reader.fieldnames) not in (COLUMNS, LEGACY_COLUMNS):
                    raise ValueError(f'Unexpected columns in {report}')
                rows = {}
                for item in reader:
                    item = {key: value.strip() for key, value in item.items()}
                    rows[item['Rate_filename']] = item
                for item in rows.values():
                    for key in COLUMNS:
                        item.setdefault(key, 'unknown')
                    if item['MAD_reduction_percent'] == 'unknown':
                        before, after = float(item['MAD_sigma_before']), float(item['MAD_sigma_after'])
                        item['MAD_reduction_percent'] = 100 * (1 - after / before) if before > 0 else float('nan')
        rows[row['Rate_filename']] = row
        temporary = report.with_suffix(report.suffix + '.tmp')
        with temporary.open('w', newline='') as stream:
            stream.write('# Exposure-level RATE comparison; mask=1 is non-illuminated background.\n')
            stream.write('# RMS about each median; MAD_sigma=1.4826022185*MAD; identical finite/DQ-valid pixels.\n')
            stream.write('# Positive RMS/MAD reduction indicates improvement. Noise: 4 significant digits; percentages: 2 decimals.\n')
            _write_aligned_table(stream, [rows[name] for name in sorted(rows)])
        temporary.replace(report)
    return row


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('before_file')
    parser.add_argument('after_file')
    parser.add_argument('mask_file')
    parser.add_argument('--report-file')
    parser.add_argument('--fit-method', choices=('median', 'fft'), help='Actual cleaning method; unknown if omitted')
    args = parser.parse_args()
    print(write_noise_report(**vars(args)))

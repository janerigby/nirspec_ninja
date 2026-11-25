"""Shared spatial profiles, trace masks, and non-flatness measurements; no plotting."""
import numpy as np
import warnings
from scipy.special import ndtri
from nirspec_ninja.subtract_background import fs_background as bg


def profile(slit, clip_nan=5):
    good, excluded = bg.spatial_pixel_mask(slit, clip_nan=clip_nan)
    rows = np.flatnonzero(~excluded)
    if len(rows) < 4:
        raise ValueError('Too few usable spatial rows')
    lo = max(np.min(slit.wavelength[r][good[r]]) for r in rows)
    hi = min(np.max(slit.wavelength[r][good[r]]) for r in rows)
    if hi <= lo:
        raise ValueError('No common wavelength range across usable rows')
    # Preserve missing rows between valid bands so plot lines break at gaps.
    display_rows = np.arange(rows[0], rows[-1] + 1)
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        full = np.nanmedian(np.where(good, slit.brightness, np.nan), axis=1)[display_rows]
        common = good & (slit.wavelength >= lo) & (slit.wavelength <= hi)
        matched = np.nanmedian(np.where(common, slit.brightness, np.nan), axis=1)[display_rows]
    full[excluded[display_rows]] = np.nan
    matched[excluded[display_rows]] = np.nan
    # Coordinate fractions are only a visual comparison, not a sky registration.
    fraction = (display_rows - rows[0]) / (rows[-1] - rows[0])
    return fraction, full, matched, (float(lo), float(hi))


def spatial_flatness(fraction, values):
    """Largest 3-row-median departure; a descriptive score, not significance."""
    finite = np.isfinite(fraction) & np.isfinite(values)
    if finite.sum() < 4:
        return dict(status='unavailable', reason='Fewer than four valid rows')
    baseline = float(np.nanmedian(values[finite]))
    # Only contiguous, fully valid windows; do not join missing spatial rows.
    windows = np.lib.stride_tricks.sliding_window_view(values, 3)
    centers = fraction[1:-1]
    valid = np.all(np.isfinite(windows), axis=1) & np.isfinite(centers)
    steps = np.diff(fraction)
    positive = steps[np.isfinite(steps) & (steps > 0)]
    if positive.size:
        valid &= (steps[:-1] <= positive.min() * 1.01) & (steps[1:] <= positive.min() * 1.01)
    if valid.sum() < 2:
        return dict(status='unavailable', reason='Fewer than two contiguous three-row windows')
    smoothed = np.median(windows[valid], axis=1)
    departure = smoothed - baseline
    index = int(np.argmax(np.abs(departure)))
    delta = float(departure[index])
    return dict(status='ok', baseline_MJy_sr=baseline, peak_departure_MJy_sr=delta,
                peak_fraction=float(centers[valid][index]), window_rows=3,
                fractional_departure=None if abs(baseline) <= 1e-12 else abs(delta) / abs(baseline))


def diagnostic_trace(slit, threshold=2.5, enclosed_flux=.99, clip_nan=5, use_error=True):
    """Detect the strongest continuum trace and enclose 99% of Gaussian flux."""
    from nirspec_ninja.fixed_slit import source_detection as detector
    good, excluded = bg.spatial_pixel_mask(slit, clip_nan=clip_nan)
    usable_rows = int(np.count_nonzero(~excluded))
    if usable_rows < 4:
        passing_before_margin = int(np.count_nonzero(2 * np.sum(good, axis=1) > good.shape[1]))
        return dict(status='unavailable', reason=(
            f'Insufficient usable data: {passing_before_margin} rows pass the >50% valid-pixel rule; '
            f'{usable_rows} remain after the {clip_nan}-row exclusion margin (need at least 4)'),
            detected=False, significance=None, raw_significance=None, center=None, fwhm=None,
            excluded_rows=np.flatnonzero(excluded).tolist(), usable_rows=usable_rows,
            rows_before_margin=passing_before_margin, threshold_sigma=threshold,
            enclosed_source_flux_fraction=enclosed_flux)
    result = detector.detect_source_s2d(
        np.where(good & ~excluded[:, None], slit.brightness, np.nan),
        error_2d=slit.error if use_error else None, dq_2d=slit.dq, allowed_source_center_rows=~excluded,
        slit_length_pixels=detector.NIRSPEC_FIXED_SLIT_LENGTH_ARCSEC[slit.name]
        / detector.NIRSPEC_SPATIAL_PIXEL_SCALE_ARCSEC)
    record = {key: result.get(key) for key in
              ['status', 'reason', 'significance', 'raw_significance', 'center', 'fwhm']}
    for key in ['significance', 'raw_significance', 'center', 'fwhm']:
        value = record[key]
        record[key] = float(value) if value is not None and np.isfinite(value) else None
    record['threshold_sigma'] = threshold
    record['detected'] = (record['status'] == 'ok' and record['significance'] is not None
                          and record['significance'] >= threshold)
    record['excluded_rows'] = np.flatnonzero(excluded).tolist()
    record['blocks_above_3sigma'] = result.get('blocks_above_3sigma')
    record['enclosed_source_flux_fraction'] = enclosed_flux
    if record['detected']:
        sigma = record['fwhm'] / np.sqrt(8 * np.log(2))
        half_width = sigma * ndtri((1 + enclosed_flux) / 2)
        record['bounds_flux_rows'] = [record['center'] - half_width, record['center'] + half_width]
    return record


def trace_row_mask(trace, nrows):
    """Exclude row centers within the detected source's enclosed-flux bounds."""
    if not trace.get('detected'):
        return np.zeros(nrows, dtype=bool)
    lo, hi = trace['bounds_flux_rows']
    rows = np.arange(nrows)
    return (rows >= lo) & (rows <= hi)

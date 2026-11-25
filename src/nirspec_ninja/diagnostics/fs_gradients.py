"""Reusable L2/L3 fixed-slit gradient diagnostics; no source detection."""
from pathlib import Path
import json
import csv
import re
import warnings
import numpy as np
from scipy.interpolate import PchipInterpolator, make_smoothing_spline
from astropy.io import fits
from astropy.time import Time
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from nirspec_ninja.subtract_background import fs_background as bg
from nirspec_ninja.fixed_slit.spatial_diagnostics import (
    profile, spatial_flatness, diagnostic_trace, trace_row_mask,
)

# Diagnostic selection is independent of the slits allowed for sky subtraction.
# Left to right in the FOV, matching the finder charts.
GRADIENT_SLITS = ('S200A2', 'S400A1', 'S1600A1', 'S200A1')


def good_slit_mark(mark):
    return mark is True or mark == -1


def create_bkg_override(path, object_name, gratings):
    """Create an all-use starter once; never overwrite previous user choices."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        with path.open('x') as handle:
            handle.write('# Starter evaluations: review the slit diagnostics before using as sky.\n'
                         '# 1 = use with trace masking; 0 = DO NOT USE; -1 = use whole valid slit.\n')
            handle.write('target grating ' + ' '.join(GRADIENT_SLITS) + '\n')
            for grating in sorted(set(gratings)):
                handle.write(f'{object_name} {grating} 1 1 1 1\n')
    except FileExistsError:
        return False
    return True


def read_crowding_notes(path, target, grating=None, fallback_target=None):
    """Read target/grating/slit 0/1 evaluations (also accepts legacy target-only tables).

    Target matching ignores case and periods (Cosmic.eye == COSMIC.EYE).
    Missing targets/slit columns remain unmarked; malformed tables fail loudly.
    """
    if path is None:
        return {}, ''
    normalize = lambda value: value.casefold().replace('.', '')
    exact_target = normalize(target)
    base_target = normalize(re.sub(r'[-_.]?IFU[-_.]?POS\d+$', '', target, flags=re.IGNORECASE))
    header = None
    selected = None
    fallback = None
    seen = set()
    with Path(path).expanduser().open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split('#', 1)[0].split()
            if not fields:
                continue
            if header is None:
                header = fields
                if header[0].casefold() != 'target':
                    raise ValueError(f'{path}:{number}: expected target table header')
                by_grating = len(header) > 1 and header[1].casefold() == 'grating'
                offset = 2 if by_grating else 1
                if by_grating and grating is None:
                    raise ValueError(f'{path}: grating is required for this table')
                columns = header[offset:-1] if header[-1].casefold() == 'notes' else header[offset:]
                if not columns or len(set(columns)) != len(columns) or any(
                        name not in GRADIENT_SLITS for name in columns):
                    raise ValueError(f'{path}:{number}: invalid slit columns')
                continue
            count = len(columns)
            if len(fields) < count + offset or any(value not in ('-1', '0', '1') for value in fields[offset:count + offset]):
                raise ValueError(f'{path}:{number}: each slit must be marked -1, 0 or 1')
            key = (normalize(fields[0]), fields[1].casefold() if by_grating else None)
            if key in seen:
                raise ValueError(f'{path}:{number}: duplicate target/grating {key}')
            seen.add(key)
            if not by_grating or key[1] == str(grating).casefold():
                evaluation = (dict(zip(columns, (-1 if value == '-1' else value == '1'
                                                for value in fields[offset:count + offset]))),
                              ' '.join(fields[count + offset:]))
                if key[0] == exact_target:
                    selected = evaluation
                elif key[0] == base_target:
                    fallback = evaluation
    if header is None:
        raise ValueError(f'{path}: empty crowding table')
    if selected is not None:
        return selected
    if fallback is not None:
        return fallback
    if fallback_target is not None and normalize(fallback_target) != exact_target:
        return read_crowding_notes(path, fallback_target, grating=grating)
    return {}, ''



def median_good_profiles(profiles, marks):
    """Median across good slits, then across normalized rows: a constant sky level."""
    good = [(name, x, values) for name, x, values in profiles if good_slit_mark(marks.get(name))]
    grid = np.linspace(0, 1, 101)
    if not good:
        return grid, None, []
    sampled = []
    for name, x, values in good:
        finite = np.isfinite(x) & np.isfinite(values)
        if finite.sum() >= 2:
            sampled.append((name, interpolate_valid_runs(x, values, grid)))
    if not sampled:
        return grid, None, []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        row_median = np.nanmedian([values for _, values in sampled], axis=0)
        median = np.full(grid.shape, np.nanmedian(row_median))
        median[~np.isfinite(row_median)] = np.nan
    return grid, median, [name for name, _ in sampled]




def interpolate_valid_runs(wave, values, grid):
    """Interpolate contiguous valid runs without bridging masked spectral gaps."""
    result = np.full(grid.shape, np.nan)
    indices = np.flatnonzero(np.isfinite(wave) & np.isfinite(values))
    for run in np.split(indices, np.flatnonzero(np.diff(indices) > 1) + 1):
        if not len(run):
            continue
        order = np.argsort(wave[run])
        x, index = np.unique(wave[run][order], return_index=True)
        y = values[run][order][index]
        if len(x) < 2:
            selected = np.isclose(grid, x[0], rtol=0, atol=1e-12)
            result[selected] = y[0]
        else:
            selected = (grid >= x[0]) & (grid <= x[-1])
            result[selected] = np.interp(grid[selected], x, y)
    return result


def spatial_median_spectrum(slit, extra_excluded_rows=None):
    """Align usable S2D rows in wavelength, then take their spatial median."""
    good, excluded = bg.spatial_pixel_mask(slit)
    if extra_excluded_rows is not None:
        excluded |= np.asarray(extra_excluded_rows, dtype=bool)
    rows = np.flatnonzero(~excluded)
    if not len(rows):
        raise ValueError('No usable spatial rows for spectrum')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        grid = np.nanmedian(slit.wavelength[rows], axis=0)
        grid = np.unique(grid[np.isfinite(grid)])
        aligned = [interpolate_valid_runs(slit.wavelength[row],
                   np.where(good[row], slit.brightness[row], np.nan), grid) for row in rows]
        values = np.nanmedian(aligned, axis=0)
    if not np.any(np.isfinite(values)):
        raise ValueError('No usable spectrum')
    return grid, values


def combined_median_spectrum(spectra):
    """Equal slit weights on the coarsest input grid; require all slits' coverage."""
    if not spectra:
        return np.array([]), np.array([])
    reference = max(spectra, key=lambda item: np.median(np.diff(item[1])))
    lo = max(wave[np.isfinite(values)].min() for _, wave, values in spectra)
    hi = min(wave[np.isfinite(values)].max() for _, wave, values in spectra)
    grid = reference[1][(reference[1] >= lo) & (reference[1] <= hi)]
    aligned = np.array([interpolate_valid_runs(wave, values, grid) for _, wave, values in spectra])
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        median = np.nanmedian(aligned, axis=0)
    median[~np.all(np.isfinite(aligned), axis=0)] = np.nan
    return grid, median


def interpolate_background_gaps(wave, values):
    """Fill interior display gaps with shape-preserving cubic interpolation."""
    finite = np.isfinite(wave) & np.isfinite(values)
    filled = values.copy()
    interpolated = np.zeros(values.shape, dtype=bool)
    if finite.sum() < 2:
        return filled, interpolated
    x, y = wave[finite], values[finite]
    interpolated = np.isfinite(wave) & ~np.isfinite(values) & (wave > x.min()) & (wave < x.max())
    filled[interpolated] = PchipInterpolator(x, y, extrapolate=False)(wave[interpolated])
    return filled, interpolated


def smooth_background_spline(wave, values, scale_um=.5):
    """Robust, wavelength-based continuum fit across measured chips and gaps.

    Median bins of scale/5 suppress spikes. A curvature-penalized cubic spline
    uses a ~scale half-response period; Huber reweighting reduces anchor outliers.
    Missing pixels supply no fit constraints. Never extrapolate beyond coverage.
    """
    if not np.isfinite(scale_um) or scale_um <= 0:
        raise ValueError('Grating smoothing scale must be positive and finite')
    finite = np.isfinite(wave) & np.isfinite(values)
    fitted = np.full(values.shape, np.nan)
    gaps = np.zeros(values.shape, dtype=bool)
    if finite.sum() < 2:
        return fitted, gaps
    x, y = wave[finite], values[finite]
    bins = np.floor((x - x.min()) / (scale_um / 5)).astype(int)
    ids = np.unique(bins)
    anchors_x = np.array([np.median(x[bins == index]) for index in ids])
    anchors_y = np.array([np.median(y[bins == index]) for index in ids])
    covered = np.isfinite(wave) & (wave >= x.min()) & (wave <= x.max())
    if len(anchors_x) < 5:
        fitted[covered] = PchipInterpolator(x, y, extrapolate=False)(wave[covered])
    else:
        spacing = np.median(np.diff(anchors_x))
        penalty = (scale_um / (2 * np.pi))**4 / spacing
        weights = np.ones(len(anchors_x))
        for _ in range(4):
            spline = make_smoothing_spline(anchors_x, anchors_y, w=weights, lam=penalty)
            residual = anchors_y - spline(anchors_x)
            scatter = max(1.4826 * np.median(np.abs(residual - np.median(residual))),
                          1e-8 * max(np.median(np.abs(anchors_y)), 1e-6))
            weights = np.minimum(1., 1.5 * scatter / np.maximum(np.abs(residual), scatter * 1e-8))
        fitted[covered] = spline(wave[covered])
    gaps = covered & ~finite
    return fitted, gaps


def predicted_jbt_background(header, cache_dir, models):
    """JBT total at the target and UTC observation midpoint, Jan 1 = day 0."""
    ra, dec = float(header['TARG_RA']), float(header['TARG_DEC'])
    if 'EXPMID' in header:
        midpoint = Time(float(header['EXPMID']), format='mjd', scale='utc')
    elif 'EXPSTART' in header and 'EXPEND' in header:
        midpoint = Time((float(header['EXPSTART']) + float(header['EXPEND'])) / 2,
                        format='mjd', scale='utc')
    else:
        midpoint = Time(header['DATE-OBS'] + 'T' + header.get('TIME-OBS', '00:00:00'), scale='utc')
    day = midpoint.to_datetime().timetuple().tm_yday - 1
    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache = cache_dir / f'jbt_{ra:.8f}_{dec:.8f}_day{day:03d}.npz'
    if cache.exists():
        with np.load(cache, allow_pickle=False) as saved:
            wave, total = saved['wave'].copy(), saved['total'].copy()
            version = str(saved['version'].item())
    else:
        from jwst_backgrounds import jbt
        key = (ra, dec)
        if key not in models:
            models[key] = jbt.background(ra, dec, 2.0)
        model = models[key]
        matching = np.flatnonzero(np.asarray(model.bkg_data['calendar']) == day)
        if not matching.size:
            raise ValueError(f'Observation day {day} (Jan 1=0) unavailable in JBT')
        wave = np.asarray(model.bkg_data['wave_array'], dtype=float)
        total = np.asarray(model.bkg_data['total_bg'][matching[0]], dtype=float)
        version = str(model.cache_version)
        np.savez(cache, wave=wave, total=total, version=version)
    return wave, total, dict(ra_deg=ra, dec_deg=dec, observation_midpoint_utc=midpoint.isot,
                             day_zero_based=day, cache_version=version, cache_file=str(cache),
                             units='MJy/sr', status='ok')






def make_fs_gradient_report(fs_root, output_dir, target, crowding_notes=None, jbt_background=True,
                            include_l2=False, detect_traces=True, use_bkg_override=None,
                            grating_smoothing_um=.5):
    """Plot all calibrated fixed slits; human bad marks exclude only candidate sky."""
    selections = {}
    def selection(grating):
        key = str(grating).upper()
        if key not in selections:
            selections[key] = read_crowding_notes(crowding_notes, ROOT.parent.name, grating=key,
                                                  fallback_target=target)
            marks, _ = selections[key]
            if crowding_notes is not None:
                good_names = [name for name in GRADIENT_SLITS if good_slit_mark(marks.get(name))]
                bad_names = [name for name in GRADIENT_SLITS if marks.get(name) is False]
                print(f'Crowding selection {ROOT.parent.name}/{key}: good={good_names}, DO NOT USE={bad_names}', flush=True)
        return selections[key]
    ROOT, OUT = Path(fs_root), Path(output_dir)
    OUT.mkdir(parents=True, exist_ok=True)
    l2 = {}
    for filename in sorted(ROOT.glob('L2b/*_s2d.fits')) if include_l2 else []:
        header = fits.getheader(filename)
        for slit in bg.read_slit_pixels(filename):
            if slit.name in GRADIENT_SLITS:
                if not slit.calibrated:
                    raise ValueError(f'Expected MJy/sr brightness in {filename}, {slit.name}')
                key = (header['GRATING'], header['FILTER'], slit.name)
                try:
                    measured = profile(slit)
                except ValueError as error:
                    print(f'Gradient diagnostic skipped {filename.name}/{slit.name}: {error}')
                    continue
                l2.setdefault(key, []).append((filename, measured))
    all_candidates = sorted(p for p in ROOT.glob(f'L3/Level3_{target.replace(".", "")}_FS_*_s2d.fits')
                            if any('_' + n.lower() + '_' in p.name for n in GRADIENT_SLITS))
    records = []
    candidates = all_candidates
    if not all_candidates:
        (OUT / 'gradient_summary.json').write_text(json.dumps(records, indent=2) + '\n')
        print(f'No eligible L3 S2Ds for {target}: no gradient plot produced')
        return None
    override_path = OUT / 'FS_bkg_override.txt'
    created = create_bkg_override(override_path, ROOT.parent.name,
                                  [fits.getheader(p)['GRATING'] for p in all_candidates])
    print(f'{"Created" if created else "Preserved"} background override table: {override_path.resolve()}', flush=True)
    if use_bkg_override is True or (use_bkg_override is None and crowding_notes is None):
        crowding_notes = override_path
    print('Review FS_bkg_override.txt; edit 1/0/-1 choices and rerun with --use-bkg-overrides '
          '(or use_bkg_override=True).', flush=True)
    # Read bad L3 profiles too, solely for the de-emphasized comparison overlay.
    comparisons = {}
    background_comparisons = {}
    spectra_by_configuration = {}
    traces = {}
    for filename in all_candidates:
        header = fits.getheader(filename)
        marks, human_notes = selection(header['GRATING'])
        slit = next(bg.read_slit_pixels(filename))
        if not slit.calibrated:
            continue
        try:
            x, _, values, _ = profile(slit)
        except ValueError:
            continue
        comparisons.setdefault((header['GRATING'], header['FILTER']), []).append((slit.name, x, values))
        trace = diagnostic_trace(slit) if detect_traces else dict(status='disabled', detected=False)
        override = marks.get(slit.name) == -1
        trace['whole_slit_override'] = override
        mask = np.zeros(slit.brightness.shape[0], dtype=bool) if override else trace_row_mask(trace, slit.brightness.shape[0])
        trace['masked_source_rows'] = np.flatnonzero(mask).tolist()
        traces[filename] = trace
        _, excluded = bg.spatial_pixel_mask(slit)
        usable = np.flatnonzero(~excluded)
        masked_values = values.copy()
        masked_values[mask[usable[0]:usable[-1] + 1]] = np.nan
        background_comparisons.setdefault((header['GRATING'], header['FILTER']), []).append(
            (slit.name, x, masked_values))
        if good_slit_mark(marks.get(slit.name)):
            try:
                wave, spectrum = spatial_median_spectrum(slit, extra_excluded_rows=mask)
            except ValueError as error:
                print(f'Spectrum diagnostic skipped {filename.name}: {error}')
            else:
                spectra_by_configuration.setdefault((header['GRATING'], header['FILTER']), []).append(
                    (slit.name, wave, spectrum))
    n_panels = len(candidates) + len(comparisons)
    if not n_panels:
        (OUT / 'gradient_summary.json').write_text(json.dumps(records, indent=2) + '\n')
        return None
    files_by_configuration = {}
    for filename in candidates:
        header = fits.getheader(filename)
        marks, human_notes = selection(header['GRATING'])
        configuration = (header['GRATING'], header['FILTER'])
        files_by_configuration.setdefault(configuration, []).append(filename)
    configurations = sorted(set(comparisons) | set(files_by_configuration))
    nrows = len(configurations)
    ncols = max((len(files) for files in files_by_configuration.values()), default=0) + 2
    figure = Figure(figsize=(4.7 * ncols, 4 * nrows))
    FigureCanvasAgg(figure)
    axes = figure.subplots(nrows, ncols, squeeze=False)
    slit_axes = {}
    comparison_axes = {}
    for row, configuration in enumerate(configurations):
        files = files_by_configuration.get(configuration, [])
        files = sorted(files, key=lambda filename: next(
            index for index, name in enumerate(GRADIENT_SLITS)
            if '_' + name.lower() + '_' in filename.name))
        for column, filename in enumerate(files):
            slit_axes[filename] = axes[row, column]
        for column in range(len(files), ncols - 2):
            axes[row, column].set_visible(False)
        comparison_axes[configuration] = axes[row, -2]
        if configuration not in comparisons:
            axes[row, -2].set_visible(False)
    axes_by_configuration = {}
    orange_by_configuration = {}
    for filename in candidates:
        axis = slit_axes[filename]
        header = fits.getheader(filename)
        marks, human_notes = selection(header['GRATING'])
        slit = next(bg.read_slit_pixels(filename))
        if not slit.calibrated:
            raise ValueError(f'Expected MJy/sr brightness in {filename}, {slit.name}')
        key = (header['GRATING'], header['FILTER'], slit.name)
        try:
            fraction, full, matched, interval = profile(slit)
        except ValueError as error:
            axis.set(title='/'.join(key))
            axis.text(.5, .5, str(error), ha='center', transform=axis.transAxes)
            records.append(dict(file=str(filename), status='skipped', reason=str(error)))
            continue
        upstream = []
        for index, (source, (x, _, values, wave_range)) in enumerate(l2.get(key, [])):
            axis.plot(x, values, color='gray', alpha=.35, lw=.8,
                      label='L2 exposures: common wavelength range' if index == 0 else None)
            upstream.append(dict(file=str(source), range_um=wave_range,
                                 end_difference=float(np.nanmedian(values[x >= .75]) - np.nanmedian(values[x <= .25]))))
        bad = marks.get(slit.name) is False
        axis.plot(fraction, matched, color='crimson' if bad else 'tab:orange',
                  label='L3: common wavelength range')
        override = marks.get(slit.name) == -1
        human_mark = 'bad' if bad else ('whole_slit' if override else ('good' if good_slit_mark(marks.get(slit.name)) else 'unmarked'))
        title = '/'.join(key) + (' — DO NOT USE' if bad else '')
        if bad:
            axis.set_facecolor('#fff0f0')
            axis.title.set_color('crimson')
            for spine in axis.spines.values():
                spine.set_edgecolor('crimson')
        trace = traces.get(filename, dict(status='unavailable', detected=False))
        metric_values = next(values for name, _, values in background_comparisons[key[:2]]
                             if name == slit.name)
        flatness = spatial_flatness(fraction, metric_values)
        flatness['trace_masked'] = bool(trace.get('detected')) and not override
        score = flatness.get('fractional_departure')
        text = f'Spatial non-flatness: {100 * score:.1f}%' if score is not None else 'Spatial non-flatness: undefined'
        # Keep the score outside the data/legend area; put it on its own
        # line below the spatial-coordinate label.
        if trace['detected']:
            # Both axes use a linear transform of original S2D row numbers,
            # including missing rows; trace coordinates are never compressed.
            _, excluded = bg.spatial_pixel_mask(slit)
            rows = np.flatnonzero(~excluded)
            origin, span = rows[0], rows[-1] - rows[0]
            center = (trace['center'] - origin) / span
            lower, upper = [(bound - origin) / span for bound in trace['bounds_flux_rows']]
            axis.axvline(center, color='teal', ls=':', lw=1.8,
                         label=f"Trace center ({trace['significance']:.1f}σ)")
            if not override:
                axis.axvspan(lower, upper, color='teal', alpha=.12,
                             label='Trace: 99% of source flux')
            else:
                title += ' — whole slit override'
        elif detect_traces and trace['status'] != 'ok':
            title += '\nTrace fit failed'
        axis.set(title=title, xlabel='Linear spatial row fraction\n' + text,
                 ylabel='Median surface brightness [MJy/sr]')
        axis.xaxis.label.set_fontsize(12)
        axis.xaxis.label.set_fontweight('normal')
        axis.set_xlim(0, 1)
        axis.legend(fontsize=7)
        # All slit panels share the stretch selected from good/unmarked data.
        axes_by_configuration.setdefault(key[:2], []).append(axis)
        if not bad:
            orange_by_configuration.setdefault(key[:2], []).append(matched[np.isfinite(matched)])
        records.append(dict(file=str(filename), range_um=interval,
                            grating=key[0], filter=key[1], slit=key[2],
                            human_mark=human_mark, human_notes=human_notes,
                            excluded_from_background=bad, spatial_flatness=flatness,
                            trace=trace,
                            l3_end_difference=float(np.nanmedian(matched[fraction >= .75]) - np.nanmedian(matched[fraction <= .25])),
                            l2=upstream))
    comparison_records = []
    for configuration, profiles in sorted(comparisons.items()):
        marks, human_notes = selection(configuration[0])
        axis = comparison_axes[configuration]
        profiles = sorted(background_comparisons[configuration], key=lambda item: GRADIENT_SLITS.index(item[0]))
        grid, median, good_names = median_good_profiles(background_comparisons[configuration], marks)
        for name, x, values in profiles:
            if good_slit_mark(marks.get(name)):
                axis.plot(x, values, lw=1.2, label=name)
            elif marks.get(name) is False:
                axis.plot(x, values, color='crimson', alpha=.25, lw=.8, ls='--',
                          label=f'{name}: DO NOT USE', zorder=1)
        if median is not None:
            label = 'Median bkg (traces masked)' if detect_traces else 'Median bkg'
            axis.plot(grid, median, color='black', lw=2.5, label=label, zorder=5)
            # Only explicitly good slits determine the comparison stretch.
            good_values = np.concatenate([values[np.isfinite(values)] for name, _, values in profiles
                                          if name in good_names])
            lower, upper = float(good_values.min()), float(good_values.max())
            padding = .05 * (upper - lower) if upper > lower else max(abs(lower) * .05, 1e-6)
            axis.set_ylim(lower - padding, upper + padding)
        else:
            marked_good = [name for name in GRADIENT_SLITS if good_slit_mark(marks.get(name))]
            message = ('No matching good slit evaluations\nCheck target/grating in crowding table'
                       if not marked_good else 'Good slits have insufficient usable rows\nafter pixel/trace masking')
            axis.text(.5, .5, message, ha='center', transform=axis.transAxes)
        axis.set(title='/'.join(configuration) + ': median bkg',
                 xlabel='Linear spatial row fraction',
                 ylabel='Median surface brightness [MJy/sr]')
        handles, _ = axis.get_legend_handles_labels()
        if handles:
            axis.legend(fontsize=7)
        comparison_records.append(dict(configuration=list(configuration), good_slits=good_names,
                                       marked_good_slits=[name for name in GRADIENT_SLITS if good_slit_mark(marks.get(name))],
                                       remaining_spatial_rows={name: int(np.isfinite(values).sum())
                                                               for name, _, values in profiles},
                                       sky_level=None if median is None else float(np.nanmedian(median)),
                                       fraction=grid.tolist(),
                                       median_good=None if median is None else
                                       [float(v) if np.isfinite(v) else None for v in median],
                                       note='Median across good slits, then across normalized spatial rows. '
                                            'Regions enclosing 99% of fitted Gaussian source flux are excluded. '
                                            'Equal slit weights; interpolated normalized rows, not sky registration. '
                                            'Each slit uses its own common wavelength interval.'))
    jbt_models = {}
    jbt_records = []
    for row, configuration in enumerate(configurations):
        axis = axes[row, -1]
        spectra = sorted(spectra_by_configuration.get(configuration, []),
                         key=lambda item: GRADIENT_SLITS.index(item[0]))
        is_grating = configuration[0].upper() != 'PRISM'
        for name, wave, values in spectra:
            if is_grating:
                fitted, _ = smooth_background_spline(wave, values, grating_smoothing_um)
                displayed = np.where(np.isfinite(values), fitted, np.nan)
            else:
                displayed = values
            axis.plot(wave, displayed, lw=1, label=name)
        wave, median = combined_median_spectrum(spectra)
        if len(wave) and np.any(np.isfinite(median)):
            if is_grating:
                filled, interpolated = smooth_background_spline(wave, median, grating_smoothing_um)
                displayed = np.where(np.isfinite(median), filled, np.nan)
            else:
                displayed = median
                filled, interpolated = interpolate_background_gaps(wave, displayed)
            label = 'Median bkg (traces masked)' if detect_traces else 'Median bkg'
            axis.plot(wave, displayed, color='black', lw=2.5, label=label)
            # Include the adjoining valid endpoints in each dashed segment.
            gap_indices = np.flatnonzero(interpolated)
            for index, run in enumerate(np.split(gap_indices, np.flatnonzero(np.diff(gap_indices) > 1) + 1)):
                if not run.size:
                    continue
                segment = slice(max(0, run[0] - 1), min(len(wave), run[-1] + 2))
                axis.plot(wave[segment], filled[segment], color='black', lw=2.5, ls='--',
                          label='Interpolated gap' if index == 0 else None)
            np.savetxt(OUT / ('_'.join(configuration) + '_median_bkg.csv'),
                       np.column_stack((wave, median)), delimiter=',',
                       header='wavelength_um,median_bkg_MJy_sr', comments='')
            spectrum_path = OUT / ('_'.join(configuration) + '_smoothed_bkg.csv')
            np.savetxt(spectrum_path, np.column_stack((wave, filled, interpolated.astype(int))),
                       delimiter=',', header='wavelength_um,median_bkg_MJy_sr,interpolated', comments='')
            spectrum_path.with_suffix('.json').write_text(json.dumps(dict(
                object_directory=ROOT.parent.name, product_target=target,
                grating=configuration[0], filter=configuration[1], units='MJy/sr',
                wavelength_units='micron', contributing_slits=[name for name, _, _ in spectra],
                smoothing='robust cubic spline' if is_grating else 'PRISM unchanged; PCHIP gap interpolation',
                grating_smoothing_um=grating_smoothing_um if is_grating else None,
                crowding_table=None if crowding_notes is None else str(Path(crowding_notes).expanduser().resolve()),
                interpolated_bins=int(interpolated.sum()), uncertainty='not estimated for fitted continuum'),
                indent=2) + '\n')
        else:
            # A rerun with no usable sky must not leave an old spectrum available
            # for subtraction under the current selection.
            stale = OUT / ('_'.join(configuration) + '_smoothed_bkg.csv')
            if stale.exists():
                stale.replace(stale.with_suffix('.csv.unavailable'))
            stale.with_suffix('.json').write_text(json.dumps(dict(status='unavailable',
                grating=configuration[0], filter=configuration[1], reason='No common good-slit spectrum'), indent=2) + '\n')
            axis.text(.5, .5, 'No common good-slit spectrum', ha='center', transform=axis.transAxes)
        if jbt_background:
            try:
                # All slit products in a configuration share the L3 observation metadata.
                source = next(filename for filename in all_candidates
                              if (fits.getheader(filename)['GRATING'], fits.getheader(filename)['FILTER']) == configuration)
                predicted_wave, predicted_total, provenance = predicted_jbt_background(
                    fits.getheader(source), ROOT / 'diagnostics' / 'jbt_cache', jbt_models)
                # Crop the model to plotted FS coverage so JBT's mid-IR tail
                # cannot expand the axes and hide the NIRSpec spectra.
                if spectra:
                    lo = min(x.min() for _, x, _ in spectra)
                    hi = max(x.max() for _, x, _ in spectra)
                    keep = (predicted_wave >= lo) & (predicted_wave <= hi)
                    axis.plot(predicted_wave[keep], predicted_total[keep], color='purple',
                              ls='--', lw=2, label='JBT predicted total')
                provenance.update(configuration=list(configuration), input=str(source))
                jbt_records.append(provenance)
            except Exception as error:
                print(f'JBT prediction unavailable for {configuration}: {error}')
                axis.text(.02, .02, 'JBT unavailable (see jbt_predictions.json)', fontsize=7,
                          transform=axis.transAxes)
                jbt_records.append(dict(configuration=list(configuration), status='unavailable', reason=str(error)))
        suffix = f'\nSmooth continuum: {grating_smoothing_um:.2f} μm scale' if is_grating else ''
        axis.set(title='/'.join(configuration) + ': median bkg spectrum' + suffix,
                 xlabel='Wavelength [μm]', ylabel='Median surface brightness [MJy/sr]')
        if spectra:
            axis.legend(fontsize=7)
    for configuration, group in axes_by_configuration.items():
        reference_values = orange_by_configuration.get(configuration, [])
        if not reference_values:
            continue
        values = np.concatenate(reference_values)
        if not values.size:
            continue
        lower, upper = float(values.min()), float(values.max())
        padding = .05 * (upper - lower) if upper > lower else max(abs(lower) * .05, 1e-6)
        for axis in group[1:]:
            axis.sharey(group[0])
        group[0].set_ylim(lower - padding, upper + padding)
    figure.suptitle(f'{ROOT.parent.name} — fixed-slit gradient diagnostics',
                    fontsize=18, fontweight='bold')
    figure.tight_layout(rect=(0, 0, 1, .95))
    figure.savefig(OUT / 'L2_L3_spatial_gradients.png', dpi=150)
    (OUT / 'gradient_summary.json').write_text(json.dumps(records, indent=2) + '\n')
    columns = ['target', 'object_directory', 'grating', 'filter', 'slit', 'human_mark', 'human_value',
               'excluded_from_background', 'non_flatness_percent', 'trace_masked', 'baseline_MJy_sr',
               'peak_departure_MJy_sr', 'peak_fraction', 'window_rows', 'status', 'reason',
               'trace_status', 'trace_detected', 'trace_sigma', 'trace_center_row', 'trace_fwhm_rows', 'file']
    with (OUT / 'spatial_flatness.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for record in records:
            metric = record.get('spatial_flatness', {})
            fraction = metric.get('fractional_departure')
            row = {key: record.get(key, '') for key in columns}
            row.update({key: metric.get(key, '') for key in
                        ['trace_masked', 'baseline_MJy_sr', 'peak_departure_MJy_sr', 'peak_fraction', 'window_rows']})
            row.update(target=target, object_directory=ROOT.parent.name,
                       human_value={'good': 1, 'bad': 0, 'whole_slit': -1}.get(record.get('human_mark'), ''),
                       non_flatness_percent='' if fraction is None else 100 * fraction,
                       status=metric.get('status', record.get('status', 'unavailable')),
                       reason=metric.get('reason', record.get('reason', '')))
            trace = record.get('trace', {})
            row.update(trace_status=trace.get('status', ''), trace_detected=trace.get('detected', ''),
                       trace_sigma=trace.get('significance', ''), trace_center_row=trace.get('center', ''),
                       trace_fwhm_rows=trace.get('fwhm', ''))
            writer.writerow(row)
    (OUT / 'background_candidate_profiles.json').write_text(json.dumps(comparison_records, indent=2) + '\n')
    (OUT / 'jbt_predictions.json').write_text(json.dumps(jbt_records, indent=2) + '\n')
    (OUT / 'diagnostic_settings.json').write_text(json.dumps(dict(
        crowding_table=None if crowding_notes is None else str(Path(crowding_notes).expanduser().resolve()),
        override_table=str(override_path.resolve()), grating_smoothing='robust cubic smoothing spline',
        grating_smoothing_um=grating_smoothing_um, smoothed_curve_used_when_FS_selected=True), indent=2) + '\n')
    figure.clear()
    print(f'Fixed-slit gradient diagnostics: {OUT / "L2_L3_spatial_gradients.png"}')
    return str(OUT / 'L2_L3_spatial_gradients.png')

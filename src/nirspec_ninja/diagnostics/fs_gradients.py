"""Reusable L2/L3 fixed-slit gradient diagnostics; no source detection."""
from pathlib import Path
import json
import warnings
import numpy as np
from scipy.ndimage import gaussian_filter1d
from astropy.io import fits
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg
from nirspec_ninja.subtract_background import fs_background as bg

# Diagnostic selection is independent of the slits allowed for sky subtraction.
# Left to right in the FOV, matching the finder charts.
GRADIENT_SLITS = ('S200A2', 'S400A1', 'S1600A1', 'S200A1')


def read_crowding_notes(path, target):
    """Read a whitespace table: target, slit columns containing 0/1, Notes.

    Target matching ignores case and periods (Cosmic.eye == COSMIC.EYE).
    Missing targets/slit columns remain unmarked; malformed tables fail loudly.
    """
    if path is None:
        return {}, ''
    normalize = lambda value: value.casefold().replace('.', '')
    header = None
    selected = None
    with Path(path).expanduser().open() as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split('#', 1)[0].split()
            if not fields:
                continue
            if header is None:
                header = fields
                if header[0].casefold() != 'target':
                    raise ValueError(f'{path}:{number}: expected target table header')
                columns = header[1:-1] if header[-1].casefold() == 'notes' else header[1:]
                if not columns or len(set(columns)) != len(columns) or any(
                        name not in GRADIENT_SLITS for name in columns):
                    raise ValueError(f'{path}:{number}: invalid slit columns')
                continue
            count = len(columns)
            if len(fields) < count + 1 or any(value not in ('0', '1') for value in fields[1:count + 1]):
                raise ValueError(f'{path}:{number}: each slit must be marked 0 or 1')
            if normalize(fields[0]) == normalize(target):
                if selected is not None:
                    raise ValueError(f'{path}:{number}: duplicate target {target}')
                selected = (dict(zip(columns, (value == '1' for value in fields[1:count + 1]))),
                            ' '.join(fields[count + 1:]))
    if header is None:
        raise ValueError(f'{path}: empty crowding table')
    return selected if selected is not None else ({}, '')

def profile(slit):
    good, excluded = bg.spatial_pixel_mask(slit)
    rows = np.flatnonzero(~excluded)
    if len(rows) < 4:
        raise ValueError('Too few usable spatial rows')
    lo = max(np.min(slit.wavelength[r][good[r]]) for r in rows)
    hi = min(np.max(slit.wavelength[r][good[r]]) for r in rows)
    if hi <= lo:
        raise ValueError('No common wavelength range across usable rows')
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        full = np.nanmedian(np.where(good, slit.brightness, np.nan), axis=1)[rows]
        common = good & (slit.wavelength >= lo) & (slit.wavelength <= hi)
        matched = np.nanmedian(np.where(common, slit.brightness, np.nan), axis=1)[rows]
    # Coordinate fractions are only a visual comparison, not a sky registration.
    fraction = (rows - rows[0]) / (rows[-1] - rows[0])
    return fraction, full, matched, (float(lo), float(hi))


def median_good_profiles(profiles, marks):
    """Median across good slits, then across normalized rows: a constant sky level."""
    good = [(name, x, values) for name, x, values in profiles if marks.get(name) is True]
    grid = np.linspace(0, 1, 101)
    if not good:
        return grid, None, []
    sampled = []
    for name, x, values in good:
        finite = np.isfinite(x) & np.isfinite(values)
        if finite.sum() >= 2:
            sampled.append((name, np.interp(grid, x[finite], values[finite], left=np.nan, right=np.nan)))
    if not sampled:
        return grid, None, []
    with warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        row_median = np.nanmedian([values for _, values in sampled], axis=0)
        median = np.full(grid.shape, np.nanmedian(row_median))
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


def spatial_median_spectrum(slit):
    """Align usable S2D rows in wavelength, then take their spatial median."""
    good, excluded = bg.spatial_pixel_mask(slit)
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


def smooth_grating_spectrum(values, fwhm_samples=25):
    """Display-only Gaussian smoothing, without joining masked spectral gaps."""
    result = np.full(values.shape, np.nan)
    indices = np.flatnonzero(np.isfinite(values))
    sigma = fwhm_samples / np.sqrt(8 * np.log(2))
    for run in np.split(indices, np.flatnonzero(np.diff(indices) > 1) + 1):
        if len(run):
            result[run] = gaussian_filter1d(values[run], sigma, mode='nearest')
    return result


def make_fs_gradient_report(fs_root, output_dir, target, crowding_notes=None):
    """Plot four calibrated fixed slits, optionally excluding human-marked bad slits."""
    marks, human_notes = read_crowding_notes(crowding_notes, target)
    eligible = tuple(name for name in GRADIENT_SLITS if marks.get(name, True))
    ROOT, OUT = Path(fs_root), Path(output_dir)
    OUT.mkdir(parents=True, exist_ok=True)
    l2 = {}
    for filename in sorted(ROOT.glob('L2b/*_s2d.fits')):
        header = fits.getheader(filename)
        for slit in bg.read_slit_pixels(filename):
            if slit.name in eligible:
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
    records = [dict(file=str(p), status='excluded', reason='Human-marked bad slit',
                    human_mark='bad', human_notes=human_notes)
               for p in all_candidates if any('_' + n.lower() + '_' in p.name
                                               for n in GRADIENT_SLITS if not marks.get(n, True))]
    candidates = [p for p in all_candidates if any('_' + n.lower() + '_' in p.name for n in eligible)]
    if not all_candidates:
        (OUT / 'gradient_summary.json').write_text(json.dumps(records, indent=2) + '\n')
        print(f'No eligible L3 S2Ds for {target}: no gradient plot produced')
        return None
    # Read bad L3 profiles too, solely for the de-emphasized comparison overlay.
    comparisons = {}
    spectra_by_configuration = {}
    for filename in all_candidates:
        header = fits.getheader(filename)
        slit = next(bg.read_slit_pixels(filename))
        if not slit.calibrated:
            continue
        try:
            x, _, values, _ = profile(slit)
        except ValueError:
            continue
        comparisons.setdefault((header['GRATING'], header['FILTER']), []).append((slit.name, x, values))
        if marks.get(slit.name) is True:
            try:
                wave, spectrum = spatial_median_spectrum(slit)
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
                                 end_difference=float(np.median(values[x >= .75]) - np.median(values[x <= .25]))))
        axis.plot(fraction, full, '--', label='L3: all valid wavelengths')
        axis.plot(fraction, matched, label='L3: common wavelength range')
        human_mark = 'good' if marks.get(slit.name) is True else 'unmarked'
        title = '/'.join(key) + (' [human: good]' if human_mark == 'good' else '')
        axis.set(title=title, xlabel='Fraction along usable rows (not sky-aligned)',
                 ylabel='Median surface brightness [MJy/sr]')
        axis.legend(fontsize=7)
        axes_by_configuration.setdefault(key[:2], []).append(axis)
        orange_by_configuration.setdefault(key[:2], []).append(matched[np.isfinite(matched)])
        records.append(dict(file=str(filename), range_um=interval,
                            human_mark=human_mark, human_notes=human_notes,
                            l3_end_difference=float(np.median(matched[fraction >= .75]) - np.median(matched[fraction <= .25])),
                            l2=upstream))
    comparison_records = []
    for configuration, profiles in sorted(comparisons.items()):
        axis = comparison_axes[configuration]
        profiles = sorted(profiles, key=lambda item: GRADIENT_SLITS.index(item[0]))
        grid, median, good_names = median_good_profiles(profiles, marks)
        for name, x, values in profiles:
            if marks.get(name) is True:
                axis.plot(x, values, lw=1.2, label=f'{name}: good')
            elif marks.get(name) is False:
                axis.plot(x, values, color='crimson', alpha=.25, lw=.8, ls='--',
                          label=f'{name}: bad', zorder=1)
        if median is not None:
            axis.plot(grid, median, color='black', lw=2.5, label='Sky: median over good slits and rows', zorder=5)
            # Only explicitly good slits determine the comparison stretch.
            good_values = np.concatenate([values[np.isfinite(values)] for name, _, values in profiles
                                          if name in good_names])
            lower, upper = float(good_values.min()), float(good_values.max())
            padding = .05 * (upper - lower) if upper > lower else max(abs(lower) * .05, 1e-6)
            axis.set_ylim(lower - padding, upper + padding)
        else:
            axis.text(.5, .5, 'No usable human-marked good slits', ha='center', transform=axis.transAxes)
        axis.set(title='/'.join(configuration) + ': median bkg',
                 xlabel='Fraction along usable rows (not sky-aligned)',
                 ylabel='Median surface brightness [MJy/sr]')
        handles, _ = axis.get_legend_handles_labels()
        if handles:
            axis.legend(fontsize=7)
        comparison_records.append(dict(configuration=list(configuration), good_slits=good_names,
                                       sky_level=None if median is None else float(median[0]),
                                       fraction=grid.tolist(),
                                       median_good=None if median is None else
                                       [float(v) if np.isfinite(v) else None for v in median],
                                       note='Median across good slits, then across normalized spatial rows. '
                                            'Equal slit weights; interpolated normalized rows, not sky registration. '
                                            'Each slit uses its own common wavelength interval.'))
    for row, configuration in enumerate(configurations):
        axis = axes[row, -1]
        spectra = sorted(spectra_by_configuration.get(configuration, []),
                         key=lambda item: GRADIENT_SLITS.index(item[0]))
        is_grating = configuration[0].upper() != 'PRISM'
        for name, wave, values in spectra:
            displayed = smooth_grating_spectrum(values) if is_grating else values
            axis.plot(wave, displayed, lw=1, label=name)
        wave, median = combined_median_spectrum(spectra)
        if len(wave) and np.any(np.isfinite(median)):
            displayed = smooth_grating_spectrum(median) if is_grating else median
            axis.plot(wave, displayed, color='black', lw=2.5, label='Combined median bkg')
            np.savetxt(OUT / ('_'.join(configuration) + '_median_bkg.csv'),
                       np.column_stack((wave, median)), delimiter=',',
                       header='wavelength_um,median_bkg_MJy_sr', comments='')
        else:
            axis.text(.5, .5, 'No common good-slit spectrum', ha='center', transform=axis.transAxes)
        suffix = '\nGaussian FWHM: 25 wavelength samples' if is_grating else ''
        axis.set(title='/'.join(configuration) + ': median bkg spectrum' + suffix,
                 xlabel='Wavelength [μm]', ylabel='Median surface brightness [MJy/sr]')
        if spectra:
            axis.legend(fontsize=7)
    for configuration, group in axes_by_configuration.items():
        values = np.concatenate(orange_by_configuration[configuration])
        if not values.size:
            continue
        lower, upper = float(values.min()), float(values.max())
        padding = .05 * (upper - lower) if upper > lower else max(abs(lower) * .05, 1e-6)
        for axis in group[1:]:
            axis.sharey(group[0])
        group[0].set_ylim(lower - padding, upper + padding)
    figure.tight_layout()
    figure.savefig(OUT / 'L2_L3_spatial_gradients.png', dpi=150)
    (OUT / 'gradient_summary.json').write_text(json.dumps(records, indent=2) + '\n')
    (OUT / 'background_candidate_profiles.json').write_text(json.dumps(comparison_records, indent=2) + '\n')
    figure.clear()
    print(f'Fixed-slit gradient diagnostics: {OUT / "L2_L3_spatial_gradients.png"}')
    return str(OUT / 'L2_L3_spatial_gradients.png')

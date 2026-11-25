"""Fixed-slit median extraction, moved from the supplied jrr/jrjwst.py.

Original function names and extraction calculations are preserved.
Source detection uses only the migrated source_detection module.
"""
import glob
import warnings
from os.path import basename
from pathlib import Path
from re import sub
import numpy as np
import pandas
from astropy.io import fits
from ._extraction_helpers import gethead, find_science_extensions, get_nirspec_good_fixedslit_names, get_keywords_from_jwst_header2, rebin_spec_new
from .spatial_diagnostics import diagnostic_trace, trace_row_mask, profile, spatial_flatness
from nirspec_ninja.subtract_background.fs_background import SlitPixels, spatial_pixel_mask


def measure_extraction_diagnostics(infile, thisslit, intype, median_2D,
                                   median_wave_1D, sci_extension, clip_nan=5,
                                   detect_traces=True):
    """Measure CSV diagnostics independently of the 1D extraction calculation."""
    use_error = False
    if intype == 'S2D':
        with fits.open(infile, memmap=False) as hdus:
            sci = hdus[sci_extension]
            version = sci.header.get('EXTVER', 1)
            data = np.array(sci.data, dtype=float, copy=True)
            wave = np.array(hdus['WAVELENGTH', version].data, dtype=float, copy=True)
            use_error = ('ERR', version) in hdus
            error = np.array(hdus['ERR', version].data, dtype=float, copy=True) if use_error else np.ones_like(data)
            dq = np.array(hdus['DQ', version].data, copy=True) if ('DQ', version) in hdus else None
    else:
        data = np.array(median_2D, dtype=float, copy=True)
        wave = np.broadcast_to(median_wave_1D, data.shape).copy()
        error, dq = np.ones_like(data), None
    slit = SlitPixels(thisslit, data, wave, error, dq)
    trace = diagnostic_trace(slit, clip_nan=clip_nan, use_error=use_error) if detect_traces else dict(
        status='skipped', reason='trace detection disabled', detected=False)
    flatness = dict(status='unavailable', reason='No usable spatial profile')
    try:
        x, _, values, _ = profile(slit, clip_nan=clip_nan)
        _, excluded = spatial_pixel_mask(slit, clip_nan=clip_nan)
        rows = np.flatnonzero(~excluded)
        mask = trace_row_mask(trace, data.shape[0])
        values[mask[rows[0]:rows[-1] + 1]] = np.nan
        flatness = spatial_flatness(x, values)
        flatness['trace_masked'] = bool(mask.any())
        if flatness['status'] == 'ok' and flatness.get('fractional_departure') is None:
            flatness.update(status='unavailable', reason='Near-zero background baseline')
    except ValueError as error:
        flatness['reason'] = str(error)
    return trace, flatness

def find_science_extensions_nirspecFS(infile):
    # JWST level 3 calfiles have extensions for each exposure.  Group the science extensions
    # w their fixed slits

    with fits.open(infile) as fitsfile:
        N_extensions = len(fitsfile)

    out_dict = {} 
    fixed_slit_names = get_nirspec_good_fixedslit_names()
    for thisslit in fixed_slit_names:
        out_dict[thisslit] = []
        sci_extensions = find_science_extensions(infile, kind='SCI')
        for ii in sci_extensions:
            sltname = gethead(infile, 'SLTNAME', extension=ii)
            if sltname == thisslit:
                out_dict[thisslit].append(ii)
    return(out_dict)

def _expand_row_mask(anchor_rows, radius, row_count):
    """Return a mask covering each anchor row and all rows within a radius."""
    expanded_row_mask = np.zeros(row_count, dtype=bool)
    for anchor_row in anchor_rows:
        start = max(0, anchor_row - radius)
        # The upper slice bound is exclusive, hence the extra +1.
        stop = min(row_count, anchor_row + radius + 1)
        expanded_row_mask[start:stop] = True
    return expanded_row_mask

def _nirspecfs_nan_mask_rows(sci_image, clip_nan):
    nan_counts = np.count_nonzero(np.isnan(sci_image), axis=1)
    # Integer arithmetic implements "at least half" without a fractional
    # threshold when the wavelength-axis length is odd.
    mostly_nan_rows = np.flatnonzero(2 * nan_counts >= sci_image.shape[1])
    # Edge rows are always anchors, even when their pixels are all finite.
    masking_anchor_rows = np.unique(np.concatenate((
        mostly_nan_rows, np.array([0, sci_image.shape[0] - 1]))))
    excluded_row_mask = _expand_row_mask(
        masking_anchor_rows, clip_nan, sci_image.shape[0])
    neighbor_masked_rows = np.flatnonzero(
        excluded_row_mask &
        ~np.isin(np.arange(sci_image.shape[0]), mostly_nan_rows))
    return mostly_nan_rows, neighbor_masked_rows, excluded_row_mask

def median_combine_level3_nirspecFS2(infile, thisslit, outdir, sci_to_wave_off='Default', \
                write2D=False, clip_nan=5, debug=False, detect_traces=True):
    """  Makes a spatial (and exposure level, if CAL file) median combine of a level 3 JWST
    input spectral file.  Should handle L3 CAL or S2D files.  For making background measurements.
    v2.0.1 of JWST pipeline now makes a level 3 calfile that has an extension for every exposure.
    In principle, using this file may be best, since it ignores WCS.   In practice, outputs are v similar.
    sci_to_wave_off is the offset of extension nubmers of SCI and WAVELENGTH. Changed from edge-clipping
    in v1 to clipping within clip_nan rows of the edges OR rows that are mostly nans. """
    if not isinstance(clip_nan, (int, np.integer)) or clip_nan < 0:
        raise ValueError('clip_nan must be a non-negative integer')

    if ('s2d' in infile) or ('S2D' in infile):
        intype = 'S2D'
    elif ('cal' in infile) or ('CAL' in infile):
        intype = 'CAL'
    else:
        raise ValueError('Could not determine whether input is S2D or CAL: ' + infile)

    offsets = {'CAL': 3, 'S2D': 2}
    if sci_to_wave_off == 'Default':
        sci_to_wave_off = offsets[intype]
    fixed_slit_names = get_nirspec_good_fixedslit_names()
    if thisslit not in fixed_slit_names:
        raise Exception('ERROR: fixed slit name', thisslit, 'not in', fixed_slit_names)
    sci_extensions = find_science_extensions_nirspecFS(infile)
    if len(sci_extensions[thisslit]) == 0:
        raise Exception('ERROR: no SCI extensions for this slit', thisslit, infile)

    sci_images = []
    wave_images = []
    for sci_extension in sci_extensions[thisslit]:
        with fits.open(infile, mode='readonly', memmap=False) as sfile:
            sci_image = np.array(
                sfile[sci_extension].data, dtype=float, copy=True)
            wave_image = np.array(
                sfile[sci_extension + sci_to_wave_off].data,
                dtype=float, copy=True)

        mostly_nan_rows, neighbor_masked_rows, excluded_row_mask = \
            _nirspecfs_nan_mask_rows(sci_image, clip_nan)

        if debug:
            print(
                f'SCI extension {sci_extension} mostly-NaN rows: '
                f'{mostly_nan_rows.tolist()}')
            print(
                f'SCI extension {sci_extension} other rows masked within '
                f'{clip_nan} rows '
                f'of a mostly-NaN row or image edge: '
                f'{neighbor_masked_rows.tolist()}')

        sci_image[excluded_row_mask, :] = np.nan
        wave_image[excluded_row_mask, :] = np.nan
        sci_images.append(sci_image)
        wave_images.append(wave_image)

    with warnings.catch_warnings():
        warnings.filterwarnings(
            'ignore',
            message=r'All-NaN slice encountered',
            category=RuntimeWarning,
        )

        image_shapes = [image.shape for image in wave_images]

        if len(set(image_shapes)) == 1:
            # Preserve the existing behavior when every exposure has the
            # same dimensions.
            wavelength_cube = np.stack(wave_images)
            median_wave_1D = np.nanmedian(
                wavelength_cube, axis=(0, 1)
            )
        else:
            # Some CAL products have exposure cutouts that differ by one
            # spatial or wavelength pixel. Construct a common wavelength
            # grid without first forcing the images into a rectangular cube.
            wavelength_length = min(
                image.shape[1] for image in wave_images
            )
            common_pixel_coordinate = np.linspace(
                0.0, 1.0, wavelength_length
            )

            normalized_wave_rows = []

            for wave_image in wave_images:
                native_pixel_coordinate = np.linspace(
                    0.0, 1.0, wave_image.shape[1]
                )

                for wave_row in wave_image:
                    finite = np.isfinite(wave_row)

                    if np.count_nonzero(finite) < 2:
                        continue

                    normalized_wave_rows.append(
                        np.interp(
                            common_pixel_coordinate,
                            native_pixel_coordinate[finite],
                            wave_row[finite],
                            left=np.nan,
                            right=np.nan,
                        )
                    )

            if not normalized_wave_rows:
                raise ValueError(
                    'No usable wavelength rows found in ' + infile
                )

            median_wave_1D = np.nanmedian(
                np.stack(normalized_wave_rows), axis=0
            )

        # Allocate a rectangular output only after every exposure has been
        # placed onto the common wavelength grid. Spatially shorter images
        # are padded with NaNs.
        maximum_spatial_rows = max(
            image.shape[0] for image in sci_images
        )

        resampled_science_cube = np.full(
            (
                len(sci_images),
                maximum_spatial_rows,
                median_wave_1D.size,
            ),
            np.nan,
            dtype=float,
        )

        for exposure_index, (science_image, wavelength_image) in enumerate(
                zip(sci_images, wave_images)):

            # Approximately center cutouts whose spatial dimensions differ.
            row_offset = (
                maximum_spatial_rows - science_image.shape[0]
            ) // 2

            for spatial_row in range(science_image.shape[0]):
                science_row = science_image[spatial_row, :]
                wavelength_row = wavelength_image[spatial_row, :]

                if (
                    not np.any(np.isfinite(science_row))
                    or not np.any(np.isfinite(wavelength_row))
                ):
                    continue

                resampled_science_cube[
                    exposure_index,
                    row_offset + spatial_row,
                    :,
                ] = rebin_spec_new(
                    wavelength_row,
                    science_row,
                    median_wave_1D,
                )

        median_2D = np.nanmedian(
            resampled_science_cube, axis=0
        )

        if write2D:
            fits.writeto(
                outdir + 'FS_2Dmedian_from'
                + intype + '_' + thisslit + '.fits',
                median_2D,
                overwrite=True,
            )

        median_sci_1D = np.nanmedian(
            resampled_science_cube, axis=(0, 1)
        )

        with warnings.catch_warnings():
            warnings.filterwarnings(
                'ignore',
                message=r'Degrees of freedom <= 0 for slice\.',
                category=RuntimeWarning,
            )
            std_1D = np.nanstd(
                resampled_science_cube, axis=(0, 1)
            )

        mad = np.nanmedian(
            np.absolute(
                resampled_science_cube - median_sci_1D
            ),
            axis=(0, 1),
        )

    df = pandas.DataFrame({
        'wave': median_wave_1D, 'fnu': median_sci_1D,
        'stddev': std_1D, 'mad': mad})
    fileroot = infile.split('/')[-1].split('_s')[0]
    outfile    = sub('.fits', '_extract.csv', basename(infile))

    # Prep the header     
    myheaderdict =  get_keywords_from_jwst_header2(infile)
    try:
        trace, flatness = measure_extraction_diagnostics(
            infile, thisslit, intype, median_2D, median_wave_1D,
            sci_extensions[thisslit][0], clip_nan=clip_nan, detect_traces=detect_traces)
    except Exception as error:
        reason = f'{type(error).__name__}: {error}'
        trace = dict(status='failed', reason=reason, detected=False)
        flatness = dict(status='unavailable', reason=reason)
    if trace['status'] == 'ok' and trace.get('significance') is not None:
        myheaderdict['trace_sig2'] = trace['significance']
    fraction = flatness.get('fractional_departure')
    myheaderdict['non-flat'] = 100 * fraction if fraction is not None else np.nan
    myheaderdict['non-flat_status'] = flatness['status']
    if flatness['status'] != 'ok':
        myheaderdict['non-flat_reason'] = flatness.get('reason', 'unavailable')
    header = '# Custom background from NIRSpec fixed slit, from median combine of all exposures and\n'
    header += '# spatial direction, from Level 3 CAL or S2D file\n'
    header += '# filename ' + basename(infile) + '\n'
    header += '# Diagnostics: trace_sig2 is corrected sigma; non-flat is spatial non-flatness (%)\n'
    for key, value in myheaderdict.items():
        header += '##' + key + ' ' + str(value) + '\n'
    with open(outdir + outfile, "w") as handle:
        handle.write(header)    # write the header
        df.to_csv(handle, index=False) # write the file
    df.attrs['extraction_report'] = {
        'input': str(Path(infile).resolve()), 'slit': thisslit,
        'output': str(Path(outdir + outfile).resolve()),
        'status': 'ok', 'trace': trace, 'spatial_flatness': flatness,
    }
    if trace['status'] == 'failed':
        print(f"TRACE FAILED | {Path(infile).resolve()} | slit={thisslit} | "
              f"{trace['reason']} | CSV written; trace_sig2 unavailable", flush=True)
    return df

def wrap_median_combine_level3_nirspecFS(indir, outdir, clip_nan=5, return_report=False, detect_traces=True):
    """Return spectra, or compact per-file reports while continuing after errors.

    The default preserves the historical dictionary of DataFrames and raises
    extraction errors. Reporting mode records errors and tries the next file.
    """
    df = {}
    reports = []
    for kind in ('s2d', 'cal'):
        infiles = sorted(glob.glob(str(Path(indir) / ('*' + kind + '*.fits'))))
        for infile in infiles:
            label = basename(infile)
            thisslit = None
            try:
                thisslit = gethead(infile, 'SLTNAME')
                if thisslit == 'S200B1':
                    continue  # ignore this uncalibrated slit
                frame = median_combine_level3_nirspecFS2(
                    infile, thisslit, outdir, clip_nan=clip_nan, detect_traces=detect_traces)
                if return_report:
                    reports.append(frame.attrs['extraction_report'])
                else:
                    df[label] = frame
            except Exception as error:
                if not return_report:
                    raise
                reason = f'{type(error).__name__}: {error}'
                reports.append({'input': str(Path(infile).resolve()),
                                'slit': thisslit, 'status': 'failed',
                                'reason': reason, 'output': None})
                print(f'EXTRACTION FAILED | {Path(infile).resolve()} | '
                      f'slit={thisslit} | {reason}', flush=True)
    return reports if return_report else df

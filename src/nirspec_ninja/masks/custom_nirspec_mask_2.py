"""NIRSpec WCS masks with optional RATE refinement for JWST 3.0 flicker cleaning.

Default grating geometry uses the pipeline's nominal science WCS footprints.
Optional mask_grating_extrapolate continues full apertures to the red detector edge,
with an optical IFU fallback when no nominal footprint falls on the detector.
These extensions are coverage estimates, not order identifications.
See docs/masks/HIGHER_ORDER_MASK_CONTEXT.md in the package repository for evidence and next steps.
"""
import os
import json
import hashlib
import tempfile
from contextlib import ExitStack
from pathlib import Path

# STScI helpdesk says these os commands need to come BEFORE jwst pipeline packages are imported
home = "/Users/jrrigby1/Ref_files/"
os.environ.setdefault("CRDS_PATH", home + "crds_cache/jwst_ops")
os.environ.setdefault("CRDS_SERVER_URL", "https://jwst-crds.stsci.edu")
mycontext = 'jwst_1584.pmap'
os.environ.setdefault("CRDS_CONTEXT", mycontext)   # v3.0.0

import numpy as np
from astropy.io import fits
import gwcs
from gwcs.utils import to_index
from jwst import datamodels
from jwst.assign_wcs import AssignWcsStep, nirspec
from jwst.assign_wcs.util import NoDataOnDetectorError
from jwst.msaflagopen import MSAFlagOpenStep
from stdatamodels.jwst.datamodels import dqflags, read_metadata
from .aperture_extrapolation import extend_aperture


def refine_mask_from_rate(background_mask, rate, nsigma=3.0, pixel_dq=None):
    """Intersect a custom mask with global rejection on a calibrated RATE.

    Arrays are never modified in place. Supply the uncleaned combined 2D rate;
    this performs no ramp fitting, WCS assignment, or pipeline orchestration.
    Optional calibrated pixel_dq excludes DO_NOT_USE, NONLINEAR, NO_LIN_CORR,
    and NO_SAT_CHECK. SATURATED remains a per-group cleaning exclusion.
    """
    from jwst.clean_flicker_noise.clean_flicker_noise import clip_to_background
    if isinstance(nsigma, (bool, np.bool_)) or not np.isfinite(nsigma) or nsigma <= 0:
        raise ValueError('nsigma must be finite and positive')
    data = np.asarray(rate)
    mask = np.asarray(background_mask, dtype=bool).copy()
    if data.ndim != 2 or mask.shape != data.shape:
        raise ValueError('Custom mask and calibrated RATE must have matching 2D shapes')
    if pixel_dq is not None:
        dq = np.asarray(pixel_dq)
        if dq.shape != data.shape or not np.issubdtype(dq.dtype, np.integer):
            raise ValueError('Calibrated pixel DQ must be an integer array matching RATE')
        excluded = (dqflags.pixel['DO_NOT_USE'] | dqflags.pixel['NONLINEAR']
                    | dqflags.pixel['NO_LIN_CORR'] | dqflags.pixel['NO_SAT_CHECK'])
        mask &= (dq & excluded) == 0
    mask &= np.isfinite(data) & (data != 0)
    clip_to_background(data.copy(), mask, sigma_upper=nsigma)
    if not np.any(mask):
        raise ValueError('No usable background remains after custom-mask rejection')
    return mask


def refine_mask_file_from_rate(mask_file, rate, nsigma=3.0, pixel_dq=None):
    """Atomically retain the globally refined exposure mask; return its array."""
    target = Path(mask_file).expanduser().resolve()
    with datamodels.open(target) as model:
        mask = refine_mask_from_rate(model.data, rate, nsigma=nsigma, pixel_dq=pixel_dq)
        model.data = mask.astype(np.float32)
        with tempfile.NamedTemporaryFile(dir=target.parent, suffix='.fits', delete=False) as handle:
            temporary = Path(handle.name)
        try:
            model.save(str(temporary))
            os.replace(temporary, target)
        finally:
            temporary.unlink(missing_ok=True)
    return mask


def _ifu_cache_path(image, directory, shape, options):
    """Resolve references before lookup; nod/pointing/velocity are not cache keys.

    IFU footprints can differ slightly at a wavelength boundary across dithers.
    This cache deliberately shares the first geometry within an instrument
    configuration; per-exposure illumination refinement remains independent.
    """
    import jwst
    references = {}
    modes = ['NRS_IFU'] + (['NRS_FIXEDSLIT'] if options['fs_wcs'] else [])
    with image.copy() as model:
        for mode in modes:
            model.meta.exposure.type = mode
            steps = [AssignWcsStep()]
            if mode == 'NRS_IFU' and options['failed_open']:
                steps.append(MSAFlagOpenStep())
            for step in steps:
                for kind in step.reference_file_types:
                    reference = step.get_reference_file(model, kind)
                    identity = [reference]
                    if reference and reference != 'N/A':
                        stat = Path(reference).stat()
                        identity += [stat.st_size, stat.st_mtime_ns]
                    references[f'{mode}:{kind}'] = identity
    payload = dict(
        version=jwst.__version__, context=os.environ.get('CRDS_CONTEXT'),
        code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        instrument=image.meta.instrument.instance,
        subarray=image.meta.subarray.instance, shape=shape,
        options=options, references=references,
    )
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, default=str).encode()).hexdigest()
    cache_dir = Path(directory).expanduser().resolve()
    cache_dir.mkdir(parents=True, exist_ok=True)
    return cache_dir / (digest + '.npz')


def _add_trace_search_region(region, footprint, xstart, ystart, dispersion, cross):
    """Freeze a narrow search corridor around one detector-x-dispersed trace.

    Preserve the pixel footprint and shutter gaps. Extend each wavelength
    endpoint along the centerline slope fitted to its last 20 occupied columns.
    Endpoint path length is bounded by dispersion; vertical adjustments are
    bounded by cross detector pixels (conservative for nearly horizontal traces).
    No newly detected illumination or DQ exclusion can seed further growth.
    """
    from scipy.ndimage import maximum_filter1d

    columns = np.flatnonzero(footprint.any(axis=0))
    if not len(columns):
        return
    height, width = region.shape
    # Padding only defines eligibility; it never directly excludes pixels.
    padded = np.pad(footprint, ((cross, cross), (0, 0)))
    padded = maximum_filter1d(padded, size=2 * cross + 1, axis=0,
                             mode='constant', cval=0)
    lo, hi = max(0, ystart - cross), min(height, ystart + footprint.shape[0] + cross)
    region[lo:hi, xstart:xstart + footprint.shape[1]] |= padded[
        lo - (ystart - cross):hi - (ystart - cross)]
    for direction, end, sample in ((-1, columns[0], columns[:20]), (1, columns[-1], columns[-20:])):
        centers = np.array([np.median(np.flatnonzero(footprint[:, col])) for col in sample])
        slope = np.polyfit(sample - end, centers, 1)[0] if len(sample) > 1 else 0.0
        rows = np.flatnonzero(footprint[:, end]) + ystart
        for step in range(1, int(np.floor(dispersion / np.hypot(1, slope))) + 1):
            x = xstart + end + direction * step
            if not 0 <= x < width:
                continue
            shift = int(np.rint(slope * direction * step))
            if np.hypot(step, shift) > dispersion:
                continue
            for offset in range(-cross, cross + 1):
                y = rows + shift + offset
                region[y[(y >= 0) & (y < height)], x] = True


def _mask_slits_by_wcs(model, mask, trace_region=None, limits=(30, 3), extrapolation=None):
    """Exclude slit footprints without MOS extraction margins (JWST 3.0)."""
    height, width = mask.shape
    slits = model.meta.wcs.get_transform('gwa', 'slit_frame').slits
    boxes = []
    for slit in slits:
        wcs = nirspec.nrs_wcs_set_input(model, slit.name)
        bbox = wcs.bounding_box
        boxes.append(bbox)
        # Evaluate integer detector pixel centers, including both bbox ends,
        # and clip before indexing so subarrays cannot wrap around the mask.
        xstart = max(0, int(np.ceil(bbox[0][0])))
        xstop = min(width, int(np.floor(bbox[0][1])) + 1)
        ystart = max(0, int(np.ceil(bbox[1][0])))
        ystop = min(height, int(np.floor(bbox[1][1])) + 1)
        if xstart >= xstop or ystart >= ystop:
            continue
        y, x = np.mgrid[ystart:ystop, xstart:xstop]
        _, slit_y, wavelength = wcs.get_transform('detector', 'slit_frame')(x, y)
        # get_open_msa_slits adds 0.5 at each outer end for extraction.
        # Fixed slits (quadrant 5) do not receive that MOS margin.
        margin = 0.0 if slit.quadrant == 5 else 0.5
        illuminated = (np.isfinite(wavelength) & np.isfinite(slit_y)
                       & (slit_y >= slit.ymin + margin)
                       & (slit_y <= slit.ymax - margin))
        mask[ystart:ystop, xstart:xstop][illuminated] = False
        if trace_region is not None:
            _add_trace_search_region(trace_region, illuminated, xstart, ystart, *limits)
        if extrapolation is not None:
            report = extend_aperture(mask, trace_region, illuminated, xstart, ystart,
                                     search_padding=limits[1], label=f'slit {slit.name}')
            if report is not None:
                extrapolation.append(report)
    return boxes


def _extend_mask_from_rate(mask, rate, valid, boxes, padding, nsigma, trace_region):
    """Add measured illumination only within frozen trace corridors."""
    from scipy.ndimage import median_filter

    if trace_region.shape != mask.shape or trace_region.dtype != np.dtype(bool):
        raise ValueError('Trace search region must be boolean and match the science grid')
    samples = mask & valid
    if np.count_nonzero(samples) < 100:
        raise ValueError('Too few valid background pixels in illumination reference')
    baseline = np.median(rate[samples])
    # Suppress isolated hot pixels and noise while retaining narrow sky bands.
    smooth = median_filter(np.where(valid, rate, baseline), size=(3, 5), mode='nearest')
    values = smooth[samples]
    background = np.median(values)
    sigma = 1.4826 * np.median(np.abs(values - background))
    if not np.isfinite(sigma) or sigma <= 0:
        raise ValueError('Cannot estimate background scatter from illumination reference')
    candidates = np.zeros(mask.shape, dtype=bool) if padding is not None else trace_region.copy()
    height, width = mask.shape
    padx, pady = (0, 0) if padding is None else padding
    for bbox in (() if padding is None else boxes):
        xstart = max(0, int(np.ceil(bbox[0][0])) - padx)
        xstop = min(width, int(np.floor(bbox[0][1])) + 1 + padx)
        ystart = max(0, int(np.ceil(bbox[1][0])) - pady)
        ystop = min(height, int(np.floor(bbox[1][1])) + 1 + pady)
        if xstart < xstop and ystart < ystop:
            candidates[ystart:ystop, xstart:xstop] = True
    candidates &= trace_region
    mask[candidates & ((smooth > background + nsigma * sigma) | ~valid)] = False


def _ramp_illumination_proxy(input_data, background_mask, meta):
    """Sparse raw group differences in DN/s; no Detector1 or ramp fitting."""
    import warnings

    def median(values, **kwargs):
        # All-invalid regions remain NaN and are excluded by the caller.
        with warnings.catch_warnings():
            warnings.filterwarnings('ignore', message='All-NaN slice encountered', category=RuntimeWarning)
            return np.nanmedian(values, **kwargs)

    with ExitStack() as stack:
        scale, offset, integer_raw = 1.0, 0.0, False
        if isinstance(input_data, (str, Path)):
            hdus = stack.enter_context(fits.open(input_data, memmap=True, do_not_scale_image_data=True))
            data = hdus['SCI'].data
            scale = hdus['SCI'].header.get('BSCALE', 1)
            offset = hdus['SCI'].header.get('BZERO', 0)
            integer_raw = hdus['SCI'].header['BITPIX'] == 16
        else:
            data = input_data.data
            integer_raw = data.dtype.kind == 'u' and data.dtype.itemsize == 2
        if data.ndim != 4 or data.shape[1] < 3:
            raise ValueError('Ramp illumination requires at least three groups')
        nints, ngroups, ny, nx = data.shape
        group_time = meta.exposure.group_time
        if group_time is None or not np.isfinite(group_time) or group_time <= 0:
            raise ValueError('Ramp illumination requires a finite positive TGROUP')
        slowaxis = abs(meta.subarray.slowaxis or 0) - 1
        if slowaxis not in (0, 1):
            raise ValueError('Ramp illumination requires SLOWAXIS = +/-1 or +/-2')
        rows = np.ones(ny, dtype=bool)
        if (ny, nx) == (3200, 2048) and 'IRS2' in (meta.exposure.readpatt or '').upper():
            from jwst.refpix.irs2_subtract_reference import make_irs2_mask
            rows = make_irs2_mask(3200, 2048, meta.exposure.nrs_normal or 16,
                                  meta.exposure.nrs_reference or 4)
            if meta.instrument.detector == 'NRS2':
                rows = rows[::-1]
        if (int(rows.sum()), nx) != background_mask.shape:
            raise ValueError('Ramp illumination science grid does not match mask')
        integrations = np.unique(np.linspace(0, nints - 1, min(nints, 4), dtype=int))
        middle = (ngroups + 1) // 2
        pairs = [(1, middle), (middle, ngroups - 1)] if ngroups >= 7 else [(1, ngroups - 1)]
        samples = []
        for integration in integrations:
            for first, last in pairs:
                a = np.asarray(data[integration, first, rows, :], dtype=np.float32)
                b = np.asarray(data[integration, last, rows, :], dtype=np.float32)
                valid = np.isfinite(a) & np.isfinite(b)
                if integer_raw:
                    # Reject raw ADC clipping; this is not calibrated saturation flagging.
                    valid &= (a * scale + offset < 65530) & (b * scale + offset < 65530)
                difference = (b - a) * scale / ((last - first) * group_time)
                difference[~valid] = np.nan
                background = valid & background_mask
                count = difference.shape[slowaxis]
                channel_size = 512 if count == 2048 else count
                for start in range(0, count, channel_size):
                    region = ((slice(start, start + channel_size), slice(None)) if slowaxis == 0
                              else (slice(None), slice(start, start + channel_size)))
                    block = difference[region]
                    baseline = median(np.where(background[region], block, np.nan),
                                      axis=slowaxis, keepdims=True)
                    block -= np.nan_to_num(baseline)
                baseline = median(np.where(background, difference, np.nan),
                                  axis=1 - slowaxis, keepdims=True)
                difference -= np.nan_to_num(baseline)
                samples.append(difference)
        return median(np.stack(samples), axis=0)


def make_custom_NIRSpec_mask(
    input_data,
    output_file,
    *,
    exp_type_override=None,
    mask_IFU=None,
    mask_MSA=None,
    mask_MSA_by_wcs=False,
    mask_FS_by_wcs=False,
    mask_grating_extrapolate=None,
    mask_by_illumination=False,
    mask_illumination_reference=None,
    mask_illumination_nsigma=5.0,
    mask_illumination_padding=None,
    mask_illumination_dispersion_pixels=30,
    mask_illumination_cross_dispersion_pixels=3,
    mask_FS_strip=None,
    mask_FSprism=False,
    mask_failed_open=None,
    mask_box=False,
    mask_corners=None,
    mask_origin=1,
    single_mask=True,
    input_dir=None,
    overwrite=False,
    ifu_geometry_cache_dir=None,
    wcs_reference_overrides=None,
):
    """Write user_mask FITS; return True=background, False=exclude.

    Accept UNCAL/RATE/RATEINTS FITS paths or unextracted NIRSpec datamodels.
    FITS inputs are read for headers/metadata (and existing failed-open DQ
    flags), never SCI values unless mask_by_illumination=True or an explicit
    mask_illumination_reference is supplied.
    Build WCSs on small synthetic 2D ImageModels:
    no Detector1 calibration, reference-pixel correction, or ramp fitting.
    Raw IRS2 inputs use the science grid specified by SUBSIZE1/SUBSIZE2,
    checked against the raw shape; interleaved reference samples are not part
    of the output. Boxes always use RATE/science Image coordinates.

    None switches follow EXP_TYPE or exp_type_override: IFU on for NRS_IFU,
    MSA on for NRS_MSASPEC, failed-open flagging on for either. Existing
    MSA_FAILED_OPEN DQ flags are honored unless mask_failed_open=False.
    mask_MSA=True forces MSA geometry regardless of the header/override;
    it requires MSAMETFL, MSAMETID, PATT_NUM and the exposure's MSA file.
    input_dir defaults to the input FITS parent; it locates that MSA file.
    IFU and default MSA exclusions use the JWST 3.0 mask helper algorithms.
    mask_MSA_by_wcs=True replaces MSA bounding-box exclusions with curved
    footprints evaluated at detector pixel centers, removing the extra 0.5
    slit-coordinate margin at both ends of regular MOS slitlets. Fixed slits
    retain their nominal WCS limits. This option selects the MSA algorithm;
    mask_MSA/EXP_TYPE still control whether MSA masking runs. It uses no SCI
    values or empirical edge fitting, and retains the nominal +/-0.55 shutter
    extent and any gaps inside a slitlet's outer envelope as excluded regions.
    Other exclusions (boxes, fixed-slit strip, failed-open flags) still apply.
    Failed shutters are projected by MSAFlagOpenStep using CRDS msaoper
    references.
    Missing metadata or references raises an error instead of guessing.

    mask_grating_extrapolate=None (default) enables aperture continuation for
    IFU G235M/F170LP, G235H/F170LP, and G140H/F100LP, the configurations
    with recorded real-data extension checks. False disables it explicitly;
    True enables it for other grating configurations/modes as well.
    mask_grating_extrapolate=True extends each full WCS aperture redward only
    to the high-x detector edge using degree 1--3 boundary fits chosen by withheld-end
    prediction errors. Errors set conservative vertical padding; unreliable
    fits warn and widen. It requires curved WCS masks for MOS/fixed slits.
    No order identification is implied. Prism geometry is unchanged.
    Geometry caching is bypassed and fit diagnostics are saved in an adjacent
    .extrapolation.json file. Existing local illumination refinement can
    add wings around these extended aperture bands. Polynomial continuation never adds
    blueward exclusions; existing exclusions and local measured-illumination
    refinement remain in effect.
    When no nominal IFU WCS falls on the detector, this option instead projects
    IFU orders 1--3 over the wavelength-reference blue bound through 5.3 microns
    and records a warning. This is provisional optical geometry, not measured
    filter transmission or validated higher-order calibration.

    wcs_reference_overrides optionally supplies explicit AssignWcs reference
    paths (and msaoper for failed-open flagging). Geometry caching is bypassed
    when overrides are supplied. Otherwise normal CRDS selection applies.

    With extrapolation disabled, IFU/MSA/fixed-slit geometry follows the
    pipeline's nominal science WCS domain. No separate higher-order footprints
    are constructed. Higher-order traces or first-order emission beyond that
    domain can therefore remain unmasked and affect the flicker-noise fit.
    The default 30-pixel illumination corridor is only local refinement; it
    cannot recover missing orders hundreds of pixels beyond those footprints.
    See HIGHER_ORDER_MASK_CONTEXT.md for the F170LP/G235M example and plan.

    mask_FS_by_wcs=True independently masks individual fixed-slit footprints
    using a separate NRS_FIXEDSLIT WCS, preserving dark gaps between slits.
    Full-frame inputs include all fixed slits projected on the detector;
    subarrays use JWST's corresponding fixed-slit selection. Nominal slit
    limits are retained, without subtracting the MOS extraction margin.
    mask_FS_strip covers every column in inclusive DS9 Image rows 950–1108.
    None enables it for IFU, or MOS without quadrant-5 fixed slits, unless
    mask_FS_by_wcs=True replaces this automatic strip. Explicit
    mask_FS_strip=True, mask_FSprism=True, and mask_box still apply and can
    cover the gaps; disable these broad exclusions to recover those pixels.
    mask_FSprism covers inclusive origin-1 corners (846, 952, 1571, 1105).
    mask_box adds an independent box (x1,y1,x2,y2), inclusive, with origin 1
    by default or mask_origin=0. Reversed/fractional corners are supported:
    exclude integer pixel centers inside the box, clipped to the image.
    These coordinates remain local to the science image for subarrays.

    mask_by_illumination=True refines WCS slit masks using the input itself.
    For UNCAL/4D ramps, it reads selected group pairs from up to four evenly
    spaced integrations, skips group zero, strips IRS2 reference samples,
    removes provisional per-channel/broad median offsets, and combines the
    differences into a DN/s illumination proxy. No Detector1 steps or ramp
    fitting run. At least three groups and positive TGROUP are required;
    short ramps and faint sky may be less reliable. For a 2D RATE input, it
    uses that image directly. RATEINTS/3D inputs are not supported for this
    option. Single/cube masks use the same combined illumination proxy.
    This is the tested ramp refinement with no additional endpoint padding.
    By default no NaN/zero exclusions or image-statistics clipping occur.
    Alternatively, mask_illumination_reference supplies a separate 2D RATE
    FITS for the same exposure/grid; do not combine it with mask_by_illumination.
    With WCS slit masking enabled, either illumination method adds measured
    illumination outside the nominal footprints: a 3-row by 5-column median
    filter is compared with background median + mask_illumination_nsigma times
    the robust MAD scatter, estimated from currently usable, finite reference
    pixels (excluding DO_NOT_USE flags for RATE images). The ramp proxy is not
    gain calibrated; its background threshold is estimated in DN/s.
    Search is bounded around each actual IFU/MSA/fixed-slit WCS footprint:
    mask_illumination_dispersion_pixels=30 permits endpoint continuation along
    the local trace slope, and mask_illumination_cross_dispersion_pixels=3
    permits small vertical adjustments. NIRSpec dispersion is along detector x.
    These are eligibility limits, not automatic exclusions. The original
    footprints remain fixed; distant glow and unrecorded slits cannot seed growth.
    Optional mask_illumination_padding=(x, y) further restricts this corridor to
    expanded WCS bounding boxes. Nonfinite or DO_NOT_USE reference pixels in
    the search region are also excluded.
    Existing exclusions are never removed. This is empirical refinement,
    not a guarantee of covering arbitrarily faint sky; inspect the result.
    Use this mask with clean_flicker_noise fit_method='fft', autoparam=False,
    background_method='median' or 'model', n_sigma>0 to clip group differences
    in Detector1. background_method=None disables that additional clipping.

    single_mask=True writes a 2D mask for all integrations; False writes
    identical geometric planes (plus per-plane existing flags if supplied).
    Input values/metadata are unchanged, saved metadata retain the original
    EXP_TYPE, and existing outputs require overwrite=True.

    ifu_geometry_cache_dir optionally caches IFU geometry on disk across
    exposures/workers. Reference selection, instrument/subarray metadata,
    mask options, CRDS context and software identify an entry. Pointing,
    nod number and velocity corrections are deliberately omitted: this
    reuses the first footprint, including possible small boundary differences.
    Illumination refinement and existing input DQ are never cached. MOS
    geometry is not cached. None disables caching for exact comparisons.
    """
    output = Path(output_file).expanduser().resolve()
    if wcs_reference_overrides is not None:
        if not isinstance(wcs_reference_overrides, dict):
            raise TypeError('wcs_reference_overrides must be a dictionary')
        unknown = set(wcs_reference_overrides) - set(AssignWcsStep.reference_file_types) - {'msaoper'}
        if unknown:
            raise ValueError(f'Unknown WCS reference types: {sorted(unknown)}')
        wcs_reference_overrides = {k: str(v) for k, v in wcs_reference_overrides.items()}
    if output.exists() and not overwrite:
        raise FileExistsError(output)
    supported = {'NRS_IFU', 'NRS_MSASPEC', 'NRS_FIXEDSLIT', 'NRS_BRIGHTOBJ'}
    if exp_type_override is not None:
        if not isinstance(exp_type_override, str):
            raise TypeError('exp_type_override must be a string or None')
        exp_type_override = exp_type_override.strip().upper()
        if exp_type_override not in supported:
            raise ValueError(f'Unsupported exp_type_override: {exp_type_override!r}')
    for value in (mask_IFU, mask_MSA, mask_FS_strip, mask_failed_open, mask_grating_extrapolate):
        if value is not None and not isinstance(value, (bool, np.bool_)):
            raise TypeError('Mask switches must be True, False, or None')
    mask_failed_open = None if mask_failed_open is None else bool(mask_failed_open)
    for value in (mask_FSprism, mask_box, single_mask, mask_MSA_by_wcs, mask_FS_by_wcs, mask_by_illumination):
        if not isinstance(value, (bool, np.bool_)):
            raise TypeError('Mask algorithm switches and single_mask must be boolean')
    if mask_by_illumination and mask_illumination_reference is not None:
        raise ValueError('Use mask_by_illumination or mask_illumination_reference, not both')
    if mask_illumination_reference is not None:
        if not isinstance(mask_illumination_reference, (str, Path)):
            raise TypeError('mask_illumination_reference must be a RATE FITS path')
    if mask_by_illumination or mask_illumination_reference is not None:
        for name, value in (('mask_illumination_dispersion_pixels', mask_illumination_dispersion_pixels),
                            ('mask_illumination_cross_dispersion_pixels', mask_illumination_cross_dispersion_pixels)):
            if (isinstance(value, (bool, np.bool_)) or not np.isscalar(value)
                    or not np.isfinite(value) or value < 0 or value != int(value)):
                raise ValueError(f'{name} must be a nonnegative integer')
        if (isinstance(mask_illumination_nsigma, (bool, np.bool_))
                or not np.isfinite(mask_illumination_nsigma) or mask_illumination_nsigma <= 0):
            raise ValueError('mask_illumination_nsigma must be finite and positive')
        padding = None
        if mask_illumination_padding is not None:
            padding = np.asarray(mask_illumination_padding, dtype=float)
            if (padding.shape != (2,) or not np.all(np.isfinite(padding))
                    or np.any(padding < 0) or np.any(padding != np.floor(padding))):
                raise ValueError('mask_illumination_padding must contain two nonnegative integers (x, y)')
            padding = tuple(int(value) for value in padding)
    boxes = [(845, 951, 1570, 1104)] if mask_FSprism else []
    if mask_box:
        if isinstance(mask_origin, (bool, np.bool_)) or mask_origin not in (0, 1):
            raise ValueError('mask_origin must be 0 or 1')
        corners = np.asarray(mask_corners, dtype=float)
        if corners.shape != (4,) or not np.all(np.isfinite(corners)):
            raise ValueError('mask_corners must contain four finite coordinates')
        x1, y1, x2, y2 = corners - mask_origin
        boxes.append((min(x1, x2), min(y1, y2), max(x1, x2), max(y1, y2)))

    with ExitStack() as stack:
        def own(model):
            stack.callback(model.close)
            return model

        if isinstance(input_data, (str, Path)):
            filename = Path(input_data).expanduser().resolve()
            if filename == output:
                raise ValueError('Output must differ from the science input')
            metadata = read_metadata(filename, flatten=False)
            image = own(datamodels.ImageModel())
            image.update({'meta': metadata['meta']})
            # No SCI access; non-memmapped DQ supports FITS unsigned BZERO.
            hdus = stack.enter_context(fits.open(filename, memmap=False))
            header = hdus['SCI'].header
            ndim = header['NAXIS']
            raw_shape = (header['NAXIS2'], header['NAXIS1'])
            nints = header.get('NAXIS4', 1) if ndim == 4 else header.get('NAXIS3', 1)
            directory = str(filename.parent) if input_dir is None else str(Path(input_dir).expanduser())
            existing_dq = None
            if mask_failed_open is not False:
                for name in ('DQ', 'PIXELDQ'):
                    if name in hdus:
                        existing_dq = hdus[name].data
                        break
        else:
            if not isinstance(input_data, (datamodels.ImageModel, datamodels.IFUImageModel,
                                           datamodels.CubeModel, datamodels.RampModel,
                                           datamodels.Level1bModel)):
                raise TypeError('Input must be an unextracted NIRSpec image, cube or ramp')
            image = own(datamodels.ImageModel())
            image.update(input_data)
            shape = input_data.data.shape
            ndim, raw_shape = len(shape), shape[-2:]
            nints = shape[0] if ndim > 2 else 1
            directory = '' if input_dir is None else str(Path(input_dir).expanduser())
            existing_dq = None
            if mask_failed_open is not False:
                existing_dq = getattr(input_data, 'dq', None)
                if existing_dq is None:
                    existing_dq = getattr(input_data, 'pixeldq', None)
        if ndim not in (2, 3, 4):
            raise ValueError('SCI must be a 2D image, 3D rateints cube, or 4D ramp')
        header_type = (image.meta.exposure.type or '').upper()
        selected = header_type if exp_type_override is None else exp_type_override
        if selected not in supported:
            raise ValueError(f'Unsupported NIRSpec EXP_TYPE: {selected!r}')
        detector = image.meta.instrument.detector
        if detector not in ('NRS1', 'NRS2'):
            raise ValueError('NIRSpec DETECTOR must be NRS1 or NRS2')

        science_shape = raw_shape
        if max(raw_shape) > 2048:
            # IRS2 DMS raw arrays include interleaved samples along rows.
            # Only metadata are needed to describe the post-refpix grid.
            science_shape = (image.meta.subarray.ysize, image.meta.subarray.xsize)
            if (ndim != 4 or 'IRS2' not in (image.meta.exposure.readpatt or '').upper()
                    or raw_shape != (3200, 2048) or science_shape != (2048, 2048)):
                raise ValueError('Unrecognized expanded raw grid; expected full-frame NIRSpec IRS2')
            if existing_dq is not None and existing_dq.shape[-2:] == raw_shape:
                from jwst.refpix.irs2_subtract_reference import make_irs2_mask
                normal = image.meta.exposure.nrs_normal
                reference = image.meta.exposure.nrs_reference
                rows = make_irs2_mask(3200, 2048, 16 if normal is None else normal,
                                     4 if reference is None else reference)
                if detector == 'NRS2':
                    rows = rows[::-1]
                existing_dq = existing_dq[..., rows, :]
        image.data = np.zeros(science_shape, dtype=np.float32)
        image.dq = np.zeros(science_shape, dtype=np.uint32)
        image.meta.wcs = None
        # Update sizes to describe the synthetic science image, not raw samples.
        image.meta.subarray.ysize, image.meta.subarray.xsize = science_shape
        use_ifu = selected == 'NRS_IFU' if mask_IFU is None else mask_IFU
        use_msa = selected == 'NRS_MSASPEC' if mask_MSA is None else mask_MSA
        grating = (image.meta.instrument.grating or '').upper()
        filter_name = (image.meta.instrument.filter or '').upper()
        if mask_grating_extrapolate is None:
            extrapolate = bool(use_ifu) and (grating, filter_name) in {
                ('G235M', 'F170LP'), ('G235H', 'F170LP'), ('G140H', 'F100LP')}
        else:
            extrapolate = bool(mask_grating_extrapolate) and grating != 'PRISM'
        if extrapolate and use_msa and not mask_MSA_by_wcs:
            raise ValueError('Grating extrapolation for MOS requires mask_MSA_by_wcs=True')
        if extrapolate and not (use_ifu or use_msa or mask_FS_by_wcs):
            raise ValueError('Grating extrapolation requires WCS aperture masks')
        extrapolation_reports = []
        if (mask_by_illumination or mask_illumination_reference is not None) and not (use_ifu or mask_FS_by_wcs or (use_msa and mask_MSA_by_wcs)):
            raise ValueError('Illumination refinement requires an enabled WCS slit mask')
        if mask_by_illumination and ndim == 3:
            raise ValueError('mask_by_illumination requires a 4D ramp or a 2D RATE input')
        models = {}

        def wcs_step():
            if wcs_reference_overrides is None:
                return AssignWcsStep()
            return AssignWcsStep(**{f'override_{k}': v for k, v in wcs_reference_overrides.items()
                                   if k != 'msaoper'})

        def geometry(mode):
            if mode not in models:
                model = own(image.copy())
                model.meta.exposure.type = mode
                step = wcs_step()
                step.input_dir = directory
                model = own(step.run(model))
                if getattr(model.meta, 'wcs', None) is None:
                    raise RuntimeError(f'WCS assignment did not produce a WCS for {mode}')
                models[mode] = model
            return models[mode]

        mask = np.ones(science_shape, dtype=bool)
        trace_region = np.zeros(science_shape, dtype=bool)
        limits = (int(mask_illumination_dispersion_pixels), int(mask_illumination_cross_dispersion_pixels))
        slit_boxes = []
        cache_path = None
        cache_hit = False
        if ifu_geometry_cache_dir is not None and selected == 'NRS_IFU' and use_ifu and not use_msa and not extrapolate and wcs_reference_overrides is None:
            cache_path = _ifu_cache_path(image, ifu_geometry_cache_dir, science_shape, dict(
                fs_wcs=bool(mask_FS_by_wcs), fs_strip=mask_FS_strip,
                failed_open=mask_failed_open is not False, boxes=boxes, illumination_limits=limits))
            # Lock only shared geometry: workers wait for its first complete
            # publication, then refine their own arrays independently.
            import fcntl
            lock = stack.enter_context(open(str(cache_path) + '.lock', 'a'))
            fcntl.flock(lock, fcntl.LOCK_EX)
            if cache_path.exists():
                with np.load(cache_path, allow_pickle=False) as cached:
                    mask = cached['mask'].astype(bool)
                    slit_boxes = cached['boxes'].tolist()
                    trace_region = cached['trace_region'].astype(bool)
                if mask.shape != science_shape:
                    raise ValueError(f'Cached IFU geometry has wrong shape: {cache_path}')
                cache_hit = True
                print(f'IFU geometry cache hit: {cache_path.name}', flush=True)
        if use_ifu and not cache_hit:
            # clean_flicker_noise.mask_ifu_slices algorithm, without its DQ map.
            try:
                ifu_wcs_list = nirspec.nrs_ifu_wcs(geometry('NRS_IFU'))
            except NoDataOnDetectorError:
                if not extrapolate:
                    raise
                import warnings
                from .nirspec_order_geometry import mask_ifu_optical_orders
                warnings.warn('No nominal IFU footprint on detector; using provisional '
                              'optical projections of orders 1--3', RuntimeWarning)
                step = wcs_step()
                references = {k: step.get_reference_file(image, k)
                              for k in ('camera', 'collimator', 'fpa', 'disperser',
                                        'ifuslicer', 'ifupost', 'wavelengthrange')}
                extrapolation_reports.append(dict(method='optical fallback', warning=True,
                                                    references=references))
                extrapolation_reports.extend(mask_ifu_optical_orders(
                    image, references, mask, trace_region, search_padding=limits[1]))
                ifu_wcs_list = []
                # Failed-open flagging builds its own slit WCS; provide metadata
                # and references even though the nominal IFU WCS is absent.
                fallback = own(image.copy())
                fallback.meta.exposure.type = 'NRS_IFU'
                for k in step.reference_file_types:
                    getattr(fallback.meta.ref_file, k).name = step.get_reference_file(fallback, k)
                models['NRS_IFU'] = fallback
            for wcs in ifu_wcs_list:
                x, y = gwcs.wcstools.grid_from_bounding_box(wcs.bounding_box, step=(1, 1), center=True)
                valid = ~np.isnan(wcs(x, y)[2])
                mask[y[valid].astype(int), x[valid].astype(int)] = False
                # Keep each IFU slice separate so neighboring dark strips survive.
                xx, yy = x[valid].astype(int), y[valid].astype(int)
                if len(xx):
                    xlo, xhi, ylo, yhi = xx.min(), xx.max(), yy.min(), yy.max()
                    footprint = np.zeros((yhi - ylo + 1, xhi - xlo + 1), bool)
                    footprint[yy - ylo, xx - xlo] = True
                    _add_trace_search_region(trace_region, footprint, xlo, ylo, *limits)
                    if extrapolate:
                        extrapolation_reports.append(extend_aperture(
                            mask, trace_region, footprint, xlo, ylo,
                            search_padding=limits[1], label=f'IFU slice {len(slit_boxes)}'))
                    slit_boxes.append(((float(xlo), float(xhi)), (float(ylo), float(yhi))))
        if use_msa:
            model = geometry('NRS_MSASPEC')
            if mask_MSA_by_wcs:
                if extrapolate:
                    slit_boxes.extend(_mask_slits_by_wcs(model, mask, trace_region, limits, extrapolation_reports))
                else:
                    slit_boxes.extend(_mask_slits_by_wcs(model, mask, trace_region, limits))
            else:
                # Original clean_flicker_noise.mask_slits algorithm.
                slits = model.meta.wcs.get_transform('gwa', 'slit_frame').slit_ids
                for slit in slits:
                    bbox = model.meta.wcs.bounding_box[slit]
                    xlo, xhi = to_index(bbox[0])
                    ylo, yhi = to_index(bbox[1])
                    mask[ylo:yhi, xlo:xhi] = False
        if mask_FS_by_wcs and not cache_hit:
            if extrapolate:
                slit_boxes.extend(_mask_slits_by_wcs(geometry('NRS_FIXEDSLIT'), mask, trace_region, limits, extrapolation_reports))
            else:
                slit_boxes.extend(_mask_slits_by_wcs(geometry('NRS_FIXEDSLIT'), mask, trace_region, limits))
        if not cache_hit and (mask_failed_open is True or (mask_failed_open is None and selected in ('NRS_IFU', 'NRS_MSASPEC'))):
            step = (MSAFlagOpenStep(override_msaoper=wcs_reference_overrides['msaoper'])
                    if wcs_reference_overrides and 'msaoper' in wcs_reference_overrides
                    else MSAFlagOpenStep())
            step.input_dir = directory
            try:
                flagged = own(step.run(geometry(selected)))
            except json.JSONDecodeError as exc:
                reference = getattr(step, 'reference_name', None)
                size = 'unknown'
                if reference:
                    try:
                        size = Path(reference).stat().st_size
                    except OSError:
                        pass
                exc.add_note(
                    f'Failed-open shutter msaoper reference: {reference!r}; '
                    f'file size: {size} bytes. Validate this file as JSON; '
                    'if invalid, quarantine it and re-download it through CRDS '
                    'in a single process before restarting the worker pool.'
                )
                raise
            if flagged.meta.cal_step.msa_flagging != 'COMPLETE':
                raise RuntimeError('Failed-open shutter flagging did not complete')
            mask[(flagged.dq & dqflags.pixel['MSA_FAILED_OPEN']) != 0] = False
        strip = mask_FS_strip
        if strip is None:
            strip = selected == 'NRS_IFU' and not mask_FS_by_wcs
            if selected == 'NRS_MSASPEC' and not mask_FS_by_wcs:
                slits = geometry('NRS_MSASPEC').meta.wcs.get_transform('gwa', 'slit_frame').slits
                strip = not any(s.quadrant == 5 for s in slits)
        if strip:
            mask[949:1108, :] = False
        height, width = science_shape
        for xlo, ylo, xhi, yhi in boxes:
            xstart, xstop = np.clip((np.ceil(xlo), np.floor(xhi) + 1), 0, width).astype(int)
            ystart, ystop = np.clip((np.ceil(ylo), np.floor(yhi) + 1), 0, height).astype(int)
            mask[ystart:ystop, xstart:xstop] = False
        if cache_path is not None and not cache_hit:
            # Atomic publication; incomplete entries are never read.
            with tempfile.NamedTemporaryFile(dir=cache_path.parent, suffix='.npz', delete=False) as temp:
                temporary = Path(temp.name)
                try:
                    np.savez(temp, mask=mask, boxes=np.asarray(slit_boxes, dtype=float), trace_region=trace_region)
                    temp.flush()
                    os.replace(temporary, cache_path)
                finally:
                    temporary.unlink(missing_ok=True)
            print(f'IFU geometry cache created: {cache_path.name}', flush=True)
        if cache_path is not None:
            # np.load/astype created private arrays. Geometry is now immutable;
            # raw/rate illumination, DQ refinement and saving must not serialize.
            fcntl.flock(lock, fcntl.LOCK_UN)
        if mask_by_illumination and ndim == 4:
            proxy = _ramp_illumination_proxy(input_data, mask, image.meta)
            _extend_mask_from_rate(mask, proxy, np.isfinite(proxy), slit_boxes,
                                   padding, mask_illumination_nsigma, trace_region)
        elif mask_illumination_reference is not None or mask_by_illumination:
            source = input_data if mask_by_illumination else mask_illumination_reference
            reference_path = Path(source).expanduser().resolve() if isinstance(source, (str, Path)) else None
            if reference_path == output:
                raise ValueError('Output must differ from illumination reference')
            # Copy a supplied model so closing the owned reference cannot close input_data.
            reference = own(datamodels.open(reference_path) if reference_path is not None else source.copy())
            reference_shape = getattr(reference.data, 'shape', None)
            if not isinstance(reference, (datamodels.ImageModel, datamodels.IFUImageModel)) or reference_shape != science_shape:
                raise ValueError(
                    'Illumination reference must be a 2D RATE image on the same science grid; '
                    f'got {type(reference).__name__} with shape {reference_shape}, '
                    f'expected shape {science_shape}')
            for group, names in (
                ('instrument', ('name', 'detector', 'filter', 'grating')),
                ('subarray', ('xstart', 'ystart', 'xsize', 'ysize')),
                ('exposure', ('start_time',)),
            ):
                for name in names:
                    expected = getattr(getattr(image.meta, group), name)
                    actual = getattr(getattr(reference.meta, group), name)
                    if expected is not None and actual != expected:
                        raise ValueError(f'Illumination reference has mismatched {group}.{name}')
            valid = np.isfinite(reference.data)
            if reference.dq is not None:
                valid &= (reference.dq & dqflags.pixel['DO_NOT_USE']) == 0
            _extend_mask_from_rate(mask, reference.data, valid, slit_boxes,
                                   padding, mask_illumination_nsigma, trace_region)
        if not single_mask:
            mask = np.broadcast_to(mask, (nints, height, width)).copy()
        if existing_dq is not None:
            if existing_dq.shape[-2:] != science_shape:
                raise ValueError('Existing DQ does not match the science grid')
            failed = (existing_dq & dqflags.pixel['MSA_FAILED_OPEN']) != 0
            if failed.ndim == 3 and single_mask:
                failed = np.any(failed, axis=0)
            if failed.ndim == 3 and not single_mask and failed.shape[0] != nints:
                raise ValueError('Existing DQ integration count does not match SCI')
            mask &= ~failed
        cls = datamodels.ImageModel if single_mask else datamodels.CubeModel
        saved = own(cls(data=mask.astype(np.float32)))
        saved.update(image)
        saved.save(str(output))
        if extrapolate:
            output.with_suffix('.extrapolation.json').write_text(json.dumps(
                dict(method='WCS aperture extrapolation', detector=detector,
                     crds_context=os.environ.get('CRDS_CONTEXT'),
                     reference_overrides=wcs_reference_overrides,
                     fits=extrapolation_reports), indent=2) + '\n')
        return mask

"""Diagnostic IFU aperture masks from order-aware NIRSpec optical models.

The wavelength bounds are coverage policy, not a filter throughput model.
Higher-order and extrapolated geometry require validation on each setup.
"""
import numpy as np
from jwst.assign_wcs import nirspec


def mask_ifu_optical_orders(image, references, mask, region=None, *, orders=(1, 2, 3),
                            physical_bounds=None, padding=1, search_padding=3):
    """Union projected full IFU apertures, including detectors lacking nominal WCS.

    Defaults to the filter/grating reference's blue bound through the NIRSpec
    5.3-micron model limit. No filter leakage or diffraction efficiency is
    inferred. Order sign follows the nominal JWST reference convention.
    """
    reports = []
    height, width = mask.shape
    with image.copy() as model:
        nominal_order, nominal_range = nirspec.get_spectral_order_wrange(
            model, references['wavelengthrange'])
        if nominal_order == 0:
            raise ValueError('Optical order projection requires a grating')
        bounds = (nominal_range[0], 5.3e-6) if physical_bounds is None else physical_bounds
        model.meta.wcsinfo.waverange_start, model.meta.wcsinfo.waverange_end = nominal_range
        disperser = nirspec.get_disperser(model, references['disperser'])
        detector = model.meta.instrument.detector
        det2gwa = (nirspec._fix_slit_name(nirspec.dms_to_sca(model), 0) |
                   nirspec._fix_slit_name(nirspec.detector_to_gwa(references, detector, disperser), 0))
        wave = np.linspace(*bounds, 2400)
        for order in orders:
            model.meta.wcsinfo.spectral_order = int(np.sign(nominal_order)) * order
            group = nirspec.gwa_to_ifuslit(list(range(30)), model, disperser,
                                         references, (-.55, .55))
            for name in group.slit_ids:
                forward = det2gwa | group.get_model(name)
                inverse = group.get_model(name).inverse | det2gwa.inverse
                edges = [inverse(np.zeros_like(wave), np.full_like(wave, y), wave)
                         for y in (-.55, .55)]
                xy = np.asarray(edges)
                xlo = max(0, int(np.floor(np.nanmin(xy[:, 0])))-3)
                xhi = min(width, int(np.ceil(np.nanmax(xy[:, 0])))+4)
                border = padding + search_padding + 1
                ylo = max(0, int(np.floor(np.nanmin(xy[:, 1])))-border)
                yhi = min(height, int(np.ceil(np.nanmax(xy[:, 1])))+border+1)
                if xlo >= xhi or ylo >= yhi:
                    continue
                yy, xx = np.mgrid[ylo:yhi, xlo:xhi]
                _, sy, wavelength = forward(xx, yy)
                footprint = (np.isfinite(wavelength) & np.isfinite(sy) &
                             (sy >= -.55) & (sy <= .55) &
                             (wavelength >= bounds[0]) & (wavelength <= bounds[1]))
                # Narrow vertical padding, without expanding wavelength admission.
                if padding:
                    from scipy.ndimage import maximum_filter1d
                    footprint = maximum_filter1d(footprint, 2*padding+1, axis=0,
                                                  mode='constant', cval=0)
                mask[ylo:yhi, xlo:xhi][footprint] = False
                if region is not None:
                    from scipy.ndimage import maximum_filter1d
                    region[ylo:yhi, xlo:xhi] |= maximum_filter1d(
                        footprint, 2*search_padding+1, axis=0, mode='constant', cval=0)
                reports.append(dict(method='optical IFU projection', order=order,
                                    slice=int(name), pixels=int(footprint.sum()),
                                    physical_bounds_meters=list(bounds), padding_pixels=padding))
    return reports

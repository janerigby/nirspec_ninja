"""Conservative redward detector-x continuation of full aperture footprints.

This is geometry, not an order identification or wavelength calibration.
"""
import warnings
import numpy as np
from numpy.polynomial import Polynomial


def extend_aperture(mask, region, footprint, xstart, ystart, *, padding=1,
                    search_padding=3, label='aperture'):
    """Fit aperture edges, validate withheld endpoints, mask beyond the red x endpoint.

    Select the simplest degree (1--3) within half a pixel of the best
    bidirectional withheld-end error. Pixel quantization gets a one-pixel
    floor. Failed validation increases padding and emits a warning. Never
    removes existing exclusions. Continuation adds pixels only at larger
    detector x than the nominal footprint; no blueward continuation is made.
    Returns serializable fit diagnostics.
    """
    columns = np.flatnonzero(footprint.any(axis=0))
    if not len(columns):
        return None
    lo = np.array([np.flatnonzero(footprint[:, x])[0] for x in columns]) + ystart
    hi = np.array([np.flatnonzero(footprint[:, x])[-1] for x in columns]) + ystart
    x = columns + xstart
    # Wavelength clipping can truncate aperture heights at spectral endpoints.
    widths = hi - lo + 1
    typical = float(np.median(widths))
    width_scatter = 1.4826 * float(np.median(np.abs(widths - typical)))
    # Some nominal IFU WCS footprints have extra endpoint islands that make
    # a few columns much wider. Retain those original exclusions, but do not
    # let their artificial envelope drive detector-edge extrapolation.
    keep = ((widths >= .8 * typical) &
            (widths <= typical + max(3., 6. * width_scatter)))
    if not keep.any():
        keep = np.ones(len(widths), dtype=bool)
    rejected_columns = int((~keep).sum())
    x, lo, hi = x[keep], lo[keep], hi[keep]
    edges = np.stack((lo, hi))
    scores, models = {}, {}
    n = len(x)
    if n >= 30:
        cut = max(6, n // 5)
        for degree in (1, 2, 3):
            errors = []
            for train, test in ((slice(cut, None), slice(None, cut)),
                                (slice(None, -cut), slice(-cut, None))):
                for edge in edges:
                    fit = Polynomial.fit(x[train], edge[train], degree)
                    errors.extend(np.abs(fit(x[test]) - edge[test]))
            scores[degree] = float(np.percentile(errors, 99))
            models[degree] = [Polynomial.fit(x, edge, degree) for edge in edges]
        best = min(scores.values())
        degree = min(d for d in scores if scores[d] <= best + .5)
        fits = models[degree]
        error = scores[degree]
    else:
        degree, error = 0, float(np.ptp(edges, axis=1).max())
        fits = [Polynomial([float(np.median(edge))]) for edge in edges]
    height, width = mask.shape
    xx = np.arange(width)
    outside = xx > columns[-1] + xstart
    distance = np.maximum(x[0] - xx, np.maximum(xx - x[-1], 0))
    span = max(float(x[-1] - x[0]), 1.)
    # Endpoint validation covers 20% of the measured span; inflate for longer
    # extrapolations, without spreading vertically merely because x is long.
    margin = float(padding) + error * np.maximum(1., distance / (.2 * span))
    lower, upper = fits[0](xx), fits[1](xx)
    unreliable = n < 30 or error > 2 or np.any(lower[outside] > upper[outside])
    if unreliable:
        margin += max(2., error)
        warnings.warn(f'{label}: unreliable aperture extrapolation; wider mask '
                      f'(n={n}, withheld error={error:.2f} px)', RuntimeWarning)
    for col in xx[outside]:
        a, b = sorted((lower[col], upper[col]))
        start = max(0, int(np.ceil(a - margin[col])))
        stop = min(height, int(np.floor(b + margin[col])) + 1)
        mask[start:stop, col] = False
        if region is not None:
            region[max(0, start-search_padding):min(height, stop+search_padding), col] = True
    return dict(label=label, direction='redward', degree=degree, withheld_error_pixels=error,
                candidate_errors=scores, fit_columns=[int(x[0]), int(x[-1])],
                nominal_columns=[int(columns[0]+xstart), int(columns[-1]+xstart)],
                rejected_width_columns=rejected_columns,
                max_padding_pixels=float(margin[outside].max()) if outside.any() else 0.,
                warning=bool(unreliable))

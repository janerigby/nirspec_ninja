""" Trace detection for rectified JWST S2D spectra."""
# THIS IS Codex's improved version of Rosalia O'Briens algorithm
#GPT-5.6-Sol, Oct 2026

import warnings

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import ImageNormalize, ZScaleInterval
from scipy.special import ndtri_exp
from scipy.stats import norm


NIRSPEC_FIXED_SLIT_LENGTH_ARCSEC = {
    "S200A1": 3.2, "S200A2": 3.2, "S200B1": 3.2,
    "S400A1": 3.65, "S1600A1": 1.6,
}
NIRSPEC_SPATIAL_PIXEL_SCALE_ARCSEC = 0.1


def _failure(reason, **extra):
    """Return the common structured result for an unsuccessful measurement."""
    result = {
        "status": "failed", "reason": str(reason), "measurement": None,
        "significance": np.nan, "raw_significance": np.nan,
        "center": np.nan, "sigma": np.nan, "fwhm": np.nan,
        "amplitude": np.nan,
    }
    result.update(extra)
    return result


def _joint_good_mask(data, error=None, dq=None, dq_mask_bits=1):
    """Return one mask used consistently for SCI, ERR, and DQ."""
    good = np.isfinite(data)
    if error is not None:
        good &= np.isfinite(error) & (error > 0)
    if dq is not None:
        good &= (np.asarray(dq, dtype=np.uint64) & np.uint64(dq_mask_bits)) == 0
    return good


def _collapsed_profile(data, error=None, dq=None, dq_mask_bits=1):
    """Median spatial profile and approximate error on each row median."""
    good = _joint_good_mask(data, error, dq, dq_mask_bits)
    masked_data = np.where(good, data, np.nan)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        profile = np.nanmedian(masked_data, axis=1)
    count = np.sum(good, axis=1)
    if error is None:
        residual = masked_data - profile[:, None]
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            row_center = np.nanmedian(residual, axis=1)[:, None]
            scatter = 1.4826 * np.nanmedian(np.abs(residual - row_center), axis=1)
        profile_error = 1.2533 * scatter / np.sqrt(np.maximum(count, 1))
    else:
        variance_sum = np.sum(np.where(good, error**2, 0.0), axis=1)
        profile_error = 1.2533 * np.sqrt(variance_sum) / np.maximum(count, 1)
    valid_error = np.isfinite(profile_error) & (profile_error > 0) & (count > 0)
    if np.any(valid_error):
        fallback = np.nanmedian(profile_error[valid_error])
        profile_error[(count > 0) & ~valid_error] = fallback
    profile_error[count == 0] = np.nan
    return profile, profile_error, good


def _fit_template(profile, profile_error, center, sigma):
    """Fit background plus a Gaussian using stable weighted least squares."""
    y = np.arange(profile.size, dtype=float)
    template = np.exp(-0.5 * ((y - center) / sigma) ** 2)
    good = np.isfinite(profile) & np.isfinite(profile_error) & (profile_error > 0)
    if good.sum() < 4:
        return None
    yg = y[good]
    scale = max(profile.size - 1, 1)
    design = np.column_stack((np.ones(yg.size), (yg - y.mean()) / scale, template[good]))
    weighted_design = design / profile_error[good, None]
    weighted_profile = profile[good] / profile_error[good]
    try:
        coefficients, _, rank, singular_values = np.linalg.lstsq(
            weighted_design, weighted_profile, rcond=None
        )
    except np.linalg.LinAlgError:
        return None
    if rank < design.shape[1] or not np.all(np.isfinite(coefficients)):
        return None
    if singular_values[-1] <= 0 or singular_values[0] / singular_values[-1] > 1e12:
        return None
    fitted = design @ coefficients
    chi2 = np.sum(((profile[good] - fitted) / profile_error[good]) ** 2)
    dof = good.sum() - design.shape[1]
    error_scale = max(1.0, np.sqrt(chi2 / dof))
    try:
        covariance = np.linalg.inv(weighted_design.T @ weighted_design)
    except np.linalg.LinAlgError:
        return None
    amplitude_variance = covariance[2, 2]
    if not np.isfinite(amplitude_variance) or amplitude_variance <= 0:
        return None
    amplitude_error = np.sqrt(amplitude_variance) * error_scale
    raw_significance = coefficients[2] / amplitude_error
    if not np.isfinite(raw_significance):
        return None
    background = coefficients[0] + coefficients[1] * (y - y.mean()) / scale
    model = background + coefficients[2] * template
    return raw_significance, coefficients[2], background, model


def _profile_fraction_on_detector(center, sigma, ny):
    """Integrated fraction of a Gaussian falling within pixel boundaries."""
    return norm.cdf((ny - 0.5 - center) / sigma) - norm.cdf((-0.5 - center) / sigma)


def _corrected_significance(raw_significance, trials):
    """Bonferroni-correct a one-sided significance without tail underflow."""
    log_false_alarm = min(0.0, norm.logsf(raw_significance) + np.log(trials))
    return max(0.0, float(-ndtri_exp(log_false_alarm)))


def _detect_source_s2d_detailed(
    spectrum_2d, error_2d=None, dq_2d=None, widths=None,
    slit_length_pixels=None, max_slit_length_fraction=0.75,
    min_profile_fraction=0.30, dq_mask_bits=1, n_blocks=6,
    plot=False, outfile=None, allowed_source_center_rows=None,
):
    """Measure the strongest continuum trace and return a structured result.

    allowed_source_center_rows is an optional 1D Boolean mask of spatial rows.
    True permits a fitted source center on that row. It restricts candidate
    centers, not the pixels contributing to the fit. None preserves the full
    original search, including centers just outside the detector.
    """
    data = np.asarray(spectrum_2d, dtype=float)
    error = None if error_2d is None else np.asarray(error_2d, dtype=float)
    dq = None if dq_2d is None else np.asarray(dq_2d)
    if data.ndim != 2:
        return _failure("SCI data are not a two-dimensional array")
    if error is not None and error.shape != data.shape:
        return _failure("ERR and SCI arrays have different shapes")
    if dq is not None and dq.shape != data.shape:
        return _failure("DQ and SCI arrays have different shapes")
    if data.shape[0] < 4 or data.shape[1] < 1:
        return _failure("SCI array is too small to fit a trace")
    if n_blocks < 1:
        return _failure("n_blocks must be at least one")

    ny = data.shape[0]
    allowed_centers = None
    if allowed_source_center_rows is not None:
        allowed_centers = np.asarray(allowed_source_center_rows)
        if allowed_centers.shape != (ny,) or allowed_centers.dtype != np.dtype(bool):
            return _failure("allowed_source_center_rows must be a 1D Boolean mask matching the spatial rows")
        if not np.any(allowed_centers):
            return _failure("allowed_source_center_rows permits no source centers")
    # Widths are Gaussian sigma, not FWHM. Enforce an additional physical
    # FWHM cap so a fitted source cannot be wider than the slit itself.
    fwhm_per_sigma = 2.35482
    width_cap = ny / fwhm_per_sigma
    if slit_length_pixels is not None:
        if not np.isfinite(slit_length_pixels) or slit_length_pixels <= 0:
            return _failure("slit_length_pixels must be finite and positive")
        width_cap = min(width_cap, max_slit_length_fraction * slit_length_pixels,
                        slit_length_pixels / fwhm_per_sigma)
    if widths is None:
        if width_cap < 0.7:
            return _failure("slit-width limit leaves no allowed Gaussian widths")
        widths = np.geomspace(0.7, width_cap, 14)
    widths = np.unique(np.asarray(widths, dtype=float))
    widths = widths[np.isfinite(widths) & (widths > 0) & (widths <= width_cap)]
    if widths.size == 0:
        return _failure("no requested Gaussian widths satisfy the width limits")

    profile, profile_error, good_pixels = _collapsed_profile(data, error, dq, dq_mask_bits)
    valid_rows = np.isfinite(profile) & np.isfinite(profile_error) & (profile_error > 0)
    if valid_rows.sum() < 4:
        return _failure(
            "fewer than four spatial rows contain usable SCI/ERR/DQ pixels",
            profile=profile, profile_error=profile_error,
        )

    best = None
    trials = 0
    rejected_geometry = 0
    for sigma in widths:
        # Search slightly beyond both edges so a trace whose centroid is just
        # off-detector can still be found, but retain only templates having the
        # requested minimum integrated area on the detector.
        center_margin = int(np.ceil(2.0 * sigma))
        centers = np.arange(-center_margin, ny + center_margin, dtype=float)
        for center in centers:
            if allowed_centers is not None:
                row = int(center)
                if row < 0 or row >= ny or not allowed_centers[row]:
                    continue
            if _profile_fraction_on_detector(center, sigma, ny) < min_profile_fraction:
                rejected_geometry += 1
                continue
            trials += 1
            fitted = _fit_template(profile, profile_error, center, sigma)
            if fitted is None or fitted[1] <= 0:
                continue
            if best is None or fitted[0] > best[0]:
                best = (*fitted, center, sigma)
    if trials == 0:
        return _failure("all templates failed the on-detector profile requirement")
    if best is None:
        return _failure(
            "no finite positive-amplitude template fit was found",
            profile=profile, profile_error=profile_error, trials=trials,
        )

    raw_sig, amplitude, background, model, center, sigma = best
    significance = _corrected_significance(raw_sig, trials)
    block_sigs = []
    for columns in np.array_split(np.arange(data.shape[1]), min(n_blocks, data.shape[1])):
        block_profile, block_error, _ = _collapsed_profile(
            data[:, columns],
            None if error is None else error[:, columns],
            None if dq is None else dq[:, columns], dq_mask_bits,
        )
        fitted = _fit_template(block_profile, block_error, center, sigma)
        block_sigs.append(np.nan if fitted is None else fitted[0])
    block_sigs = np.asarray(block_sigs)

    result = {
        "status": "ok", "reason": "continuum trace fit completed",
        "measurement": significance, "significance": significance,
        "raw_significance": float(raw_sig), "center": float(center),
        "sigma": float(sigma), "fwhm": float(2.35482 * sigma),
        "amplitude": float(amplitude), "background": background,
        "model": model, "profile": profile, "profile_error": profile_error,
        "block_significances": block_sigs,
        "blocks_above_3sigma": int(np.sum(block_sigs > 3)),
        "profile_fraction_on_detector": float(_profile_fraction_on_detector(center, sigma, ny)),
        "valid_pixel_fraction": float(np.mean(good_pixels)),
        "slit_length_pixels": None if slit_length_pixels is None else float(slit_length_pixels),
        "width_cap": float(width_cap), "trials": trials,
        "templates_rejected_by_geometry": rejected_geometry,
    }
    if plot:
        try:
            plot_detection(data, result, outfile=outfile)
        except Exception as exc:
            result["plot_status"] = "failed"
            result["plot_reason"] = f"{type(exc).__name__}: {exc}"
        else:
            result["plot_status"] = "ok"
            result["plot_reason"] = "diagnostic plot completed"
    return result


def detect_source_s2d(*args, **kwargs):
    """Measure a source in arrays, returning status, reason, and significance."""
    try:
        return _detect_source_s2d_detailed(*args, **kwargs)
    except Exception as exc:
        return _failure(f"unexpected {type(exc).__name__}: {exc}")


def plot_detection(spectrum_2d, result, outfile=None):
    """Plot the 2-D spectrum and fitted one-dimensional spatial profile."""
    fig = plt.figure(figsize=(12, 7))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.5])
    image_axis = fig.add_subplot(grid[0, :])
    profile_axis = fig.add_subplot(grid[1, 0])
    normalization = ImageNormalize(spectrum_2d, interval=ZScaleInterval(contrast=0.5))
    image_axis.imshow(spectrum_2d, aspect="auto", origin="lower", cmap="Greys_r", norm=normalization)
    image_axis.axhline(result["center"], color="red", linestyle="--", alpha=0.8)
    image_axis.axhspan(result["center"] - result["fwhm"] / 2,
                       result["center"] + result["fwhm"] / 2, color="red", alpha=0.12)
    image_axis.set(
        title=(f"2D Spectrum: {result['significance']:.1f} sigma corrected "
               f"({result['raw_significance']:.1f} sigma raw)"),
        xlabel="Dispersion Axis [Pixels]", ylabel="Spatial Axis [Pixels]",
    )
    y = np.arange(result["profile"].size)
    profile_axis.errorbar(y, result["profile"], yerr=result["profile_error"],
                          color="blue", lw=1.5, elinewidth=0.5, alpha=0.9,
                          label="Collapsed profile")
    profile_axis.plot(y, result["background"], "--", color="gray", label="Fitted background")
    profile_axis.plot(y, result["model"], color="red", label="Background + source model")
    profile_axis.axvspan(result["center"] - result["fwhm"] / 2,
                         result["center"] + result["fwhm"] / 2,
                         color="red", alpha=0.12, label=f"FWHM = {result['fwhm']:.1f} pixels")
    profile_axis.set(title="1D Spatial Profile", xlabel="Spatial Axis [Pixels]",
                     ylabel="Median Flux [MJy/sr]")
    profile_axis.legend(bbox_to_anchor=(1.02, 1), loc="upper left")
    fig.tight_layout()
    if outfile is None:
        plt.show()
    else:
        fig.savefig(outfile, dpi=180, bbox_inches="tight")
        plt.close(fig)


def wrap_detect_source_s2d(s2dfile, plot=False, outfile=None, **kwargs):
    """Measure a continuum source in a JWST S2D FITS file."""
    try:
        with fits.open(s2dfile) as hdul:
            if "SCI" not in hdul:
                return _failure("FITS file has no SCI extension")
            science = hdul["SCI"].data
            if science is None:
                return _failure("SCI extension contains no data")
            science = science.copy()
            error = hdul["ERR"].data.copy() if "ERR" in hdul else None
            dq = hdul["DQ"].data.copy() if "DQ" in hdul else None
            slit_name = str(hdul["SCI"].header.get("SLTNAME", "")).upper()
    except (OSError, ValueError, KeyError) as exc:
        return _failure(f"could not read S2D file: {type(exc).__name__}: {exc}")
    if "slit_length_pixels" not in kwargs:
        slit_length_arcsec = NIRSPEC_FIXED_SLIT_LENGTH_ARCSEC.get(slit_name)
        if slit_length_arcsec is None:
            return _failure(
                f"unrecognized or missing SLTNAME {slit_name!r}; supply slit_length_pixels explicitly"
            )
        kwargs["slit_length_pixels"] = (
            slit_length_arcsec / NIRSPEC_SPATIAL_PIXEL_SCALE_ARCSEC
        )
    result = detect_source_s2d(science, error_2d=error, dq_2d=dq,
                               plot=plot, outfile=outfile, **kwargs)
    result["slit_name"] = slit_name or None
    result["dq_used"] = dq is not None
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("infile")
    parser.add_argument("--plot", default="source_detection2_diagnostic.png")
    args = parser.parse_args()
    answer = wrap_detect_source_s2d(args.infile, plot=True, outfile=args.plot)
    if answer["status"] != "ok":
        print(f"result: failed ({answer['reason']})")
    else:
        print(f"result: {answer['measurement']:.3f} sigma")

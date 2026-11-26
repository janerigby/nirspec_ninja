# Fixed-slit spectra

These tools extract sky spectra and detect contaminating continuum sources in
NIRSpec fixed-slit data. The fixed slits are always open during MSA and IFU
observations, so their spectra can provide contemporaneous sky measurements.

## Detect sources

The detector fits spatial profiles in rectified S2D spectra and reports source
significance, fitted profiles, and reasons for failed fits.

```python
from nirspec_ninja.fixed_slit.source_detection import wrap_detect_source_s2d
result = wrap_detect_source_s2d("exposure_s2d.fits")
```

To save a diagnostic plot from the command line:

```sh
python -m nirspec_ninja.fixed_slit.source_detection exposure_s2d.fits --plot diagnostic.png
```

The optional `allowed_source_center_rows` argument is a Boolean mask of spatial
rows: True permits a fitted source center on that row. It restricts candidate
centers, not the pixels used in the fit. The default, None, searches all rows.

## Extract spectra

The extractor spatially median-combines fixed-slit S2D or CAL spectra,
combining exposures when present. It writes wavelength, surface brightness,
standard deviation, and MAD to CSV.

```python
from nirspec_ninja.fixed_slit.extraction import median_combine_level3_nirspecFS2
frame = median_combine_level3_nirspecFS2(
    "exposure_s2d.fits", "S200A1", "/path/to/output/", detect_traces=False)
```

Create the output directory before running; its string path must end with a
slash. `clip_nan` controls the margin excluded near edges or mostly-NaN rows.
`detect_traces=True` (the default) also runs source detection and records
successful detection significance as `trace_sig2` in the CSV header.

The multiprocessor notebook's `run_extract` worker explicitly enables detection.
CSV metadata use `##KEYWORD  VALUE`: `trace_sig2` is corrected significance in
sigma, and `non-flat` is spatial non-flatness in percent. The non-flatness measurement
uses the same common-wavelength profile, pixel/edge masks, and 99%-flux source
mask as the gradient plots, so the extracted spectra can be screened for strong
remaining gradients. `-1` overrides from diagnostic tables are not applied by
this notebook extractor. CAL inputs are assessed on their combined 2D image;
S2Ds use their original SCI/ERR/DQ and wavelength arrays.
CAL diagnostics retain the row exclusions already applied during extraction;
they do not expand the edge mask a second time, which would erase short S1600
cutouts. S2D diagnostics apply `clip_nan` once to the original pixels.

Diagnostics do not change the extracted spectrum. An unavailable score is
written as `##non-flat  nan`; statuses and reasons stay in worker reports,
not CSV header keywords.
A failed trace fit is reported without `trace_sig2`, while extraction continues.
If the validity/row mask leaves fewer than four usable rows, detection is
unavailable rather than a failed fit. If no finite wavelength/flux bins survive
extraction, both diagnostics are `nan` and the worker reports unusable data.
This can occur in spectra with
fragmented coverage: rows must have more than half their spectral pixels valid,
and the configured edge/neighbor exclusion margin must also leave usable rows.
`detect_traces=False` skips detection and source masking but still measures
non-flatness. The [provisional 8% review threshold](../diagnostics/README.md)
is not an automatic rejection rule.

See [background estimation and subtraction](../subtract_background/README.md)
for using fixed-slit measurements to estimate IFU sky backgrounds.

## Extract detected sources

The package PSF reference is the user-supplied
`data/prism_PSF_subsampled_x4.csv`, preserved unchanged. `psf.py` provides
`psf_fwhm_arcsec(wavelength_um, grating, filt)` and
`convolved_fwhm_arcsec(wavelength_um, intrinsic_fwhm_arcsec, grating, filt)`.
The latter assumes a Gaussian intrinsic source convolved with a Gaussian
PSF, adding FWHMs in quadrature. The table supplies widths, not full PSF images;
its author/generation method were not supplied. Interpolation is linear within
each configuration's tabulated range, with NaN outside coverage and no fallback
for configurations absent from the reference. These helpers establish the
wavelength-dependent width model; the existing collapsed-profile detector
still measures its original wavelength-collapsed width.

`source_extraction.py` extracts the strongest continuum trace in each calibrated
fixed-slit S2D, using the existing detector's corrected significance. The
notebook extraction cell offers `extract_sources=False` and
`source_significance_threshold=10.0`; selection requires strictly greater than
the threshold. Background extraction and its CSVs retain their existing behavior.
Source products go to
`/Users/jrrigby1/SCIENCE/JWST_Data/NIRSpec_FS/FS_sources/<epoch>/<target>/`.
The optional pass can also be run directly:

```python
from nirspec_ninja.fixed_slit.source_extraction import wrap_extract_sources
reports = wrap_extract_sources("/path/to/L3_destripe", "/path/to/sources",
                               threshold=10., clip_nan=5)
```

Each `<S2D-stem>_source_extract.csv` uses the background CSV's `# filename`
and `##KEYWORD  VALUE` headers, including its usual observation metadata.
Computed significance, trace center/width, and aperture bounds are truncated
to three decimal places in the header only. Extraction calculations, spectral
columns, copied FITS metadata (including `velosys` and `mjd`), and the configured
selection threshold retain their precision.
Columns start with `wave,fnu`, followed by `fnu_err,background,source_pixels,
background_pixels`. Added columns `fnu_extrapolated,fnu_extrapolated_err,
aperture_model_fraction` report a Gaussian extrapolation separately. `fnu`
always remains measured aperture flux and is the flux shown in the PNG.
The extrapolated error scales the measurement error only; it does not include
uncertainty in fitted geometry or PSF/source-model mismatch. No slit-loss
correction is applied. Wavelengths retain the S2D wavelength units (JWST microns).
`background` retains native SCI units. When SCI has surface-brightness units
and a positive `PIXAR_SR`, `fnu` and `fnu_err` are integrated Jy; otherwise they
are native pixel sums, identified by `BUNIT` and `flux_kind` in the header.
Reruns replace matching source CSVs; below-threshold files write nothing and
do not remove products from earlier runs.

Each extracted source also produces a matching `_source_extract.png` beside
its CSV. It shows the rectified background-subtracted image, fitted center, extraction aperture,
background rows and excluded rows; a spatial profile with propagated errors
over complete extraction columns; the background-subtracted spectrum with
1-sigma uncertainties; and the measured and smoothed local background.
The image uses native spectral columns to avoid misrepresenting the nonlinear
PRISM wavelength grid. The source-spectrum panel uses matching native column
positions with wavelength ticks in microns, so features align vertically with
the 2D image; its wavelength spacing is nonlinear. The background panel uses
a linear wavelength axis in microns. Coverage is reported in the summary text.
The spatial profile has spatial row on its vertical axis, shared with the 2D
image, so the fitted center, aperture, and background rows align horizontally.
For contrast, the 2D panel displays `(SCI − local background) / ERR` with a
symmetric asinh stretch spanning −5 to +5 and a linear width of 0.5. The image
is grayscale: positive signal is darker and negative signal lighter. This display uses the actual
extracted background; it does not change the extraction or CSV. These pixel
ratios exclude background uncertainty and resampling covariance and are not
the detector's corrected significance. A high collapsed-profile significance
can still correspond to a faint signal in individual spectral pixels.
Labels report corrected significance and threshold, fitted FWHM, aperture
bounds/clipping, valid-bin coverage, background-row count, and the detector's
number of spectral blocks above 3 sigma when available. Green background
rows are candidate rows; wavelength-dependent usable counts remain in the CSV.
Plots support human assessment, not a new detection criterion.
Plot failures are reported explicitly and preserve the extracted CSV.
The PNG also reports the established source-masked spatial non-flatness in
percent, saved as `##non-flat  VALUE` in the source CSV header (three decimal
places). Notebook workers reuse the matching background diagnostic; standalone
extraction computes it with the same routine. Unavailable measurements are
`nan` in the CSV and labeled unavailable in the PNG; reasons stay in reports.

The default boxcar aperture now encloses 84% of the fitted Gaussian profile with fractional
edge pixels. Its measured local background is an inverse-variance weighted mean of the
usable rows outside the aperture, independently at each wavelength. Before
subtraction, the shared LEGGOS grating continuum fitter smooths this background
along wavelength: median anchors suppress spikes, then a robust cubic spline
suppresses narrow fluctuations. It never extrapolates beyond measured coverage.
Gratings retain the LEGGOS default of 0.5 µm. Unlike the LEGGOS PRISM rule,
source extraction smooths the PRISM continuum too, because individual-slit
backgrounds have low signal-to-noise. PRISM preserves the measured background
exactly from 3.15–3.45 µm, protecting the real 3.3 µm feature. That interval is
excluded from the continuum fit. Cosine blends join measured background to
continuum across 3.10–3.15 and 3.45–3.50 µm. Missing measured feature bins remain
missing rather than being replaced by a continuum estimate. The selected
continuum scale is 0.5 µm, reviewed on four PRISM examples.
Only PRISM spectral background smoothing has been tuned in this work.
Grating smoothing has not yet been optimized for source extraction: its 0.5 µm
scale is inherited from the LEGGOS background workflow, not established here
as optimal for individual fixed slits. Grating features and noise suppression
need their own comparisons before treating this scale as validated.
`prism_smoothing_um` and `grating_smoothing_um` can be supplied to `run_extract`
or the standalone extractor. Use the comparison script described in
[scripts/README.md](../../../scripts/README.md) to review candidate PRISM scales.
The source CSV's `background` column is the curve actually subtracted;
`background_smoothing` and `background_smoothing_um` record the method and scale.
For PRISM, `background_unsmoothed_lo_um`, `background_unsmoothed_hi_um`, and
`background_blend_width_um` record the protected interval and edge blends.
The PNG overlays the measured background and the smoothed curve. At least
two background pixels and a complete usable source aperture are required per
bin; incomplete bins are NaN. DQ DO_NOT_USE, invalid SCI/ERR/wavelength, and
nonpositive WHT are excluded. The notebook's `clip_nan` row margin is retained.
Clipped apertures are recorded in `aperture_clipped`; no slit-loss or aperture
correction is applied. Uncertainties retain the unsmoothed background-estimation
term as a conservative approximation and assume independent pixels. They do
not propagate through the robust spline or model background-fit systematics,
smoothing-induced covariance, or resampling covariance.
Within the unsmoothed PRISM interval, the retained background-error term is
the direct measurement's statistical variance, so its noise remains included
in the source uncertainty under the independent-pixel assumption.
These initial extraction choices need review on real detections, especially
near slit edges or where another source contaminates the background rows.
Detection/skip/failure reasons remain in reports, outside CSV metadata.

For a Gaussian with independent, uniform background-dominated noise, boxcar
S/N peaks near a 1.4-sigma half-width (84% enclosed flux), about 15% higher S/N
than the former 99% aperture. This theoretical choice is not optimized for
each dataset's source Poisson noise or background-estimation uncertainty.
`enclosed_flux` is configurable in the extractor and `run_extract`; use 0.99
to reproduce the previous aperture. Extrapolation uses the fitted Gaussian's
integrated fraction in the actual fractionally weighted pixel rows. It remains
an approximate correction, especially for clipped/non-Gaussian sources.

Boxcar extraction is the initial choice because its aperture and pixel sums
are straightforward to inspect using the detector's fitted center and width.
This choice is not based on a comparison demonstrating better performance
than optimal extraction. Optimal extraction is worth considering for faint
traces: weighting by the source profile and variance could improve
signal-to-noise. It would require validating the spatial profile across
wavelength; the detector currently fits a wavelength-collapsed Gaussian,
which could provide a starting profile but has not been validated for that use.

See the [source-extraction development record and remaining work](../../../MIGRATION_STATUS.md#fixed-slit-source-extraction--2026-11-26)
for the reviewed PRISM comparisons, verification status, and future tasks.

## Galaxy redshift likelihoods

`redshift.py` fits measured source spectra with grizli's 12 FSPS no-line
continuum templates plus its fixed-ratio emission-line complexes. This is a
spectrum-derived redshift scan, not a fit to broadband photometry or grizli
beam products. Template coefficients are nonnegative; no additional fitted
background is introduced. The measured `fnu,fnu_err` columns in Jy are used,
not the extrapolated spectrum. Negative measurements remain in the fit.

The default screening requires `trace_sig2 > 10` and `non-flat <= 6%`.
Missing/nonfinite metrics are skipped. `source_threshold` is not used for
selection. The initial routine supports PRISM with local point-source and
uniform-illumination resolving-power fits, restricted to 0.6–5.3 µm. Other
gratings are reported skipped pending spectral-LSF calibration. These spectral
resolution fits are distinct from the spatial PSF reference. Along-slit source
width does not determine across-slit illumination.

The forward model convolves redshifted spectra with a wavelength-dependent
Gaussian LSF and averages into native bins. grizli's default intrinsic
emission-line width is 400 km/s. IGM absorption uses `eazy.igm.Inoue14`;
missing eazy causes an explicit failure, not silently absent absorption.
The initial grid spans z=0–12 in steps of 0.002 in log(1+z), configurable at
the command line. It is a grid search, without continuous redshift refinement.

```sh
python scripts/fit_source_redshifts.py /path/to/FS_sources
```

The runner finds the local grizli checkout under
`~/AI_experiments/References_and_Tools/grizli` (override with `--grizli-path`).
An unbuilt checkout (missing `grizli/version.py`) must first be installed in
the Python environment, letting grizli generate its version module and install
its dependencies:

```sh
python -m pip install -e ~/AI_experiments/References_and_Tools/grizli
```

Use your JWST environment with grizli and eazy available; neither is imported
by normal source extraction. `--max-nonflat`, `--min-significance`, `--zmin`,
`--zmax`, `--dlog1pz`, and `--resolution` configure the scan. `--output`
overrides the requested default:
`/Users/jrrigby1/SCIENCE/JWST_Data/NIRSpec_FS/FS_source/Redshift_likelyhood`.
Directory inputs retain their relative subdirectories in the output, avoiding
collisions; matching result files are replaced on rerun. Inputs are not edited.

Each eligible spectrum produces a `*_redshift_likelihood.csv`, model CSV,
JSON provenance (input/template hashes and settings), and diagnostic PNG with
measured spectrum/model, standardized residuals, and likelihood curves. The
likelihood CSV uses `##KEYWORD  VALUE` headers and contains z, χ², Δχ²,
relative likelihood, a normalized profile-likelihood density, and each spectral
resolution case's χ². When both cases are requested, the combined curve takes
the smaller χ² at each z, profiling over illumination rather than averaging.
There is no redshift prior or marginalization over template coefficients;
the normalized curve is not a calibrated posterior probability. Batch reports
include skipped inputs/reasons and failures separately.

Interpret a fit as a galaxy-template hypothesis, not confirmation that a trace
is a galaxy. Stars, Solar System objects, residual backgrounds and template
mismatch can yield misleading minima. Inspect residuals, competing peaks,
boundary solutions and χ² per used bin. Correlated resampling errors, fitted
background uncertainty and wavelength-dependent aperture losses are not modeled;
real-data/injection validation and template-family comparisons remain needed.

## Acknowledgments

The original source detection code was written by Rosalía O'Brien.

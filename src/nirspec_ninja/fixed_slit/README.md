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
CSV metadata use `# KEYWORD  VALUE`: `trace_sig2` is corrected significance in
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
written as `# non-flat  nan`; statuses and reasons stay in worker reports,
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

## Acknowledgments

The original source detection code was written by Rosalía O'Brien.

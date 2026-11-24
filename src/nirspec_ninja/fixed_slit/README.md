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

See [background estimation and subtraction](../subtract_background/README.md)
for using fixed-slit measurements to estimate IFU sky backgrounds.

## Acknowledgments

The original source detection code was written by Rosalía O'Brien.

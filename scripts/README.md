# Standalone extraction tests

In your JWST reduction Python environment, run:

```sh
python /Users/jrrigby1/AI_experiments/Coding/nirspec_ninja/scripts/test_source_extraction.py
```

Runs the new source-extraction tests plus existing background-extraction,
source-detection, and notebook-setting checks using this checkout's source.
No science data or reductions are needed. Separate suite logs, environment
versions, and `summary.json` are saved in a fresh run directory under
`~/AI_experiments/Temp/source_extraction_review/`. Each suite continues after
another fails; the script exits nonzero on any failure. `--output` selects a
new directory; `--timeout` sets seconds per suite (default 900).
Temporary synthetic FITS/CSV products are cleaned up by the tests; review logs
are retained. These checks do not establish extraction accuracy on real S2Ds.

# Compare PRISM source-background smoothing

Run `python scripts/compare_source_background_smoothing.py /path/to/old_source_extract.csv`
in the JWST environment. Supply a source CSV made **before** background smoothing
was added; that file's `background` column is the raw measurement. New smoothed
CSVs are rejected to prevent fitting an already fitted curve.
The script compares 0.1, 0.2, 0.3, 0.5, 0.8 and 1.0 µm spline scales and writes
a PNG, curve table, holdout prediction scores, and summary under
`~/AI_experiments/Temp/source_background_smoothing_review/`.
Scores are exploratory because neighboring pixels may have correlated noise.
Review broad spectral structure and repeat on other sources before adopting
a default. `--scales` selects a different set. No extraction is rerun.

# Test preservation of the PRISM 3.3 µm background feature

Run `python scripts/test_prism_background_feature.py` in the JWST environment.
It reads the four example S2Ds directly from `~/AI_experiments/Temp/`, using
production extraction masks to recover raw local backgrounds and their errors.
`--input-dir` selects another folder. Results go in a new directory under
`~/AI_experiments/Temp/prism_feature_review/`: each example has a full-range,
feature-zoom and red-end comparison PNG, raw/candidate curves, and peak/area
metrics. A JSON summary records successes and failures.

The experimental candidates exclude 3.15–3.45 µm from a 0.5 µm continuum fit,
then restore an inverse-variance Gaussian-smoothed feature residual with
0.03, 0.05 or 0.1 µm kernel FWHM, tapering continuously to the continuum.
An additional unsmoothed candidate retains the measured samples exactly at
3.15–3.45 µm, blending to the continuum across 3.10–3.15 and 3.45–3.50 µm.
It is shown as a purple dashed line. Its 100% feature retention is guaranteed
by construction; inspect the retained noise and transitions as well as the peak.
These local kernel FWHMs differ from the spline's smoothing-scale definition.
Peak/area measurements use 3.2–3.4 µm relative to the feature-excluded continuum;
the raw peak is noise-biased and continuum uncertainty is not modeled.
Review feature preservation and red-end behavior before choosing a candidate.
Existing science products and production smoothing defaults are not changed.

# Fit source redshift likelihoods

Run `python scripts/fit_source_redshifts.py /path/to/FS_sources` in the JWST
environment with grizli/eazy available. It uses the local grizli checkout under
`~/AI_experiments/References_and_Tools/grizli`. The default cuts are measured
`trace_sig2 > 10` and `non-flat <= 6%`; missing metrics are skipped.
Outputs default to
`/Users/jrrigby1/SCIENCE/JWST_Data/NIRSpec_FS/FS_source/Redshift_likelyhood`.
Use `--output` to write a sandbox-accessible review folder in `Temp/`.
See the [redshift documentation](../src/nirspec_ninja/fixed_slit/README.md#galaxy-redshift-likelihoods)
for grid/resolution settings, products and scientific limitations. Initial
support is PRISM only; grating fits require a supplied spectral-resolution
calibration. The extraction test runner now also runs offline redshift
screening, spectral-response and synthetic recovery tests.

# Live-service smoke tests

Check that your reduction environment can retrieve WISE images and download
public JWST files from both AWS and MAST. Run this script from the package root
in your installed JWST reduction Python environment:

```sh
python scripts/test_live_services.py /path/to/public_nirspec_exposure.fits
```

Supply a local public NIRSpec exposure with pointing metadata. The script:

1. Queries IRSA for WISE W1, downloads an image, and creates a finder PNG,
   DS9 regions and provenance record using the package finder tool.
2. Downloads the corresponding public RATE file once from AWS and once from
   MAST, with one worker each, using the filename-list downloader.
3. Checks file sizes and FITS structure, then compares SHA-256 hashes of both
   archive copies. Each test records pass/fail separately.

Outputs and logs go in a new nirspec_live_services_* directory under
AI_experiments/Temp. summary.json contains results. No existing reduction is
rerun or overwritten. The script does not test complete program-query workflows
or scientific accuracy of the finder footprints.

--archive-file selects a different public JWST FITS basename. By default an
UNCAL input selects its corresponding RATE to reduce download size; another
input uses its own basename. --max-download-mb defaults to 512 MB per archive
copy and checks advertised size before downloading. --timeout defaults to 900
seconds per chart/download command. --output must name a new directory.

The WISE cache is fresh for each test run so a prior image does not bypass
retrieval. Logs are retained when a service fails; the other tests still run.
Open wise_finder.png to check the chart visually after the script finishes.

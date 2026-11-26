# Glow and column artifact detection

These tools flag broad glow from possible MSA electrical shorts and coherent
column artifacts from problems with 1/f correction. Detection alone does not
establish the cause. Supported inputs are full-frame PRISM/NRS1 IFU or MOS exposures.

## Score exposures

```sh
python -m nirspec_ninja.artifacts.detect score /path/to/inputs --mask-dir /path/to/masks --output /path/to/results
```

Directory scans include `*_rate.fits`, `*_rate_after_cleaning.fits`, and
`*_clean_flicker_noise.fits`. All products for an exposure use the same mask:
`exposure_rate_after_cleaning.fits` matches `exposure_mask.fits`.
Use `--PNG` for per-exposure diagnostic plots.

Glow measurements use illumination masks. A missing mask leaves the glow
measurement unavailable and marks the exposure for review.

## Interpret column scores

`column_score` is the maximum coherent column amplitude divided by the
exposure's empirical noise scale. Columns must have same-sign support in at
least four of sixteen spatial bands. Higher scores indicate stronger column
structure relative to noise; they are not formal significances or probabilities.
Compare exposures with similar observing modes and reduction settings.

No clean-data reference is required. `run.json` records scoring version 0.3
and the measurement method. Earlier reference-normalized scores and thresholds
are incompatible with these scores; rescore older exposures before comparing.

## Plot scores

```sh
python -m nirspec_ninja.artifacts.plot_score_comparison catscores.csv --show-histograms --hide-mystery
```

The plot uses the supplied CSV. `--training-csv` adds labeled samples; those
samples must have been scored with the current method.

[Historical reference scores](baselines/README.md) and [the earlier method](METHOD.md)
are retained for provenance, not for interpreting current column scores.

## Notebook integration

Spec2 scores each supported full-frame cleaning result before fixed-slit cropping.
The notebook's `save_artifact_png = True` saves
`<exposure>_rate_after_cleaning_artifacts.png` in L2a alongside the original and
cleaned RATEs. Rendering reuses the cleaned arrays. Set this option to False to
save scores without PNGs.

Cleaned RATE and Stage 2 CAL/S2D primary headers carry `glow_score` and
`column_score`. Stage 3 CAL/S2D headers and extracted CSV metadata carry
`glow_score_med`, `glow_score_max`, `column_score_med`, and `column_score_max`:
medians and maxima of available science-input exposure scores for that detector,
filter, and grating. FITS uses `gloamed`, `gloamax`, `colamed`, and `colamax`, respectively,
with comments defining median/max glow/column artifact scores. CSV metadata
retains the full names above. These are
input-exposure statistics, not measurements on the combined spectrum. Unsupported
or failed measurements are unavailable, never zero; CSV metadata uses `nan`.
Adjacent `.artifacts.json` records retain availability and Stage 3 coverage counts.
The existing CSV metadata markers and spectral table columns are preserved.

For an existing cleaned RATE, `nirspec_ninja.artifacts.pipeline.backfill_cleaned(path)`
saves scores and a PNG without repeating calibration. To update already-existing
Stage 2 or Stage 3 products, use `annotate` with the returned exposure record and
`propagate_stage3(asn_file, output_dir)`, respectively, then rerun CSV extraction.

Short regression test (from the repository root, in the JWST Python environment):

```sh
PYTHONPATH=src python -m pytest tests/test_artifact_pipeline.py -q
```

It checks file/array score equivalence, actual default PNG creation without
mutating SCI, and Stage 3 median/max propagation with detector-specific membership.

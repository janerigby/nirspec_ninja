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

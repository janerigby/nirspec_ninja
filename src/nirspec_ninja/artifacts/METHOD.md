> Historical v0.2 documentation. Clean-reference normalization and its thresholds
> are obsolete. See README.md for current v0.3 scoring.

# PRISM/NRS1 short-artifact screening — initial development version

This fresh implementation scores each exposure independently, using a reference fitted to the supplied human-verified clean files. It produces separate glow and column scores and a four-panel diagnostic PNG for every supported input. Scores measure morphological departures from the clean reference; they are **not probabilities, formal significances, or proof of an electrical short**. No rejection threshold is enabled by default.

## Run

From the project directory, use the supplied Python environment:

```sh
PYTHON=/Users/jrrigby1/AI_experiments/References_and_Tools/JWSTDP-3.0.0-10-py314-macos-arm64/bin/python
"$PYTHON" fresh_v2/detect.py score Inputs_v2 --output fresh_v2/results --save-maps
```

Pass individual FITS paths or directories (searched recursively). The default reference is `fresh_v2/reference.npz`. For a separately supplied validation directory, change the input and output paths; **reuse this reference without refitting on validation files**.

Optional decisions, once thresholds have been selected:

```sh
"$PYTHON" fresh_v2/detect.py score VALIDATION_DIRECTORY --output VALIDATION_OUTPUT --glow-threshold GLOW_VALUE --column-threshold COLUMN_VALUE
```

If either supplied threshold is met, `decision` is `reject`. With no thresholds, decisions are `unthresholded`. An input that cannot be scored receives an error row and `review`, and the command returns a nonzero exit status. If only one threshold is supplied, only that artifact score controls the decision.

To reproduce the clean reference:

```sh
"$PYTHON" fresh_v2/detect.py fit --clean Inputs_v2/IFU_okay Inputs_v2/MSA_okay --output fresh_v2/reference.npz
```

Fitting requires at least three clean files per mode. The current reference uses all 73 supplied clean files. Bad labels and notes are not used by the scoring algorithm. File paths recorded in the reference provide provenance; scoring does not infer labels from directory names. The CSV `reference_member` field identifies files used in fitting.

## Scores and localization

Only full-frame 2048×2048 PRISM/NRS1 exposures in IFU or MultiObject mode are supported. `TEMPLATE` determines mode. `S_CLNFNS` records clean-flicker-noise status; scoring does not require `COMPLETE`. Files with other dispersers, detectors, modes, or sizes receive error rows rather than a keep decision.

Nonfinite and DO_NOT_USE SCI pixels and the four-pixel detector border are excluded. Normalization uses the median absolute adjacent-column difference divided by 0.953872. This is a robust empirical pixel-noise estimate; artifacts can inflate it.

**Columns:** Take the 20th percentile within each of 16 disjoint 128-row bands. Subtract a 101-column running median. The fourth-largest same-sign band residual gives each column's coherence statistic, requiring support in at least four bands. Both positive and negative residuals are considered. This suppresses narrow horizontal traces, but does not guarantee separation from every science structure or glow edge.

**Glow:** Take medians in 8×8 tiles, subtract the exposure's global column median, apply a seven-tile vertical median, remove each tile-row's median, and smooth with a Gaussian of sigma two tiles. This suppresses isolated sources, thin spectra, and constant-height columns while retaining some broad structures. The difference from the clean reference is scored in both signs. The nominal tile centers are half-integer FITS/DS9 coordinates (first tile center: 4.5,4.5). Glow location is approximate, with eight-pixel sampling and additional filtering.

Each mode has a separate reference center (median clean feature) and scale (95th percentile of absolute clean deviations at each location). Minimum scales are 0.1 for columns and 0.5 for glow, in empirical noise units, to avoid unstable divisions. Each final score is the maximum location score. A score of 1 is a local comparison scale, **not a file-level percentile**; maxima across many locations commonly exceed 1 in clean files.

`--candidate-level` defaults to 5 and controls diagnostic overlays only; it is not a rejection threshold. Cyan lines show candidate columns and orange contours show broad candidates. Every image contains the raw SCI percentile stretch, a signed asinh SCI display, a glow candidate map, and a column profile. `--save-maps` saves the full score arrays in NPZ files. `scores.csv` contains peak coordinates, scores, processing history, diagnostic paths, and decisions. `run.json` records the invocation settings.

## Validation limits

The supplied clean files define the baseline, so evaluating those same files is a development comparison, not an independent false-positive estimate. Many inputs share observations. The supplied bad examples can show whether known affected files stand out; the notes are non-exhaustive and cannot establish separate glow/column ground truth or region-level completeness.

The method may confuse exceptionally bright or extended science emission with glow, and may miss narrow/faint glow, broad uniform offsets, or columns that occupy fewer than four bands. A column varying strongly with row can also enter the glow score. Extreme SCI values in some bad files produce very large scores; they are retained rather than silently clipped. The visual panels use bounded stretches for readability.

All supplied files except one have `S_CLNFNS=COMPLETE`; `IFU_bad/jw01247318001_03101_00002_nrs1_rate.fits` has `FAILED`. Successful operation without the step is supported by the code path, but detection accuracy on clean exposures reduced without it has **not** been validated. The separate validation set should include that processing variant before claiming coverage.

Synthetic behavior checks are in `check_injections.py`. They verify positive/negative columns and broad patches, coordinate localization, TEMPLATE precedence, and unsupported-grating rejection. They do not establish accuracy on actual short artifacts.

## v0.2: geometry-masked glow (2026-11-05)

`detect.py score` now uses the supplied illumination masks for glow. The previous
CLI is preserved as `detect_v01.py`; its reference and all previous results remain
unchanged. Column extraction and reference normalization are unchanged. Glow
scores in v0.2 have different units and must not be compared to v0.1 glow scores.

```sh
"$PYTHON" fresh_v2/detect.py score Inputs_v2 \
  --mask-dir Inputs_v2/Example_masks --output fresh_v2/new_masked_results --save-maps
"$PYTHON" fresh_v2/check_masked_glow.py
```

Directory scans select `*_rate.fits` and `*_clean_flicker_noise.fits`, so masks are not treated as exposures.
Masks match by exact exposure basename, replacing either product suffix with `_mask.fits`, and are
validated for shape, binary values and matching exposure/geometry metadata.
SCI mask 1 means background; 0 means excluded spectral footprint. Missing masks
produce a column score, blank glow score, `glow_status=missing_mask`, and `review`
unless the column threshold already requires rejection. Invalid masks produce
error rows. Neither condition permits an unsupported keep decision.

Glow uses background-only 8x8 tile medians with at least 16 valid background
pixels per tile. No spectral pixels enter these medians, and there is no smoothing
across mask boundaries. The background baseline is its 20th-percentile SCI value;
there is no row/column background subtraction. Brightness must exceed the
90th-percentile SCI value in valid excluded footprints (both measured relative
to that baseline), with a minimum contrast of one empirical noise unit.
`--spectrum-percentile` changes this initial operational definition of spectral
brightness. It is a practical proxy, not a fitted spectrum or expert-calibrated
brightness measure; excluded footprints can contain empty sky or incidental
features.

Edge-connected eligible tiles form patches. A patch must contain bright
background pixels in at least 30 distinct detector rows (`--min-glow-rows`).
The row count is measured on original pixels, not rounded to tile height.
Tile sampling can nevertheless miss thin or boundary-aligned patches, and
masked footprints can split a single physical glow region into separate patches.
Only positive glow is targeted; positive and negative column detection remains.
Glow score is the maximum retained tile contrast divided by empirical pixel
noise. No retained patch means score zero and blank peak coordinates. CSV fields
record mask path, baseline, spectral brightness, component count and maximum
component row count. Green diagnostic boundaries show excluded footprints;
orange contours show candidates above the display level. Noise normalization
can be inflated by artifacts. Near-uniform detector glow can raise the global
baseline and be suppressed.

There are 21 exact mask/rate matches in the supplied data: 5 known clean files,
1 known bad file annotated with column artifacts, and 15 unlabeled mystery files.
`masked_results_v02/` contains these 21 plus a missing-mask example: 22 CSV rows,
PNGs and maps, zero processing errors. None of these 21 passes the current glow
scope/brightness criteria; this is not evidence that genuine extended glow is
being recovered. Known glow-positive development files do not yet have matching
masks. The existing clean reference is used only for columns; there is insufficient
masked clean coverage to rebuild a mode-specific glow reference (only one IFU
clean match). Synthetic checks verify 29/30-row behavior, faint-glow exclusion,
spectral-only exclusion and broad bright-glow recovery, not real-data accuracy.

Mask directory lookup now searches recursively, supporting the nested layout in
`Inputs_v2/Many_masks`. Duplicate matching basenames are rejected as ambiguous.
To score the expanded inventory:

```sh
"$PYTHON" fresh_v2/detect.py score Inputs_v2 --mask-dir Inputs_v2/Many_masks \
  --output fresh_v2/many_masks_expanded_results_v02 --save-maps
"$PYTHON" fresh_v2/summarize_masked_run.py fresh_v2/many_masks_expanded_results_v02
```

Expanded results are now consolidated in `many_masks_expanded_results_v02/`:
146 successfully processed exposures, 145 with masks. Identical duplicate mask
SCI arrays and checked metadata are accepted; conflicting copies fail. `run.json`
links the supplemental runs used to replace initially missing or duplicate-error
rows. `report.md` and `noted_glow_regions.csv` include known-location checks.
All 73 clean examples have zero retained glow; 21/32 masked bad examples and
1/40 unlabeled mystery files have nonzero glow. Mixed artifact labels do not
establish glow-specific accuracy. Some noted glow regions are missed, and bright
columns can enter the glow map. The sole missing mask is
`jw01373020001_03101_00001_nrs1_mask.fits`.

## External files and comparison plots

See `RUN_OUTSIDE.md` for complete commands. `score --file-list` reads absolute
RATE paths, one per line. Without `--mask-dir`, masks are located adjacent to each
RATE under `custom_masks/<exposure-stem>_*/mask.fits` (or by exact named mask
basename). Use explicit `--mask-dir Inputs_v2/Many_masks` for the local inventory.
`plot_masked_scores.py BASELINE_CSV --new-csv NEW_CSV --output PLOT.png` overlays
new blue squares on Dev Bad, Dev OK and Mystery, with matching main-panel ranges
and enlarged low-score panels for readable OK points.

`score --PNG False` now disables per-file diagnostic PNGs while still generating
the final `score_comparison.png`. Supply `--plot-baseline BASELINE_CSV` when the
baseline is outside its default location next to this checkout. PNG defaults to
True. CSV and optional NPZ maps are still generated. Scoring skips the unused
legacy glow feature; the current masked glow score is unchanged.

Default CSV reporting is compact (eight columns): file, mode, column_score,
glow_score, status, glow_status, decision, error. `--detailed-csv` optionally writes
all prior metadata to a separate scores_detailed.csv. Summary figures now include
marginal x/y count histograms by group on each main panel, shared score/count
ranges and bins, and the existing enlarged low-score panels.

## Standalone editable comparison plot

`plot_score_comparison.py` reads CSV files using Matplotlib and NumPy. Edit the
settings near the top for colors, markers, figure size, title, zoom ranges and
resolution; edit `plot_scores()` for layout changes. No detector or FITS inputs
are needed. The existing `plot_masked_scores.py` command forwards to this script.

```sh
python fresh_v2/plot_score_comparison.py baselines/reference_scores_v02.csv --output comparison.png
# Optional overlay:
python fresh_v2/plot_score_comparison.py baseline/scores.csv --new-csv new_run/scores.csv --output comparison.png
```

Without `--output`, it writes `score_comparison.png` next to the baseline CSV.
Rows with failed status, missing masks, or missing glow scores are omitted.

To compare a run directly against the bad and okay training samples:

```sh
python fresh_v2/plot_score_comparison.py Inputs_v2/scores.csv --training-csv baselines/reference_scores_v02.csv --output comparison.png
```

Training files under `IFU_bad`/`MSA_bad` are tomato triangles and those under
`IFU_okay`/`MSA_okay` are teal open circles. The run is shown as blue squares.
Use the masked training CSV above for comparison with current masked glow scores.

The standalone command defaults `--training-csv` to
`../baselines/reference_scores_v02.csv` relative to the script, independent of
the working directory. Supply your results CSV alone to include that training
comparison; pass `--training-csv` to choose a different training file. CSV inputs
are read only.

Comparison scatter markers use hollow outlines and different sizes so coincident
bad-training triangles and new-run squares remain visible. Edit `MARKER_SIZES`,
`MARKER_LINEWIDTH`, and `MARKER_ALPHA` near the top of the standalone script.
Score coordinates and histogram counts are unchanged.

Mystery samples are shown by default. Add `--hide-mystery` to exclude them from
all scatter panels, histograms, and automatic axis ranges. `--show-mystery`
explicitly enables them. This toggle does not hide the run labeled New.

Histograms now use 64 bins and larger marginal panels. Each count axis has its
own limits and uses a symlog scale (linear below 5 counts, logarithmic above),
so low-count bins remain readable alongside the concentrated okay population.
Score bin edges remain identical across IFU and MSA. Customize `HISTOGRAM_BINS`,
`HISTOGRAM_PANEL_SIZE`, and `HISTOGRAM_COUNT_LINTHRESH` near the script top.

Bad and okay training histograms have opposite diagonal hatches with transparent
interiors; New and Mystery remain plain outlines. Each top histogram has a legend.
Edit `HISTOGRAM_HATCHES` to change the training patterns.

Plot legends label the populations `Visual Bad`, `Mystery`, `Input`, and
`Visual OK`. Use `--input-label "My sample"` to rename Input in both scatter
and histogram legends. This label applies to the positional results CSV and
any additional CSVs supplied through `--new-csv`.

Marginal histograms are off by default. Use `--show-histograms` to enable them
or `--hide-histograms` to explicitly disable them. Scatter panels and reference
lines are shown with either setting.

## Cleanup — 2026-11-12

The frozen plotting baseline is now `baselines/reference_scores_v02.csv` at the
project root. Both scoring and standalone plotting use it by default. The original
consolidated run stays in place. Historical run directories and logs are under
`archive/cleanup_2026-11-12/fresh_v2/`; earlier commands mentioning moved run
folders should use those archived locations. See the project-root README.

The local development archive was subsequently deleted after confirmation of
an external backup. Historical helpers need those products restored to their
archived paths before use. Cleanup manifests remain in
`cleanup_records/2026-11-12/`; active scoring and plotting use the preserved
baseline and do not require historical run products.

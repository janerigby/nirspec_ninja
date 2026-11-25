# Restart summary — 2026-11-23

Before editing, read [USER_PREFERENCES.md](USER_PREFERENCES.md), especially the requirement to check before changing downstream file formats.

## Current state

Main nirspec_ninja migration is complete. Root:
/Users/jrrigby1/AI_experiments/Coding/nirspec_ninja
Application code: src/nirspec_ninja. README.md, INSTALL.md and MIGRATION_STATUS.md
contain documentation/history. Package installed editable in reference interpreter:
/Users/jrrigby1/AI_experiments/References_and_Tools/JWSTDP-3.0.0-10-py314-macos-arm64/bin/python
User also runs their own patched JWST 3.0.0 environment; do not assume identical.

```text
nirspec_ninja/
├── masks/                              ✓ Migrated
├── pipeline/
│   ├── patches/                        ✓ Migrated
│   ├── notebook_workers.py             ✓ Migrated
│   └── execution/multiprocessing/etc.   ✓ Shared helpers
├── fixed_slit/                          ✓ Detection + extraction
├── subtract_background/                ✓ Separate from fixed_slit
├── artifacts/                          ✓ Glow + raw column score
├── diagnostics/                        ✓ Noise/gradients + dataset scripts
├── finder/                             ✓ Finder chart
└── archive/                            ✓ Queries/downloads/retries + MSA metadata
```

LEGGOS remains its own readable script, notebook its own readable staged workflow.
Cleaning stays in Detector1 for LEGGOS and Spec2 for notebook. Do not merge them.

## User rules

- ASK before renaming functions/APIs or changing scientific behavior.
- ASK before ANY notebook edit. Preserve comments, formatting, visible settings.
- Do not read/import/migrate Scratch code.
- Audience: astronomy grad students, familiar with regression tests. Plain prose;
  avoid unexplained software jargon such as upstream.
- Git/licensing/attribution are LAST. Git previously failed due missing Apple tools.
- User accepted LEGGOS overnight validation; no further analysis/rerun needed.
- Remind user masks need real-data checks for other grating/filter combos on BOTH
  detectors. Only one G235M/F170LP exposure has predicted-aperture evidence.
  PRISM need for extra masking is unvalidated. See src/nirspec_ninja/masks/VALIDATION_COVERAGE.md.
  NRS2 IFU-mask integration issue remains described in MIGRATION_STATUS.md.

## Recent approved behavior

- Spec2 notebook worker always installs/verifies combined runtime patch with cap
  4.0. No installed JWST files edited. Mask companion *.flicker_patch.json written.
- Routine Spec2 success messages suppressed; failures/tracebacks visible.
- Spec2 writes L2a/flicker_noise_rms.txt comparing original RATE with retained
  *_rate_after_cleaning.fits, using custom mask + actual fit_method.
- Noise report: short headings, aligned TWO-SPACE separators (no expanded tabs).
  Noise 4 significant digits, percent reduction 2 decimals, timing 6 significant
  digits. RMS/MAD reductions, readout pattern/time/groups/integrations/frame counts,
  detector/subarray included. Process-safe updates; old reports readable/upgraded.
- Artifact scoring v0.3: column_score is raw exposure-noise-normalized coherent
  column amplitude, no clean reference. USER DECISION: no further normalization,
  validation or cleanup planned. Glow unchanged. Old reference scores incompatible.
- Plotter now uses ONLY positional CSV as Input by default; Visual Bad/OK require
  explicit --training-csv. Never silently add historical comparison samples.

## Compatibility cleanup completed

User updated notebook imports (confirmed saved notebook):
from nirspec_ninja.pipeline.multiprocessing import max_filesize
from nirspec_ninja.archive.download_msa_metadata import download_metadata
User removed gethead/mark_CaII_Fraunhofer_lines imports.
max_filesize and download_metadata migrated with names/behavior unchanged.

Updated 34 known Python caller/test files in Coding. Retired six unused import-only
redirects: original mask project's custom_nirspec_mask_2/aperture_extrapolation/
nirspec_order_geometry, JWST_pipeline/fs_background.py, MAST_queries/download_metrics.py,
nirspec_ninja/jrmulti.py. CLI wrappers remain deliberately for existing commands.
Removed supplied jrr copy; user said original exists elsewhere, no backup needed.
Backup of jrr removed. Other edited-file/retired-wrapper records remain in
Temp/migration_cleanup_20261123/manifest.json. Unrelated historical
Detect_trace/jrjwst.py still depends on jrr but is not the active notebook extractor.

## Tests

After cleanup: all 64 package tests and 44 original archive offline tests passed
in JWST 3.0.0. Seven worker integration and five noise-report tests passed after
recent Spec2 changes. Synthetic/model tests are not full scientific validation.
Notebook-stage test fixtures were updated to supply summarize_extract, matching
current notebook; notebook not edited. Logs in Temp/migration_cleanup_20261123/.

## Live-service checks complete — 2026-11-23

User ran scripts/test_live_services.py. Saved summary verified at:
/Users/jrrigby1/AI_experiments/Temp/nirspec_live_services_h6cfhk_e/summary.json
WISE retrieval/chart generation, AWS download, MAST download and byte comparison
all passed. Both RATE copies are 83,969,280 bytes and have identical SHA-256.
Input: jw04125018001_02101_00001_nrs1_rate.fits from Temp/COSMIC.EYE/L2a/sci.
These smoke checks are complete; do not list them as pending or rerun them.
Broader mask validation and repository setup remain. Visual chart review has
not been reported; no claim of scientifically validated aperture placement.

# Migration status and restart notes

## Current status — 2026-11-23 (caller cleanup)

- Functional migration is complete. Both workflows retain their own structure.
- User updated notebook imports for max_filesize and download_metadata; verified
  in saved notebook. No notebook edit made by the agent in this cleanup.
- Known Python callers in Coding updated to package imports. Six unused
  import-only redirects retired; public CLI entry points retained.
- Supplied copied jrr directory removed; user confirmed original exists elsewhere
  and requested no backup of that copy. Unrelated historical Detect_trace/jrjwst.py
  remains an old standalone copy; it is not the active notebook extractor.
- LEGGOS overnight results accepted by user; no further analysis required.
- Raw noise-normalized column score is the chosen measure. Further column-score
  normalization/validation/cleanup is not planned. Plotting old reference scores
  against these is inappropriate; old samples are no longer loaded by default.
- Git/licensing/attribution deliberately left until last.
- REMINDER: mask-making needs real-data testing on other grating/filter pairs,
  separately for NRS1/NRS2, including PRISM. Only one G235M/F170LP exposure has
  recorded predicted-aperture comparison evidence. NRS2 integration issue remains
  separately recorded below. Live AWS/MAST and WISE/IRSA paths were not tested
  during migration; offline tests were.

## Working rules and decisions

User wants one maintainable package with shared implementations and small task
modules. Preserve public function names and behavior; ASK before renaming or
changing APIs. User explicitly objected to unapproved renames. Do not read/import
or migrate code from Scratch; timestamp inspection there was permitted for
comparison only. Development/testing uses JWST 3.0.0 in:
/Users/jrrigby1/AI_experiments/References_and_Tools/JWSTDP-3.0.0-10-py314-macos-arm64/bin/python

Package installed editable in that environment. INSTALL.md contains instructions.
Existing callers use documented compatibility redirects (introduced 2026-11-19,
remove after all callers move to package imports). The supplied jrr, compatibility jrmulti.py, and active notebook are in the
package root; only src/nirspec_ninja is the application package. Source detector copies in supplied jrr and
Detect_trace were byte-identical. Supplied jrmulti/notebook timestamps were
newer than corresponding versions found elsewhere; user requested stopping if
any newer version is found.

Documentation audience: astronomy graduate students. Assume astronomy jargon
and knowledge of regression tests, but avoid unexplained software jargon.
Preferred opening tone: "These patches (or fixes) address bugs with the Clean
Flicker Noise function in the JWST pipeline as of v3.0.0. Importing this package
does NOT apply the patches. The patches are applied by running the commands
in your JWST reduction Python environment."

## Revised pipeline/workflow plan — agreed 2026-11-20

User revised the earlier proposal to preserve two independently readable
workflows, rather than merge them into a common workflow format.

- Retain the multiprocessor notebook's cells, formatting, and staged execution,
  and retain the LEGGOS Python script's structure.
- First migrate existing implementations and replace imports/necessary calls
  only. Expose settings currently hidden in jrmulti in a separate later edit.
  Do not add unrelated functionality or reformat the notebook.
- Keep explicit parameter dictionaries beside stage calls in each workflow.
  Shared parameter helpers apply those choices, rather than choose scientific
  settings behind the scenes.
- Separate processor scheduling and CRDS preparation in
  pipeline/multiprocessing from functions that actually run pipeline stages.
  Preserve separate stage runners where the workflows differ.
- Preserve LEGGOS clean_flicker_noise in Detector1 and notebook cleaning in
  Spec2, including their respective masks, settings, and diagnostic products.
- Essentially all operational jrmulti code can move into importable package
  modules, which can serve as spawned notebook workers. Keep jrmulti as a
  temporary compatibility redirect; preserve public names and signatures.
- Move existing implementations first, then consolidate genuinely shared
  helpers incrementally. Verify parameter handling, cleaning placement,
  outputs, and multiprocessing, including a real notebook-kernel check.

This replaces the earlier consolidation proposal. No pipeline or notebook
implementation has been changed as part of this planning discussion.

## Migration map

```text
nirspec_ninja/
├── masks/                              ✓ Migrated
├── pipeline/
│   ├── patches/                        ✓ Migrated
│   └── shared execution/multiprocessing ✓ First-pass helpers migrated
├── fixed_slit/
│   ├── source_detection.py             ✓ Migrated
│   └── extraction.py                   ✓ Migrated
├── subtract_background/                ✓ Migrated
├── artifacts/
│   ├── glow and column detection       ✓ Migrated
│   └── column-score validation/cleanup  ○ Deferred
├── diagnostics/
│   ├── noise and fixed-slit gradients  ✓ Migrated
│   └── column/stripe investigations    ✓ Migrated
├── finder/                             ✓ Migrated
├── archive/
│   ├── three query/download tools      ✓ Migrated
│   └── shared transfer/retry code      ✓ Consolidated
└── workflows/
    ├── LEGGOS script                   ✓ Shared helper imports; own runners retained
    └── FS notebook + jrmulti           ✓ Workers migrated; notebook imports updated
```

Map describes subdirectories of src/nirspec_ninja. DO NOT migrate obsolete
FS/FS_Serendip/screen_review_artifacts.py. Background subtraction belongs in
subtract_background, NOT fixed_slit; FS sky measurement is one candidate method.

## Implementations and approvals

- Masks: custom_nirspec_mask_2.py, aperture_extrapolation.py,
  nirspec_order_geometry.py. Original names restored after user objected to
  renames. make_custom_NIRSpec_mask is public; make_background_mask was removed.
  Algorithms/defaults preserved, including old CRDS environment defaults.
- Patches: jwst_3.0.0_saturated_groups.patch + PATCH.md + APPLY_PATCH_tcsh.md;
  patch_clean_flicker_noise.py; nirspec_flicker_patch.py. Patch files/scripts
  preserved; docs point to new location. Import applies nothing. Runtime patch
  requires explicit install_patch per worker, supports tested JWST 3.0.0 only,
  cap default 4.0 experimental. No installed JWST source was patched.
- Source detection: fixed_slit/source_detection.py filename approved; original
  functions preserved. Added approved allowed_source_center_rows Boolean mask
  (True allows candidate center on row; None preserves original search).
- Extraction: median_combine_level3_nirspecFS2 and related helpers migrated
  from supplied jrr.jrjwst. User directed removing old detector calls.
  Only new detector used; old trace_sig header removed, trace_sig2 retained.
  Approved detect_traces=True option; False skips detection for LEGGOS.
  LEGGOS imports extractor directly. Supplied jrr.jrjwst redirects moved funcs.
- subtract_background/fs_background.py imports migrated detector and uses
  allowed_source_center_rows. subtraction.py contains original exp_bg_sub.
  LEGGOS retains target/file choice; imports shared subtraction function.
- artifacts: detect.py, features.py, masked_glow.py, plot_score_comparison.py,
  plot_masked_scores.py. Baseline assets retained for historical use. User
  declared clean reference obsolete: current v0.3 scoring does NOT load it.
  column_score is maximum existing noise-normalized coherent column amplitude.
  Glow calculation unchanged. Automatic old-baseline plotting removed.
  --reference/--plot-baseline accepted but warn and ignored; fit errors as
  obsolete. User said good enough, fix later; issue documented in artifacts
  README. METHOD.md is historical v0.2. Old scores/thresholds not equivalent.
- Diagnostics: reusable fs_gradients.py, flicker_noise_diagnostics.py plus
  inspect_fs_gradients.py. Dataset-specific inspect_columns.py, inspect_batch.py,
  compare_ramps.py retain original dataset/output paths and execute only as
  scripts. Those data investigations were not rerun.
- Finder: finder/nirspec_finder.py copied unchanged; CLI/old redirect supported.
  Default cache is now beside moved module; docs describe it. Offline local
  image chart verified; online WISE retrieval not tested.
- Archive: original three MAST_* filenames preserved under archive; old
  locations redirect. download_metrics shared. common.py shares batching,
  source configuration, atomic JSON. transfers.py shares retry execution,
  timing, callbacks for success/failure and Astroquery transfer calls.
  Retry counts/delays/error messages/restart policies preserved. PID still
  stops after one failure. HTTP resumption vs Astroquery staging distinct.

## Latest task: flicker_noise_rms.txt

Diagnostic changed in src/nirspec_ninja/diagnostics/flicker_noise_diagnostics.py.
Original noise calculations retained. Adds actual fit method (LEGGOS passes
clean step fit_method), MAD percent reduction, and EFFEXPTM, READPATT, NGROUPS,
NINTS, NFRAMES, GROUPGAP, TGROUP, TFRAME, DETECTOR, SUBARRAY. Missing values unknown.
Numbers: noise 4 significant digits; percentages 2 decimals; timings 6 significant
digits. Return values retain full precision. Legacy report rows upgrade on write.

Human-readable short headings:
File RMS_pre RMS_post MAD_pre MAD_post RMS_% Npix Units Method MAD_% Exp_s
Readout Ngroups Nints Nframes Gap Group_s Frame_s Det Subarray

Latest format is aligned padded columns separated by TWO SPACES (not tabs).
Reader accepts old TSV and new compact format; maps short headings to original
internal keys. Tests were updated accordingly. User wanted narrower spacing.

Upgraded /Users/jrrigby1/AI_experiments/Temp/flicker_noise_rms.txt: 32 rows,
readout metadata recovered from matching RATE headers for all 32. Fit_method
unknown for all: files do not explicitly record it; did not guess from current
pipeline defaults. Original backup:
/Users/jrrigby1/AI_experiments/Temp/flicker_noise_rms.txt.before_upgrade

## Verification and remaining work

Passed at their respective migrations in JWST 3.0.0:
- 11 initial mask tests via package and redirects.
- 38 package tests after background migration; 40 legacy background/patch/workflow.
- 42 package tests after diagnostics; six synthetic masked-glow checks.
- two offline finder tests (PNG, DS9 regions, provenance; no network).
- artifact v0.3 scoring test verifies no reference load or baseline plotting.
- all 44 archive offline tests through package AND old entry points after
  transfer/retry consolidation.
- latest noise report: five diagnostic tests and 18 LEGGOS workflow tests passed;
  latest compact spacing revision: five diagnostic tests passed.

Commands:
python -m unittest discover -s tests -q          # from nirspec_ninja root
python -m nirspec_ninja.archive.tests.run_offline
python tests/check_masked_glow.py
python -m unittest test_fs_background test_nirspec_flicker_patch test_run_pipeline_LEGGOS_v3_0_0patched -q
# last command from Coding/JWST_pipeline

Full representative real-data reductions, live archive downloads, online WISE
retrieval still pending. The revised workflow first pass is recorded below.
Git initialization outstanding: /usr/bin/git fails because Apple developer tools
are unavailable; no Git at usual Homebrew paths. No commits made. Files saved on
disk. Do not claim entire current suite was rerun after latest spacing change.

Incidental environment warnings: matplotlib/fontconfig sandbox/cache warnings,
pysiaf cannot check online PRD version, some deprecations. Tests passed despite
those warnings. No Scratch code was used.

## First workflow migration — 2026-11-20

User approved the revised plan and requested its first pass. Migrated existing
jrmulti workers to pipeline/notebook_workers.py, retaining public names,
signatures, settings, and failure handling. Replaced local mask/extraction/output
helper imports with package imports. jrmulti.py is now a compatibility redirect.
The notebook changes only import text in cell 7; other cells, cell metadata,
outputs, and notebook serialization/formatting are unchanged.

Shared LEGGOS functions preflight_crds and run_pipeline_tasks moved to
pipeline/multiprocessing; apply_pipeline_parameters, check_pipeline_step_statuses,
and close_pipeline_result moved to pipeline/execution.py. Its own stage runners
and explicit parameter dictionaries remain in the LEGGOS script. Original
predownload_crds copied from supplied jrr into pipeline/multiprocessing; the
notebook uses it without changing its existing explicit-context behavior.
Detector1 RATEINTS helper moved to pipeline/detector1_outputs.py; old location
redirects. LEGGOS cleaning remains Detector1; notebook cleaning remains Spec2.

Hidden settings inventory is in src/nirspec_ninja/pipeline/README.md. None of
those settings were exposed or changed. Optional arguments for mask/Spec3
parameter visibility still require approval before the later edit.

Verification in JWST 3.0.0: 51 package tests passed (46 existing plus five
migration checks); 20 LEGGOS and cleaned-RATE checks passed with package root on
PYTHONPATH for the old jrmulti import. Spawned-process worker execution passed.
In-process ipykernel import/reload passed. Normal Jupyter kernel launch could
not be checked because sandbox socket binding was denied; test with a normal
notebook kernel remains pending, as do representative real-data reductions.

Original versions of the notebook, jrmulti.py, LEGGOS script, and
JWST_pipeline/detector1_outputs.py saved before changes in:
/Users/jrrigby1/AI_experiments/Temp/nirspec_ninja_migration_20261120
No Scratch code was read or used. No public names/signatures changed.

## Spec2 missing OPMODE fix — 2026-11-20

Reported crash jw04528013001_04101_00007_nrs1_rate_spec2_pipecrash.txt
occurred before pipeline execution: RATE meta.instrument.opmode was absent.
Notebook worker now reads this optional field with getattr(..., None). The
in-memory FS view still sets FIXEDSLIT; the cleaned diagnostic deletes the
temporary opmode when the source had no value, rather than recording FS mode.
Source RATE and public signatures are unchanged. Added an offline regression
that reproduces the missing attribute, runs mocked Spec2 cleaning, and verifies
unchanged source bytes, original diagnostic EXP_TYPE, absent diagnostic OPMODE,
and retained in-memory FS mode. Six pipeline migration tests and two legacy
cleaned-RATE diagnostic tests passed in the reference JWST 3.0.0 environment.
The user's miniforge environment and real exposure were not rerun here.

## Detector1/Spec2 settings visibility — 2026-11-20

User reported the migrated notebook ran correctly on two example directories,
then approved exposing Detector1/Spec2 settings and optional save_results
arguments. Added save_results=True after existing arguments to both notebook
workers. Detector1 exposes ramp_fit.skip=False and DET1_SAVE_RESULTS=True in its
stage cell; Spec2 exposes cleaning save_background=False and
SPEC2_SAVE_RESULTS=True. Existing SAVE_RATEINTS and scientific values unchanged.
Only current notebook source cells 19 and 21 changed; other cells, cell metadata,
and saved outputs were preserved, including user's edits since migration.

Explicit cleaning/save choices now take precedence over fallback values.
Detector1 reports a failure for cleaning enabled, ramp fitting skipped, or
save_results=False (a saved uncleaned RATE is required). Spec2 reports a failure
for cleaning skipped, rather than silently changing skip to False. Per-file mask
paths remain derived by the worker; cleaned-RATE diagnostics remain saved
independently of pipeline save_results. LEGGOS was not changed.

Seven pipeline migration tests passed, covering fallback defaults, explicit
Spec2 save flags, pipeline-level saving, stage requirement failures, missing
OPMODE, and spawned workers. All notebook code cells compiled. Executing old
and new stage setup with empty input selections confirmed identical parameter
choices after accounting for the newly explicit former defaults. Original
notebook saved before this edit under:
/Users/jrrigby1/AI_experiments/Temp/nirspec_ninja_settings_20261120
Real-data validation of this settings edit remains pending; the user's two
successful directories preceded it. Mask/Spec3/extraction visibility edits
remain pending.

Two legacy cleaned-RATE diagnostic tests also passed after this edit.

## Detector1 save-switch simplification — 2026-11-20

User requested removing the redundant DET1_SAVE_RESULTS notebook switch.
Removed only that assignment and its Pool call argument; preserved user's
comment edits and other notebook content. Restored the Detector1 worker's
original signature (save_rateints=False as its last argument). Its build_config
always uses save_results=True, with a comment explaining that masks and Spec2
require saved uncleaned RATE files. Spec2 saving options remain unchanged.

## Further notebook simplification — 2026-11-20

User approved removing explicit ramp_fit.skip=False entries and the
SPEC2_SAVE_RESULTS switch from the notebook. Both removed, together with the
Spec2 Pool argument. Preserved user's comment edits and other cell contents,
metadata, and outputs; notebook code compiles. Workers still enable ramp fitting
and save normal pipeline products by default, with explanatory comments/docs.
Spec2 save_results remains an optional API argument for diagnostic callers.
Retained cleaning skip settings (stage distinction) and save_background=False
(useful diagnostic option). No scientific defaults changed.

## Remaining notebook settings batch — 2026-11-20

User reported the notebook ran correctly after the earlier cleanup, then
approved mask, Spec3, and extraction visibility edits as one batch, explicitly
requesting no added comments. Implemented approved optional arguments:
make_rate_custom_mask(rate_file, input_dir, paramdict=None),
run_jwst_spec3(asn_file, output_dir, paramdict=None), and
run_extract(L3_destripedir, L3_extractdir, clip_nan=5).
Existing calls retain their defaults and return/error behavior.

Notebook mask_param_nested_dict exposes the seven existing scientific mask
choices; Spec3 shows outlier_detection.skip=False; extraction shows clip_nan=5.
Only source cells 20, 22, and 23 changed. No comments added; existing comments,
formatting, all other cells, metadata, and saved outputs preserved. Worker
mask paths/overwrite, Spec3 product saving/model cleanup, and extraction failure
handling stay unchanged. LEGGOS was not edited.

Verification: all notebook code compiles; 57 package tests and 20 legacy
LEGGOS/cleaned-RATE tests passed in reference JWST 3.0.0. Four new tests cover
omitted/empty parameters, explicit notebook defaults, overrides, input and
caller-dictionary preservation, output/model cleanup, extraction failures, and
actual notebook Pool task arguments. Tests use mocked scientific processing;
a representative notebook run for this latest batch remains pending.
Pre-edit notebook and workers saved in:
/Users/jrrigby1/AI_experiments/Temp/nirspec_ninja_remaining_settings_20261120

## Shared association writer and LEGGOS imports — 2026-11-20

User confirmed the settings-batch notebook ran with a restarted kernel and
approved this next batch. Migrated identical LEGGOS/supplied-jrr writel3asn to
pipeline/associations.py, preserving signature, comments/authorship, ordering,
serialization, return value, and output replacement. Both workflows now import
that function; supplied jrr.jrjwst retains a compatibility import. Notebook edit
was only `from jrr.jrjwst import writel3asn` ->
`from nirspec_ninja.pipeline.associations import writel3asn`, retaining its
existing comment. User subsequently reported possibly overwriting the notebook
while running a job; the old import remains supported, and no further notebook
writes were made in response.

LEGGOS now imports masks, runtime patches, noise diagnostics, and FS gradients
directly from package modules. Removed obsolete custom_mask_code_dir path
insertion, unused jrr_code_parent, and imports needed only by the moved writer.
Retained its CRDS environment preservation around mask import and its separate
stage runners, diagnostic hooks, and explicit parameter dictionaries.

Verification: 59 package tests and 20 LEGGOS/cleaned-RATE tests passed in JWST
3.0.0. Association output compared byte-for-byte with the original writer for
None, empty, and populated background lists. Legacy supplied-jrr import verified
identical to the shared function. Notebook compiles; cell metadata/outputs
preserved. Parameter dictionaries and untouched LEGGOS stage helpers compared
unchanged. Representative LEGGOS real-data validation still pending.
Pre-edit notebook, supplied jrjwst, LEGGOS script, and its test file saved in:
/Users/jrrigby1/AI_experiments/Temp/nirspec_ninja_associations_20261120

## LEGGOS real-data validation driver — 2026-11-20

User requested a script to run and return outputs. Created
Coding/JWST_pipeline/validate_LEGGOS.py and VALIDATE_LEGGOS.md. No notebook,
LEGGOS pipeline implementation, or existing reduction outputs edited.
Driver runs one supported target sequentially then with four workers (override
--cores), in a unique AI_experiments/Temp/leggos_validation_* directory.
Without a positional target it chooses the supported target with fewest UNCAL
files; --list-targets lists choices, --input-root overrides data root.
IFU/FS stages and FS background subtraction match run_all_IFS_targets, with
0.1 arcsec cubes, gradients/extraction on, trace screening and legacy background
subtraction off. Per-phase subprocesses inherit overrides for output and mask
geometry cache, including spawned worker imports. Existing CRDS cache retained.

Writes sequential.log, parallel.log, and validation_summary.json with metadata,
parameters, timings, inventories, and exact FITS image/association-products
comparisons. Optional --baseline-dir compares a previous target reduction.
Headers/ASDF metadata, binary tables, plots, and CSV values are not compared.
Four offline driver tests passed: output/cache redirection, fresh-process
inheritance, SCI/DQ/NaN and normalized-association comparison, and mocked two-phase
execution/report writing. CLI help checked. No real reduction was run here;
user will run the driver in the reduction environment and return reports.

## Documentation cleanup — 2026-11-20

User requested documentation cleanup while LEGGOS validation runs and an
explanation of the NRS2 grating issue. Updated 16 Markdown files across package
and JWST_pipeline docs: current summary above chronological migration history;
package import/module locations; archive test/install commands and accurate CI
status; noise-report formatting; current RATE diagnostic products; historical
reconstruction/migration labels; current mask-stage integration instructions;
and four mask-evidence links to their retained original diagnostics directory.
Replaced vague optical-fallback wording in current mask instructions with an
explanation of CRDS optical-model aperture predictions, retaining the code's
existing terminology only where identified explicitly. The separate NRS2
scientific issue remains unresolved and code behavior is unchanged.

Validated local links in all edited files. Hash comparisons confirmed no
notebook, LEGGOS driver, validation driver, or package Python source edits.
No regression suites rerun for this documentation-only change. Originals saved
under AI_experiments/Temp/nirspec_ninja_docs_20261120. User must be asked before
any subsequent notebook edits; current real-data LEGGOS outputs await review.

## Mask uses and tested grating/filter coverage — 2026-11-20

Checked mask consumers: NSClean/Clean Flicker Noise background selection,
before/after noise diagnostics, and optional artifact glow scoring. The custom
mask does not define extraction apertures or wavelength calibration. User agrees
predicted IFU illumination should be masked for cleaning. Documented desired
LEGGOS integration change without editing running code or notebook.
Added src/nirspec_ninja/masks/VALIDATION_COVERAGE.md with separate NRS1/NRS2 status for the
nine standard archive configurations. G235M/F170LP is the only recorded
real-data red-aperture comparison (one IRS2 exposure; NRS1 already-cleaned RATE,
NRS2 independent uncleaned RATE). Other grating/filter pairs are marked not
tested; PRISM extrapolation is currently disabled on both detectors, while the
need for additional illumination masking is unvalidated on each. Included evidence,
quantitative results, limits, and unresolved integration state. Evidence links
verified. No new scientific tests, Python changes, or notebook edits made.

## PRISM validation wording corrected — 2026-11-20

User challenged the earlier "not applicable" PRISM coverage label. Confirmed
that the extrapolation switch is disabled by disperser alone, for both NRS1
and NRS2. This is existing code behavior, not scientific evidence that PRISM
requires no extra masking. Updated coverage table and mask instructions to
state that need for additional illumination masking is unvalidated here for
each detector. No Python or notebook changes.

## Overnight restart prompt — 2026-11-20

User said the current LEGGOS test needs to run overnight and requested restart
instructions. RESTART_PROMPT.md captures context, constraints, current job
outputs, verification, scientific decisions, and pending work. Last read-only
check found two COSMIC.EYE validation directories, gf4le6x6 (newest) and n37qy9xj,
with no completed phases/comparisons recorded then. Do not infer failure or
completion from that snapshot; confirm current run on resume. No automation,
job restart, Python changes, or notebook edits were made.

## Notebook Spec2 weight cap — 2026-11-23

User instructed that the weight cap should always be on. run_jwst_spec2 now
explicitly installs and verifies the combined runtime patch before each cleaning
task, using DEFAULT_MAX_WEIGHT (4.0), and writes the patch status beside the mask.
The cap acts within the running worker, without editing installed JWST files.
Notebook contents were not changed. Restart its kernel before launching new
workers. Seven migration integration tests and five numerical patch tests passed
in the reference JWST 3.0.0 environment. Full real-data notebook rerun pending.
The proposed addition of the noise report to Spec2 has not been implemented.

## Quiet Spec2 success output — 2026-11-23

User requested output only for problems. Removed routine patch-verification,
cleaned-RATE-saved and Spec2-success prints from notebook_workers.py. Failure
reporting and tracebacks remain. Notebook unchanged; restart kernel/new workers
for the updated code. Processing and saved files are unchanged.

## Notebook Spec2 noise report — 2026-11-23

run_jwst_spec2 writes source.parent/flicker_noise_rms.txt (L2a) after successful
processing, comparing original RATE with retained *_rate_after_cleaning.fits
on the custom mask, passing the actual configured fit_method. Compact aligned
noise/MAD/readout format and process-safe upsert are shared with LEGGOS.
Successful reporting stays quiet; report failure uses worker traceback reporting.
Notebook unchanged. Restart kernel before launching new workers. Seven worker
integration tests (including real report creation with synthetic RATE data) and
five diagnostic tests passed in JWST 3.0.0. Real-data notebook rerun pending.
This supersedes the earlier note that Spec2 noise reporting was not implemented.

## Approved file cleanup — 2026-11-23

Executed CLEANUP_PLAN.md after user approval. Deleted 490 nominated old FITS
intermediates, 34 nominated caches, and seven partial FITS in two abandoned
runs; allocated size removed 37.28 GiB. Retained 126 representative intermediate
files and all protected scientific products, code, notebook, editor backup,
latest survey, and original reports. Archived 553 small evidence files and
recorded every deletion. Sources/retained files verified unchanged; helper
imports checked with -B. Audit records:
/Users/jrrigby1/AI_experiments/Temp/cleanup_evidence_20261123_221729
No further deletion authorization is inferred from this completed batch.

## MSA metadata migration and current decisions — 2026-11-23

Moved download_msa_metadata.py unchanged to archive; old location redirects.
Notebook not edited. Pending import is recorded in archive README. User accepts
LEGGOS overnight validation; do not schedule or propose more analysis of that run.
Compatibility cleanup means checking callers, changing imports with approval
where needed, then retiring redirects and unrelated staged jrr code. Keep the
redirects while external scripts/notebooks still rely on them.

## Remaining notebook jrr helper — 2026-11-23

User removed gethead and mark_CaII_Fraunhofer_lines imports. max_filesize moved
unchanged to pipeline.multiprocessing; supplied jrr.util redirects it. Notebook
not edited. New import: from nirspec_ninja.pipeline.multiprocessing import max_filesize.
Recursive largest-file selection and original empty-directory error checked.

## Compatibility cleanup — 2026-11-23

Updated 34 Python caller/test files in Coding to package imports; removed
import-only redirects custom_nirspec_mask_2.py, aperture_extrapolation.py,
nirspec_order_geometry.py (old mask project), JWST_pipeline/fs_background.py,
MAST_queries/download_metrics.py and nirspec_ninja/jrmulti.py. Command-line
wrappers are retained to preserve public commands. Staged jrr copy removed with
user authorization; no jrr backup retained. Edited-caller and retired-wrapper
record: Temp/migration_cleanup_20261123/manifest.json. Notebook untouched.
Historical old code not part of active workflows has not been reimplemented.
This current summary supersedes earlier pending jrr/import/column-score notes.

Cleanup verification: all 64 package tests and all 44 original archive offline
tests passed in JWST 3.0.0. Updated Python callers parsed successfully. Two
notebook-cell test fixtures now provide summarize_extract, matching the user's
current notebook; no notebook edits made. No live service calls or science
reductions performed. Caller audit is limited to Coding; unrelated scripts
outside it are not claimed checked. CLI wrappers retained deliberately.

## Live-service script — 2026-11-23

Added scripts/test_live_services.py for user-run WISE chart retrieval and one
public FITS download each from AWS/MAST, including size/FITS/hash checks and
separate logs/summary. Fresh output directories, one download worker, 512 MB
per-copy default limit; no pipeline rerun. Syntax/CLI help checked; actual live
services have not been invoked by the agent. scripts/README.md has usage.

## Live-service smoke tests passed — 2026-11-23

User ran scripts/test_live_services.py. Verified saved summary in
Temp/nirspec_live_services_h6cfhk_e/summary.json: WISE W1 retrieval/chart creation,
AWS download, MAST download and byte comparison all passed. Both RATE copies
are 83,969,280 bytes with identical SHA-256. Test exposure:
jw04125018001_02101_00001_nrs1_rate.fits. This completes the pending live-service
smoke checks; it is not a bulk query/download performance test. Broader mask
validation and repository setup remain. Finder chart is wise_finder.png in that
run directory; user visual review is not recorded.

## Finder partial image coverage — 2026-11-24

User requested support for partial NIRCam overlap. Local Cutout2D now uses
partial mode with NaN fill, retaining full chart and aperture polygons. Three
finder tests passed including partial coverage. Real copied NIRCam image checked
using an 8x coarser temporary image with corresponding WCS; output in
Temp/finder_partial_overlap_check/cosmiceye_finder_check.png. Original-resolution
render not run (about 11,600 pixels per side). WISE retrieval path unchanged.

## Finder FSzoom — 2026-11-24

Added approved --FSzoom flag: 100x60 arcsec rectangular view centered on IFU
and four adjacent fixed slits. --size overrides width and retains 5:3 aspect;
S200B1 outside this view remains in full chart and regions. Default full view
unchanged. Four finder tests passed. Original-resolution supplied NIRCam image
rendered successfully to Temp/cosmiceye_FSzoom.png and visually inspected.

## Finder documentation refreshed — 2026-11-24

Finder README consolidated: NIRSpec positional input vs NIRCam --image,
pointing headers, partial coverage, full/FSzoom/WISE example commands, all CLI
options, outputs, S200B1/region scope, warning explanations and memory limits.
Top-level README points to FSzoom instructions. Documentation-only update.

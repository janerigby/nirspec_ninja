# Pipeline processing on multiple CPUs

These tools run the JWST pipeline on many NIRSpec exposures. The fixed-slit notebook and LEGGOS
script use shared helpers while keeping their own reduction settings.

## Choose a workflow

- The [fixed-slit notebook](../../../Reduce-and-extract_fixedslit-multicores_v2.ipynb) reduces large amounts of NIRSpec fixed slit data. It   runs Detector1, mask generation, Spec2, Spec3, and extraction in separate cells.
  It applies 1/f cleaning to RATE data during Spec2.
- The [LEGGOS script](../../../../JWST_pipeline/run_pipeline_LEGGOS_v3.0.0patched.py) reduces NIRSpec IFS observations of galaxies.  It applies 1/f cleaning
  during Detector1 and has separate IFU and fixed-slit stage settings.

See [installation and notebook setup](../../../INSTALL.md) before running.
Set worker counts and scientific parameters in the workflow you use.

## Review fixed-slit background diagnostics

The LEGGOS fixed-slit branch creates gradient/continuum plots and a spatial
non-flatness CSV under each object's `FS/diagnostics/`. It also creates an
all-1 `FS_bkg_override.txt` starter table there and prints its path. Existing
override files are preserved, including user edits. Review the plots, then edit
the table if needed: `1` uses a slit with source masking, `0` excludes it from
the candidate background, and `-1` keeps the whole valid slit despite detections.

Rerun diagnostics using these local files, without repeating reduction:

```bash
python ~/AI_experiments/Coding/JWST_pipeline/rerun_LEGGOS_gradients.py --use-bkg-overrides
```

The normal rerun command retains the central crowding table. Add
`--targets OBJECT-DIRECTORY` for one object, or `--grating-smoothing 0.2` to
compare lighter grating-continuum smoothing (default 0.5 microns). PRISM is
unchanged. Starter values are for review, not automatic confirmation of clean sky.

The exact black curve is saved as `GRATING_FILTER_smoothed_bkg.csv`, including
inferred gap values and interpolation flags. In the LEGGOS driver, enable
`run_bgsub=True, background_source='FS'` to subtract this saved curve; use
`runfs_gradients=False` when reusing a reviewed CSV. Missing matching curves
cause an error. Subtraction defaults to off and retains the unsubtracted cube.

## Notebook stage settings

Parameter dictionaries beside the stage calls control the reduction. Explicit
settings take precedence over worker defaults.

| Stage | Settings and outputs |
| --- | --- |
| Detector1 | Requires ramp fitting and saves RATE. Cleaning is skipped here. `save_rateints` controls optional RATEINTS output. |
| Masks | Choose aperture masking, illumination refinement, and grating extrapolation in `mask_param_nested_dict`. Masks are overwritten when regenerated. |
| Spec2 | Requires cleaning and supplies each exposure's custom mask. Choose `fit_method` and cleaning-output flags in the step dictionary. Pipeline `save_results` defaults to True. |
| Spec3 | Configure the nested step dictionary, including outlier detection. Pipeline products are saved. |
| Extraction | `clip_nan` controls excluded edge rows; its default is 5. |

The notebook's extraction worker runs continuum source detection and writes
`trace_sig2` (corrected sigma) and `non-flat` (masked spatial non-flatness, %)
to each extracted CSV's comment header. These measurements identify spectra
that need inspection before use as sky; they do not alter extracted fluxes.
Unavailable scores are written as `nan` with a reason, and reported in the
worker summary. See [extraction diagnostics](../fixed_slit/README.md).

The mask worker defaults to WCS masking for MSA and fixed slits, illumination
refinement at 5 sigma, and grating extrapolation. Broad fixed-slit strips and
prism strips are off. These worker defaults differ from the standalone mask
builder; see [mask options and limits](../masks/README.md).

Spec2 installs the [combined 1/f patch](patches/README.md) with a weight cap of
4.0 in each worker. It saves `*_rate_after_cleaning.fits`, restores the original
observing-mode metadata in that product, and writes a patch record beside the
mask. `flicker_noise_rms.txt` in L2a compares the original and cleaned RATE
using the custom mask. See [noise report columns](../diagnostics/README.md).

Conflicting cleaning settings or skipped ramp fitting cause failures reported
through crash files and worker return values. Restart the notebook kernel
after updating the package so new workers use the updated code.

## Use the helpers in your own workflow

- [notebook_workers.py](notebook_workers.py) provides stage workers for Detector1,
  mask generation, Spec2, Spec3, and fixed-slit extraction.
- [multiprocessing](multiprocessing) provides `run_pipeline_tasks`,
  `preflight_crds`, and `predownload_crds` for process scheduling and reference
  preparation. Scientific settings remain in the calling workflow.
- [execution.py](execution.py) applies parameter dictionaries and checks step failures.
- [associations.py](associations.py) provides `writel3asn` for Level 3 science
  and background associations.
- [detector1_outputs.py](detector1_outputs.py) controls optional RATEINTS saving.

## Acknowledgments

The LEGGOS pipeline script was originally written by Brian Welch and later
updated by Jane Rigby. The multi-processor fixed-slit notebook was originally
written by Brian Welch and heavily edited by Jane Rigby.

The Level 3 association-making function, `writel3asn`, is adapted from 
[David Law's MIRI MRS notebook](https://github.com/STScI-MIRI/MRS-ExampleNB/blob/main/Flight_Notebook1/MRS_FlightNB1.ipynb).

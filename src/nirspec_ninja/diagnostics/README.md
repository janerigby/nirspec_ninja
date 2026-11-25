# Reduction diagnostics

These tools measure whether 1/f cleaning reduces detector noise and whether
spatial gradients in fixed-slit spectra could affect extraction or sky estimates.

## Compare noise before and after cleaning

Use the original RATE, the cleaned RATE, and the exposure's custom background mask:

```sh
python -m nirspec_ninja.diagnostics.flicker_noise_diagnostics \
  exposure_rate.fits exposure_rate_after_cleaning.fits exposure_mask.fits \
  --fit-method fft --report-file flicker_noise_rms.txt
```

Python callers can use `measure_noise` or `write_noise_report` from
`nirspec_ninja.diagnostics.flicker_noise_diagnostics`. Supply the actual cleaning
method (`median` or `fft`); an omitted method is recorded as unknown.

Both RATE images and the mask must have the same 2D shape, and the images
must use the same units. Measurements use the same finite background pixels
in both images, excluding DO_NOT_USE and NON_SCIENCE flags.

The notebook Spec2 worker writes `flicker_noise_rms.txt` in L2a. All exposures
in that directory share the report. LEGGOS also uses this report writer.

### Report columns

| Columns | Meaning |
| --- | --- |
| `File`, `Npix`, `Units` | Exposure filename, number of background pixels used, and SCI units. |
| `RMS_pre`, `RMS_post` | Noise RMS before and after cleaning. |
| `MAD_pre`, `MAD_post` | Scaled MAD noise estimates before and after cleaning. |
| `RMS_%`, `MAD_%` | Percent reduction; negative values indicate increased noise. |
| `Method` | Cleaning method. |
| `Exp_s`, `Readout` | Effective exposure time in seconds and readout pattern. |
| `Ngroups`, `Nints`, `Nframes`, `Gap` | Ramp sampling metadata. |
| `Group_s`, `Frame_s` | Group and frame times in seconds. |
| `Det`, `Subarray` | Detector and subarray. |

Missing metadata is recorded as unknown. Noise values have four significant
digits and percentages two decimal places; returned values retain full precision.

## Measure fixed-slit gradients

```python
from nirspec_ninja.diagnostics.fs_gradients import make_fs_gradient_report
make_fs_gradient_report(fs_root, output_dir, target, crowding_notes=None)
```

`fs_root` contains `L3/`; normally write results to `FS/diagnostics/`. The
L3-only plots show where a slit is contaminated and whether its remaining
rows can provide sky. Row coordinates are linear, with gaps left blank.
Slit order follows the finder chart: S200A2, S400A1, S1600A1, S200A1.

Review `FS_bkg_override.txt`, created with all entries set to `1`. Existing
files are never overwritten. Evaluations are per object directory and grating;
position-specific entries take precedence over the base target.

| Value | Use | Why |
| --- | --- | --- |
| `1` | Include the slit, excluding fitted source rows. | Measure sky away from a trace. |
| `0` | Plot as DO NOT USE; exclude from background estimates. | Keep contaminated or unsuitable slits out of the sky estimate. |
| `-1` | Include every otherwise valid row, even if a source is fitted. | Override a false detection caused by spatial background curvature. |

Trace shading encloses 99% of fitted Gaussian source flux. A fit can mistake
an edge gradient for a source, so inspect the S2D data before accepting its
mask. Ordinary pixel, edge, and gap masks always apply.

The spatial non-flatness score is `100 × max(|R − M|)/|M|`, where `M` is the
median spatial profile and `R` its 3-row running median. Source-masked rows
are excluded, except for `-1`. The score measures remaining spatial structure,
not source significance; it requires at least two contiguous 3-row windows,
and is undefined with too few usable rows or a
near-zero background. Scores and human choices are saved in `spatial_flatness.csv`.

In an initial LEGGOS comparison of 48 slit/grating evaluations from six objects,
the median score was 10.7% for human-rejected slits and 1.9% for accepted slits.
Eight rejected and 36 accepted evaluations had defined scores; four were
undefined. Rejected slits scored higher in 86% of rejected–accepted pairs,
but the distributions overlapped. The supplied CSVs had no whole-slit (`-1`)
overrides; one separately recomputed override scored 3.4% without source masking.

**Provisional criterion, requiring further testing:** consider rejecting a slit
as a background measurement when its masked spatial non-flatness exceeds
approximately **8%**. A strong residual gradient makes a spatially uniform
sky estimate less representative. This threshold has not been validated and
is not applied automatically. Test it across more objects, gratings, and slit
geometries, inspecting both false rejections and missed contamination. Lower
scores do not establish clean sky, and undefined scores require inspection.
Whole-slit overrides use unmasked source rows and should be assessed separately.

At each wavelength, the sky estimate is the spatial median of each good slit,
then an equal-weight median across slits. Equal weights prevent longer slits
from dominating. Grating spectra use a robust cubic smoothing spline with a
default scale of 0.5 μm: the background is expected to vary smoothly, and a
joint fit across both chips avoids independently fitted chip-edge jumps.
PRISM keeps its existing display. Inferred gap portions of the black curve
are dashed; JBT's predicted total background is overlaid when available.

`GRATING_FILTER_smoothed_bkg.csv` saves that exact black curve in MJy/sr,
including gap estimates and interpolation flags. The raw median is retained
separately. Review the choices and rerun before using the curve for subtraction.
In LEGGOS, use `rerun_LEGGOS_gradients.py --use-bkg-overrides` to apply edited
local files, or omit that flag to retain the central crowding table. Direct
calls use `use_bkg_override=True`. `grating_smoothing_um` changes the smoothing
scale; the selected file and scale are recorded in `diagnostic_settings.json`.

See [FS subtraction](../subtract_background/README.md) for using the reviewed curve.

## Exposure-specific investigations

`inspect_fs_gradients.py` and the scripts in `dataset_specific/` contain
investigations of particular exposures, including column structure and ramp
stripes. Their input and output paths are specific to those data sets; review
them before running on your own data.

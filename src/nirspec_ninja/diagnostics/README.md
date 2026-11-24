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

`fs_gradients.py` measures spatial structure in fixed-slit backgrounds.
Use crowding notes to identify slits affected by sources.

## Exposure-specific investigations

`inspect_fs_gradients.py` and the scripts in `dataset_specific/` contain
investigations of particular exposures, including column structure and ramp
stripes. Their input and output paths are specific to those data sets; review
them before running on your own data.

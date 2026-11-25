# Sky background estimation and subtraction

These tools estimate sky brightness from fixed-slit spectra or subtract a
supplied sky spectrum from NIRSpec IFU data.

## Estimate sky from the fixed slits

### Use the reviewed smooth background

First [inspect the fixed-slit diagnostics](../diagnostics/README.md), edit the
slit choices if needed, and rerun them. This matters because always-open slits
can contain sources, and smooth spatial gradients can trigger false detections.
As a provisional criterion, consider rejecting backgrounds with masked spatial
non-flatness above approximately 8%; this requires further testing and is not
an automatic cut. See the diagnostic guide for the initial comparison and limitations.
The black spectrum is saved as `FS/diagnostics/GRATING_FILTER_smoothed_bkg.csv`
in MJy/sr, including explicitly marked estimates through chip gaps.

In the LEGGOS driver, `run_bgsub=True, background_source='FS'` selects this
saved curve; `runfs_background=True` is equivalent. Set `runfs_gradients=False`
to reuse a reviewed CSV without regenerating it. Subtraction is off by default.
The matching object/grating/filter curve is evaluated on the IFU wavelength
grid and applied to each exposure before cube combination. The unsubtracted
cube is retained. Missing CSVs cause an error; wavelengths outside its coverage
are flagged or the cube is withheld. Source paths and hashes are recorded so
the subtraction can be reproduced. Fitted-background uncertainty is not yet
estimated and is not added to science ERR.

Python callers can supply `diagnostic_spectra_dir=...` to
`build_subtracted_cube`. This selects the saved curve instead of recomputing
sky from per-exposure slits.

### Measure separate per-exposure backgrounds

`fs_background.py` estimates contemporaneous sky surface brightness from
calibrated fixed-slit S2D pixels. It can exclude slits with detected sources,
checks exposure and wavelength-configuration matches, and tracks wavelength
coverage. It subtracts the estimated spectrum from calibrated IFU data before
combining exposures into a cube.

```python
from nirspec_ninja.subtract_background.fs_background import estimate_background, build_subtracted_cube
```

Use this method when the fixed slits provide a suitable sky measurement.
Background statistical errors are approximate and are not added to science
ERR because the estimate is shared by many science pixels.

See [fixed-slit source detection](../fixed_slit/README.md) for source screening.

## Subtract a supplied sky spectrum

`exp_bg_sub` interpolates an expected-background spectrum onto an existing
IFU cube's wavelength grid and subtracts it.

```python
from nirspec_ninja.subtract_background.subtraction import exp_bg_sub
exp_bg_sub("science_s3d.fits", "target_background.txt", "subtracted_s3d.fits")
```

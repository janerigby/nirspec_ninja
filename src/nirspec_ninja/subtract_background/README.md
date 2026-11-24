# Sky background estimation and subtraction

These tools estimate sky brightness from fixed-slit spectra or subtract a
supplied sky spectrum from NIRSpec IFU data.

## Estimate sky from the fixed slits

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

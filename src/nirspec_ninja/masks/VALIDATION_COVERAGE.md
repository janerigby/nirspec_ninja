# Predicted red IFU apertures: real-data validation coverage

Updated 2026-11-20. This table answers a specific question: have the predicted
IFU apertures beyond the nominal calibrated footprint been compared with
actual detector illumination? It does not describe general pipeline testing,
mask generation without a data comparison, or diffraction-order identification.

## What the mask is used for

The custom mask is supplied to Clean Flicker Noise/NSClean as `user_mask`:
1 means usable background and 0 means excluded. In the current workflows it
also supplies background pixels for the before/after noise report and can be
used by the separate artifact tools for glow scoring. It does not set spectral
extraction apertures or wavelength calibration.

The masking objective is to exclude all IFU illumination from the cleaning
background, including predicted apertures outside nominal calibrated coverage.
No nominal calibrated IFU footprint does not establish that no light falls on
the detector. CRDS optical-model predictions may therefore be used to construct
these exclusions; they do not make that light a calibrated extracted spectrum.

## Recorded comparisons

| Grating/filter | NRS1 red-aperture comparison | NRS2 predicted-aperture comparison |
| --- | --- | --- |
| G140M/F070LP | Not tested | Not tested |
| G140M/F100LP | Not tested | Not tested |
| G140H/F070LP | Not tested | Not tested |
| G140H/F100LP | Not tested | Not tested |
| G235M/F170LP | Tested on one exposure: mask compared with observed red trace bands in an already-cleaned RATE | Tested on the same exposure: predictions compared with an independently calibrated RATE made with flicker cleaning skipped |
| G235H/F170LP | Not tested | Not tested |
| G395M/F290LP | Not tested | Not tested |
| G395H/F290LP | Not tested | Not tested |
| PRISM/CLEAR | Extrapolation disabled; need for additional illumination masking not validated here | Extrapolation disabled; need for additional illumination masking not validated here |

“Not tested” means no real-data predicted-position comparison is recorded in
the evidence reviewed for this package. Successful notebook reductions, unit
tests, or generating a mask do not change this status. The table covers the
nine standard science configurations in the archive query list; it does not
claim coverage of every possible instrument configuration. For PRISM, the code
disables extrapolation on both detectors without checking which detector is
selected. That is current behavior, not evidence that no extra illuminated
regions need masking. The need for additional PRISM masking must be assessed
separately on NRS1 and NRS2 using actual data.

## G235M/F170LP evidence and limits

The tested IFU exposure is `jw04125018001_02101_00001`, full-frame IRS2, with
both detectors represented. On NRS1, the red-aperture mask excludes the observed
bands; 42.7% of pixels remain available as background. Its supplied RATE had
already been flicker-cleaned, so it is useful for comparing positions but not
for independently testing the cleaning correction.

On NRS2, an independent RATE was made with flicker cleaning skipped. The
predicted positions were not adjusted to match the data. All 30 tested slices
and all 30 tested 64-column blocks had positive central-aperture minus adjacent-
gap signal; the median excess was 0.00416 DN/s. The stacked profile peaked
3 pixels from the predicted center, and the signal occupied approximately the
predicted aperture width. The new mask retained 45.9% background pixels.
These comparisons support the predicted aperture locations for this exposure;
they are not a formal significance calculation accounting for noise covariance.

The predictions for different orders overlap. These checks do not separately
establish second- or third-order identities, filter transmission, or diffraction
efficiency. Old-versus-new NSClean corrections and spectral flux bias have not
been compared, and other grating/filter combinations and MOS data remain
unvalidated for these extensions.

Evidence is retained in the original diagnostics directory:

- [G235M/F170LP report](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/G235M_RESULTS.md)
- [Both-detector aperture projections](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/g235m_orders_both_detectors.png)
- [NRS2 data comparison](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/nrs2_order_validation.png)
- [NRS2 numerical results](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/nrs2_order_validation.json)
- [NRS1 red-mask comparison](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/g235m_mask_red_zoom.png)

## LEGGOS integration still to change

LEGGOS currently sets `mask_IFU=False` on NRS2 when nominal IFU coverage is
absent, including medium gratings and G140H/F070LP. That prevents the predicted
IFU apertures from being masked, even though `mask_grating_extrapolate=True`.
The desired change is to allow the grating predictions to reach the mask builder
while retaining PRISM handling. No code was changed during this documentation
update while the current validation job runs.

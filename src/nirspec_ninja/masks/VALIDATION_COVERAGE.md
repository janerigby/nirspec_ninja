# Predicted red IFU apertures: real-data validation coverage

Updated 2026-11-24. This table answers a specific question: have the predicted
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
| G140H/F100LP | One program 3843 IFU exposure checked: nominal red coverage reaches the edge; no admitted order-2/3 IFU pixels | One cleaning-skipped RATE checked: red illumination compared with optical predictions; all 232,786 first-order extension pixels and 43,890 order-2 pixels excluded by supplied mask |
| G235M/F170LP | Tested on one exposure: mask compared with observed red trace bands in an already-cleaned RATE | Tested on the same exposure: predictions compared with an independently calibrated RATE made with flicker cleaning skipped |
| G235H/F170LP | Two program 3843 IFU exposures checked: nominal red coverage reaches the edge; no admitted order-2/3 IFU pixels | Earlier cleaning-skipped and newer cleaned RATE checked: all 232,620 first-order extension pixels excluded by both supplied masks; no admitted order-2/3 IFU pixels |
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
been compared, and the remaining grating/filter combinations and MOS data remain
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

## Program 3843 high-grating checks, 2026-11-24

The interrupted G235H work was resumed with its now-available disperser
reference. Exposures jw03843002001_02101_00001 and
jw03843002001_04101_00001 were checked on both detectors. The earlier
NRS2 RATE has flicker cleaning skipped and gives median per-slice
aperture-minus-gap excess 0.00038758 DN/s beyond nominal endpoints,
positive in 24/30 slices. The newer RATE is already cleaned and supports
position checks only.

G140H/F100LP exposure jw03843004001_04101_00005 was checked on both
detectors after the user supplied missing optical references. Its NRS2
cleaning-skipped RATE gives median excess 0.00035863 DN/s, positive in
25/30 slices. Within admitted 0.97–5.3 micron wavelengths, projected
second-order IFU apertures intersect NRS2 and are fully excluded.
Overlapping paths prevent separate order identification from this signal.
G235H used admitted 1.66–5.3 micron wavelengths and has no projected
order-2/3 IFU pixels on either detector. Neither setup tests filter leakage
below the admitted wavelength bound or establishes diffraction efficiency.

Full optical IFU aperture pixels were checked against the actual supplied
masks, without added optical padding. These are exposure-specific results,
not general validation of all apertures, filter throughput, or correction
bias. Fixed-slit and extra-order failed-open coverage remain untested by
these full-aperture checks. No mask algorithm or pipeline wrapper changed.

[Report, plots, numerical results, and reproduction scripts](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/high_grating_validation/RESULTS.md).

## Fixed-slit red extent in the tested IFU exposures

AssignWcs fixed-slit footprints were reconstructed on both detectors for the
current program 3843 G235H/F170LP and G140H/F100LP cleaned RATE pairs.
Every on-detector nominal fixed-slit aperture reaches zero-based x=2047
(Image column 2048). Thus these exposures have no redward detector interval
beyond the fixed-slit WCS endpoints requiring continuation. This differs from
the IFU slices, which terminate early on NRS2. Blueward polynomial continuation is now disabled; these endpoint checks do not establish order identities
or independently validate all aperture widths or fixed-slit observing modes.
When mask_FS_by_wcs=True and grating extrapolation is active, the builder
already applies the same boundary-continuation algorithm to fixed slits.

[Fixed-slit coverage records and figures](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/high_grating_validation/fixed_slit_coverage.json).

## Blueward evidence review and policy correction, 2026-11-24

Polynomial aperture continuation now extends only toward larger detector x
(redward), for IFU, fixed-slit, and MOS footprints. Fit diagnostics record
direction=redward. The earlier blueward FS overlays illustrate superseded
two-sided extrapolation, not demonstrated extra flux.

Review of saved measurements found two preliminary NRS1 IFU blue-side hints:
G235H/F170LP, Image columns 9–128, median central-minus-gap 0.00009609 DN/s
with 18/30 candidate apertures positive; G140H/F100LP, columns 9–42,
0.00036442 DN/s with 18/29 positive. Both used uncleaned RATEs and approximate
center curves extrapolated from the supplied mask, not reconstructed optical
blue-side WCS predictions. No noise-covariance significance, repeated-exposure
confirmation, or physical admission below the wavelength-reference blue bound
was established. These are preliminary hints, not validated blueward masking.
No corresponding measured blueward FS excess was found in the reviewed
diagnostics. Independent confirmation would require additional exposures.

[Recorded G235H preliminary measurement](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/g235h_f170lp/edge_signal_supplied_geometry.json)
and [G140H preliminary measurement](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/high_grating_validation/g140h_f100lp/edge_signal_supplied_geometry.json).

# Custom NIRSpec masks for CleanFlickerNoise / NSClean

These tools make custom detector-level background masks for NIRSpec.  IFU, FS, and MOS modes are all supported. The purpose is to mask all illuminated pixels, so that clean flicker noise can fit 1/f noise from non-illuminated pixels. These tools are based on the mask-maker within clean_flicker_noise, but adds functionality:  
- mask emission a bit redward of the extraction range; that flux can bias the 1/f noise correction, and in particular, makes vertical artifacts at the blue and red edges of the fixed slits.
- mask the higher orders of the grating.  This is experimental
- better discriminate between illuminated and non-illuminated pixels in MSA fields. The pipelines is currently too aggressive, in our opinion, at masking MSA shutters; buffer pixels above and below our mask that need not be. This improves 1/f correction in crowded fields.

## Main interface

Use `make_custom_NIRSpec_mask` from **`nirspec_ninja.masks`**.

```python
from nirspec_ninja.masks import make_custom_NIRSpec_mask

background_mask = make_custom_NIRSpec_mask(
    '/path/to/exposure_uncal.fits',
    '/path/to/exposure_mask.fits',
    mask_MSA_by_wcs=True,
    mask_FS_by_wcs=True,
    mask_FS_strip=False,
    mask_FSprism=False,
    mask_grating_extrapolate=True,
    mask_by_illumination=True,
    input_dir='/path/to/uncal_and_msa_metadata',
    overwrite=True,
)
```

The return value is a Boolean array. The saved JWST FITS mask contains:

| Value | Meaning |
|---|---|
| `1` / `True` | Background pixel available to the noise fit |
| `0` / `False` | Excluded pixel |

Pass the saved filename to `CleanFlickerNoiseStep(user_mask=...)` or its
pipeline configuration. The mask builder does not itself run NSClean or
extract spectra. The module docstring documents the complete API.

## Inputs and aperture geometry

- Accepts UNCAL, RATE, RATEINTS, and compatible unextracted JWST datamodels.
- Aperture geometry comes from AssignWcs. Geometry-only masking uses metadata
  and relevant DQ flags, without reading science values.
- IFU/MSA selection normally follows `EXP_TYPE`. Fixed-slit WCS masking is
  enabled separately with `mask_FS_by_wcs=True`.
- Curved MOS masks require `mask_MSA_by_wcs=True`. MOS inputs also require the
  matching MSA metadata file and its exposure metadata; `input_dir` locates it.
- Handles IRS2 interleaved reference samples: output coordinates are the RATE
  science grid, not the expanded raw ramp grid.
- Existing exclusions remain excluded during illumination refinement.
  Failed-open shutter masking follows the existing pipeline/DQ machinery.
- `single_mask=True` writes a 2D mask. `False` writes geometric planes for the
  integrations, with supported per-integration DQ exclusions.

### Slitlets crossing MSA quadrants

For PRISM, mask generation automatically separates a slitlet's shutters by
quadrant when the selected metadata ID and dither contain a cross-quadrant
slitlet. This prevents JWST from interpreting quadrant-local columns such as
170, 171, and 1 as one enormous aperture. Each portion is projected in its own
quadrant, and their footprints are combined with the other mask components.
The portion without the primary shutter uses a background-only aperture with
local shutter coordinates; source offsets are not copied across quadrants.

This correction uses a private, temporary MSA metafile for mask WCS assignment
only. The original metafile, saved mask metadata, and extraction inputs remain
unchanged. Temporary files are private to each multiprocessing worker and are
removed afterward. Single-quadrant slitlets retain their existing definitions;
IFU and separately generated fixed-slit geometry use the existing path.

For gratings, crossing slitlets log a warning and retain the
existing geometry. Their masks may still exclude large regions: correction and
aperture extrapolation have not been validated for these cases. A warning is
not an assurance that such a mask is suitable for cleaning. PRISM normalization
supports primary-source and background-only crossing slitlets. Multiple primary
shutters or inconsistent shutter rows within a quadrant portion raise an error.

Real-data checks cover the three PRISM/NRS1 dithers in program 4233 observation
5: ordinary slitlet WCS definitions and separate fixed-slit masks are unchanged.
Other crossing configurations, including NRS2, still need real-data validation.

Broad fixed-slit strips and explicit boxes are additional exclusions. To keep
background gaps between individual fixed slits, use `mask_FS_by_wcs=True`,
`mask_FS_strip=False`, and `mask_FSprism=False`.

## Grating aperture extrapolation

`mask_grating_extrapolate=None` is the default. It automatically enables full
aperture continuation for IFU G235M/F170LP, G235H/F170LP, and G140H/F100LP,
which have recorded real-data extension checks. This includes redward coverage
on NRS2 and continues apertures redward to the high-x detector edge. Other configurations
and modes retain nominal masking unless explicitly enabled.

Use `mask_grating_extrapolate=False` to disable continuation, or `True` to
request it for other grating configurations/modes. Fixed-slit/MOS WCS
requirements below still apply when explicitly enabling it.
With PRISM selected, the flag is ignored and the normal prism geometry runs.
This applies to both NRS1 and NRS2. It does not establish that PRISM has no
additional illumination requiring masking; that question has not been
validated here for either detector.

For gratings, the algorithm fits each full aperture's upper and lower
boundaries and continues them redward only to the high-x detector edge. It masks the aperture,
including faint sky, rather than following only a bright object's spectrum.

Blueward polynomial continuation is disabled: it is not supported by the
validation performed here. Nominal WCS exclusions are retained. Local
illumination refinement may still exclude measured emission near nominal
apertures; this is distinct from extrapolating a band to the blue edge.
Optical fallback projections for detectors without nominal IFU coverage
remain wavelength-bounded optical predictions, not polynomial continuation.

The mask includes uncertainty margins while retaining dark gaps between
apertures where possible. Unreliable boundary fits produce a warning and a
wider mask. Fit diagnostics are saved beside the FITS mask as
`.extrapolation.json`.

When **no nominal IFU WCS falls on a detector**, such as G235M on NRS2, the
extension path instead predicts where the IFU apertures land using CRDS optical
reference files and the installed grating transforms, emitting a warning. It
projects full apertures for orders 1–3, from the wavelength-reference blue
bound to 5.3 microns.
This is not a measured filter-transmission or diffraction-efficiency model,
and these predictions are specific to IFU geometry.

Geometry caching is bypassed for extrapolation and for explicit reference
overrides. Extrapolation does not identify diffraction orders or provide a
new extraction wavelength calibration. Extra-order failed-open shutters are
not newly modeled.

## Illumination refinement

`mask_by_illumination=True` adds measured emission around the aperture masks:

- For RATE, it uses the image directly.
- For UNCAL, it builds a proxy from selected ramp group differences, with IRS2
  science-grid conversion and provisional median-offset removal. This does
  not run Detector1 calibration or ramp fitting.
- Alternatively, `mask_illumination_reference='/path/to/rate.fits'` supplies a
  matching 2D RATE. Do not combine this with `mask_by_illumination=True`.
- RATEINTS illumination refinement is unsupported.

Refinement searches near existing/predicted apertures. Defaults are a
5-sigma robust threshold, a 30-pixel local dispersion extension, and a
3-pixel cross-dispersion search margin. Without grating extrapolation, that
local search cannot recover emission hundreds of pixels beyond a footprint.
The geometric aperture mask remains excluded even where emission is faint.

A flicker-cleaned RATE can show trace locations. To assess changes to the
correction independently, start from the same UNCAL or a RATE generated with
flicker cleaning skipped.

## CRDS references

The builder respects existing `CRDS_PATH`, `CRDS_SERVER_URL`, and
`CRDS_CONTEXT` environment settings. If unset, it retains legacy defaults
pointing to `/Users/jrrigby1/Ref_files/crds_cache/jwst_ops` and
`jwst_1584.pmap`; configure the cache for your environment.

`wcs_reference_overrides` accepts explicit AssignWcs reference paths, plus
`msaoper` for failed-open masking. Use these when you need a specific set of reference files.

## Pipeline / notebook integration

The public builder automatically extends the validated IFU configurations.
For other configurations, add `mask_grating_extrapolate=True` to the mask
builder keywords, rather than Detector1 or Spec2 step parameters. The notebook
worker already explicitly enables continuation for gratings.

For LEGGOS, it belongs in the `options` dictionary inside
`make_uncal_custom_mask`. An existing "no nominal IFU coverage" branch can
force `mask_IFU=False` on NRS2, preventing optical-model IFU aperture predictions
for medium gratings and G140H/F070LP. Check the generated mask before relying
on this setting: this branch can leave predicted IFU illumination unmasked
on NRS2.

In the current notebook, this flag is visible in `mask_param_nested_dict` in
the mask-generation cell. The worker is
`nirspec_ninja.pipeline.notebook_workers.make_rate_custom_mask`. Restart the
kernel after worker-code changes. See the [fixed-slit notebook](../../../Reduce-and-extract_fixedslit-multicores_v2.ipynb).

## Assess masks for other grating/filter configurations

When assessing a new grating/filter combination, compare the mask with detector
illumination in a full-detector summary like the program 3843 high-grating figures. Include the full extent of
both NRS1 and NRS2, nominal WCS aperture boundaries, predicted extensions or
additional-order apertures, and the observed extra illumination. Include a
measured aperture-minus-adjacent-gap profile outside nominal coverage so that
geometric predictions and evidence for extra flux are clearly distinguished.
Do not substitute numerical results or cropped edge views for the full-detector
summary. Label cases with no outside-WCS region available to test.

Strong 1/f noise can obscure faint extensions. Use matching cleaned RATE files
for clearer visual comparisons when available, while retaining uncleaned RATE
measurements to assess changes caused by cleaning; these are versions of the
same exposure, not independent illumination measurements. Identify the exposure and
cleaning status, disclose display smoothing or background removal, and distinguish
slice-to-slice spread from statistical confidence. A cleaned image alone does
not establish that the cleaning preserves flux or that an order is detected.

Examples: [G235H/F170LP full-detector summary](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/high_grating_validation/g235h_f170lp_current_cleaned_summary.png)
and [G140H/F100LP full-detector summary](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/high_grating_validation/g140h_f100lp_cleaned_summary.png).

## Validation and current limits

See [the detector-by-detector validation table](VALIDATION_COVERAGE.md) for
which predicted red IFU apertures have been compared with actual illumination.
G235M/F170LP, G235H/F170LP, and G140H/F100LP have recorded real-data
position comparisons on a small number of IFU exposures. This mask is
for cleaning background selection and diagnostics, not extraction calibration.

Recorded aperture-position checks cover one G235M/F170LP, two G235H/F170LP,
and one G140H/F100LP IFU exposure on both detectors. Other grating/filter
combinations and MOS exposures remain unvalidated. Overlapping predicted paths do not establish individual
diffraction-order identities. An old-versus-new NSClean correction comparison
has not been performed.

See [the G235M report](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/G235M_RESULTS.md),
[order projections](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/g235m_orders_both_detectors.png),
[NRS2 validation](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/nrs2_order_validation.png), and
[the red-mask comparison](../../../../Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/g235m_mask_red_zoom.png).

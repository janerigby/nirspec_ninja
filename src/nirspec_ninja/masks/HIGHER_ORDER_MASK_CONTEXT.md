> This file retains the 2026-11-19 investigation and evidence.
> The integration-state section below was updated on 2026-11-20.

# Resume context: NIRSpec grating aperture masks

Updated 2026-11-19. The user stopped work for now and wants to resume later
with other grating/filter combinations, test coverage, and improve the method.
No further investigation or implementation is requested at this pause.

## User requirements

- Support IFU, fixed-slit, and MOS/MSA modes, on both NRS1 and NRS2.
- Mask full illuminated apertures, including sky, not only object spectra.
- Clean the whole detector: NSClean works better when its background fit is
  protected from illumination throughout the detector.
- Extend traces redward only to the high-x detector edge. Blueward
  polynomial continuation is disabled per the user’s corrected requirement. Avoid
  unnecessary vertical widening because the fast read direction is vertical
  and NSClean needs background pixels there.
- Batch operation must not require human inspection. Unreliable fits should
  produce a wider conservative mask and a recorded warning.
- Be terse; do not invent facts or assume uncertain instrument behavior.
  Ask about missing information. Use the Python/JWST installation and CRDS
  files in References_and_Tools. Ask about missing libraries instead of
  spending effort on workarounds.
- Representative MOS and other grating/filter data will be supplied later.

## Implemented code

Active builder: `custom_nirspec_mask_2.py`, function
`make_custom_NIRSpec_mask`. Older builder: `custom_nirspec_mask.py`.
Mask convention: True/1 = usable background; False/0 = excluded.

New opt-in `mask_grating_extrapolate=True`:

- Uses `aperture_extrapolation.py` to fit full aperture upper/lower edges with
  degree 1–3 polynomials. Degree selection uses withheld-end prediction errors.
- Continues the bands redward only to the high-x detector edge, with vertical margins derived
  from prediction errors. Poor fits warn and widen.
- Rejects anomalously wide/narrow columns only from fitting, retaining their
  original exclusions. Real IFU WCS footprints for slices 27 and 29 had
  endpoint islands that otherwise caused 30–40 pixel padding. This was fixed
  and covered by a regression check.
- Applies to IFU, curved MOS (`mask_MSA_by_wcs=True`), and fixed-slit WCS
  (`mask_FS_by_wcs=True`) footprints. It masks the aperture even if faint.
- With PRISM, the flag is ignored; normal prism geometry runs.
  This is current code behavior on both detectors. The need for additional
  PRISM illumination masking has not been validated here on NRS1 or NRS2.
- Default remains False. Geometry caching is bypassed when extrapolating.
- Writes `.extrapolation.json` beside the output mask.

For IFU detectors without nominal WCS coverage, the opt-in path warns and
uses `nirspec_order_geometry.py` to project orders 1–3. Its admitted physical
wavelength policy is the reference blue bound through 5.3 microns, with narrow
vertical padding and a local illumination search. This is provisional optical
geometry, not a filter-throughput model or independently calibrated higher-
order wavelength solution. Extra-order failed-open shutters are not modeled.

Existing illumination refinement remains available. UNCAL inputs use selected
raw group differences as a proxy; RATE inputs use the image. Full-aperture
exclusions are never reopened. Default local search alone cannot discover
missing traces hundreds of pixels beyond the nominal footprint.

Added `wcs_reference_overrides` for explicit reference paths. The builder now
respects preexisting CRDS environment settings. Legacy defaults remain if
these are unset. Reference overrides also bypass geometry caching.

## Runtime and references

Python:
`/Users/jrrigby1/AI_experiments/References_and_Tools/JWSTDP-3.0.0-10-py314-macos-arm64/bin/python`

JWST 3.0.0 runs successfully; NumPy/SciPy/Astropy/GWCS/ASDF/Matplotlib work.
Earlier context saying Python/CRDS execution was blocked is obsolete.

CRDS files:
`/Users/jrrigby1/AI_experiments/References_and_Tools/crds/`

The user added mappings, including `mappings/jwst/jwst_1584.pmap`.
Offline bestrefs selection works with `crds.core.rmap.get_best_references`
and uppercase keys from the datamodel's CRDS parameters. Many references
are in the directory root; the newly added NRS2 dark reference is in
`references/jwst_nirspec_dark_0490.fits`. The diagnostic RATE script handles
both locations. Optical selections are in `diagnostics/g235m_bestrefs.json`;
NRS2 Detector1 selections are in `diagnostics/nrs2_detector1_bestrefs.json`.
No installed JWST library was modified.

## Data and results

Example directory: `/Users/jrrigby1/AI_experiments/Temp/leggos1/`.
Exposure prefix: `jw04125018001_02101_00001`.
Configuration: G235M/F170LP, IFU, full-frame IRS2.
Both detector UNCAL files are present. NRS1 supplied RATE was already flicker-
cleaned; it remains useful for trace location, but not independent correction
validation.

New products in `diagnostics/`:

- `jw04125018001_02101_00001_nrs1_extended_user_mask.fits`
- `jw04125018001_02101_00001_nrs2_extended_user_mask.fits`
- Their adjacent `.extrapolation.json` fit/fallback records.
- `g235m_extrapolated_mask.png` and `g235m_mask_red_zoom.png`.
- `g235m_orders_both_detectors.png` / `.pdf`.
- `g235m_order_projection.json` and `g235m_projected_aperture_curves.npz`.
- `jw04125018001_02101_00001_nrs2_rate_no_flicker.fits`.
- `nrs2_order_validation.png` / `.pdf` and `.json`.
- `G235M_RESULTS.md` and `g235m_mask_validation.json`.

The new masks retain 42.7% background on NRS1 and 45.9% on NRS2. Minimum
background rows per column: 810 on NRS1, 918 on NRS2. The reported NRS1 band
near Image pixel (1850,445) is excluded. For allowed NRS1 pixels at x>=1665,
the RATE 99th percentile decreases from 0.05432 to 0.01290 DN/s. These are
mask-coverage measures, not evidence about correction bias.

NRS2 independent validation: ran Detector1 on its UNCAL with flicker cleaning
skipped, after the user supplied dark_0490. Dark subtraction, reference-pixel
correction, linearity correction, and ramp fitting record COMPLETE;
S_CLNFNS is SKIPPED. Measured central +/-8-pixel signal relative to adjacent
+/-23–28-pixel gaps, across 30 predicted slices and 30 64-column blocks.
Median excess is 0.00416 DN/s, positive in every tested slice and block. The
stacked signal occupies the predicted approximately +/-16.6-pixel aperture.
It weakens toward the right edge. These spatial checks are not a formal
noise-covariance significance calculation.

## What is established about diffraction orders

The nominal G235M/F170LP wavelength reference selects first order and
1.66–3.17 microns. JWST uses negative grating-order signs: -1, -2, -3.
Diagnostic projections change the grating equation's order, rather than
stretching detector x or bounding boxes. Lower optical transforms allow
NRS2 projections despite the nominal medium-resolution IFU WCS rejection.

At DS9 Image pixel (1850,445), slice 14 gives:
first order 3.46268 microns; second order 1.73134 microns;
third order 1.15423 microns. First and second are compatible with F170LP;
third is below its nominal passband here. This replaces the earlier rough
endpoint/dispersion estimate.

For admitted wavelengths, second-order paths reach NRS1's red portion and
NRS2; third-order paths reach NRS2. Actual NRS2 emission supports the aperture
geometry, but overlapping order paths and extended first order prevent
unique order identification. A separate third-order detection is NOT
established. Neither appearance nor geometric agreement alone establishes
higher-order identity or throughput.

Installed default CleanFlickerNoise geometry has the same nominal-WCS cutoff.
IFU masking uses finite wavelengths in nominal bounding boxes; FS/MOS masking
uses nominal boxes. Neither extends apertures to detector edges. Statistical
clipping can reject bright missed signal but cannot guarantee full aperture
coverage. IFU observations separately mask a full-width fixed-slit strip.

## Reproduction and checks

Run with the above Python:

1. `diagnostics/build_example_masks.py` creates public-builder masks from
   both UNCALs using offline-selected reference overrides and ramp illumination.
2. `diagnostics/make_uncleaned_rates.py nrs2` creates the independent RATE.
3. `diagnostics/project_orders.py` produces both-detector order projections;
   it uses the independent NRS2 RATE when present, otherwise a raw proxy.
4. `diagnostics/plot_mask_validation.py` plots the actual saved NRS1 mask.
5. `diagnostics/validate_nrs2_orders.py` produces NRS2 signal comparisons.

Eleven checks passed:
`python -m unittest test_aperture_extrapolation test_custom_nirspec_mask_2 -q`.
They cover continuation, full apertures, dark gaps, fit warnings, endpoint
islands, and existing illumination/cache behavior.

## Current integration state — 2026-11-20

Current tested/untested combinations and mask uses are documented in
[VALIDATION_COVERAGE.md](VALIDATION_COVERAGE.md). The user agrees that predicted
IFU illumination should be excluded from the cleaning background; the LEGGOS
branch change remains to be implemented after the running validation job.

The active notebook is
`Coding/nirspec_ninja/Reduce-and-extract_fixedslit-multicores_v2.ipynb`.
Its mask cell exposes `mask_grating_extrapolate=True` in
`mask_param_nested_dict`; workers are in
`nirspec_ninja.pipeline.notebook_workers`. Kernel-restarted example reductions
have succeeded. Ask before editing the notebook again.

LEGGOS also passes `mask_grating_extrapolate=True`, but its
`if not ifu_coverage:` branch sets `mask_IFU=False` on NRS2 for PRISM, medium
gratings, and G140H/F070LP. For gratings this prevents using CRDS optical-model
predictions of IFU apertures when nominal WCS is absent. These predictions
exclude possible illuminated regions from cleaning; they do not provide a new
extraction calibration or establish diffraction-order identities.

This behavior was preserved during consolidation. The desired change
would retain PRISM handling while allowing affected grating exposures to reach
the aperture-prediction code. That change has not yet been implemented.
G235M/F170LP evidence motivates a paired cleaning comparison; other affected
configurations still require representative validation.

## Plot requirement for future configurations

The user requires a full-detector summary plot whenever support or validation
expands to a new grating/filter combination. Show both NRS1 and NRS2, nominal
WCS apertures, predicted extensions/additional orders, and observed extra flux,
with an outside-WCS aperture-minus-gap profile. See the required-plots section
in [README.md](README.md) and its program 3843 examples. Matching cleaned RATEs
can improve visibility under strong 1/f noise; retain uncleaned measurements
to assess cleaning-induced changes and label cleaning status and display
processing. Cleaned and uncleaned products of one exposure are not independent
data; ask the user for additional exposures if independent confirmation is needed.

## Next work when the user resumes

1. Obtain representative exposures for other legal configurations and MOS
   examples with matching MSA metadata. Existing successful prism behavior
   must be retained. Data for other setups are currently unvetted.
2. Test fits, warnings, aperture widths, measured signal outside masks, and
   available background per column on both detectors. Improve fitting or
   illumination adjustment only where evidence requires it.
3. Compare old and new NSClean corrections from identical UNCAL inputs,
   measuring residual readout noise in independent dark gaps and checking
   bright/faint spectral flux for subtraction bias. This comparison has NOT
   been performed yet.
4. Investigate filter throughput and order discrimination if needed; current
   geometry does not establish second/third-order identities.
5. Review higher-order failed-open shutter coverage and the IFU optical
   fallback's admission policy. Do not present provisional geometry as fully
   validated across modes/configurations.

The nine standard science configurations in the user's MAST query list were
agreed: PRISM/CLEAR; G140M/H with F070LP or F100LP; G235M/H with F170LP;
G395M/H with F290LP. The list is in
`/Users/jrrigby1/AI_experiments/Coding/MAST_queries/MAST_query_large_nirspec_datasets.py`.

## Default integration, 2026-11-24

The public builder now defaults mask_grating_extrapolate to None: continuation
runs automatically when IFU masking is enabled for G235M/F170LP,
G235H/F170LP, or G140H/F100LP. False is an explicit off switch; True retains
the existing opt-in behavior for other configurations/modes. PRISM remains
unchanged. The notebook worker already passes True. This uses the existing
full-aperture boundary continuation, not an illumination-only detection or
a new extraction calibration. Real-data regeneration and public-builder
regression checks are retained in the high_grating_validation diagnostics.

Earlier descriptions of False as the default describe the prior API.

## Redward-only correction, 2026-11-24

The user rejected unsupported blueward extrapolation. extend_aperture now
adds exclusions only beyond the nominal high-x endpoint and records
direction=redward in fit diagnostics. This applies to IFU, fixed-slit, and
MOS callers. Nominal footprints, existing exclusions, and local measured
illumination refinement remain unchanged. Earlier two-sided diagnostic
figures and masks are historical products, not current continuation policy.

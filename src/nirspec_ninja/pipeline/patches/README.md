# JWST pipeline patches

These patches correct bugs in the JWST 3.0.0 Clean Flicker Noise step that
can introduce artifacts into NIRSpec spectra. Importing the package does not
apply them. Run the commands below in your JWST reduction Python environment.

## Saturated groups: protect valid spectra from contaminated noise fits

[jwst_3.0.0_saturated_groups.patch](jwst_3.0.0_saturated_groups.patch)
excludes SATURATED ramp groups from the noise fit, even with a supplied mask.
A pixel can have a valid faint slope in early groups but saturate later;
using those later groups can introduce vertical artifacts into usable spectra.
Valid early groups remain available.

This changes the installed JWST Python files. See [PATCH.md](PATCH.md),
[tcsh application instructions](APPLY_PATCH_tcsh.md), and
[JWST bug report #11062](https://github.com/spacetelescope/jwst/issues/11062).
Restart notebooks, Python sessions, and parallel reduction jobs after applying.

## Kernel centering: prevent a one-pixel shift in noise corrections

[patch_clean_flicker_noise.py](patch_clean_flicker_noise.py) corrects offsets
in both the weighting and smoothing kernels so they act at zero lag rather
than shifting their response by one pixel. See
[JWST bug report #11050](https://github.com/spacetelescope/jwst/issues/11050).

```sh
python -m nirspec_ninja.pipeline.patches.patch_clean_flicker_noise --check
python -m nirspec_ninja.pipeline.patches.patch_clean_flicker_noise --apply
python -m nirspec_ninja.pipeline.patches.patch_clean_flicker_noise --restore
```

`--check` (the default) checks without changing files. `--apply` backs up,
patches, and verifies the installed code. `--restore` restores the backup.
Restart notebooks, Python sessions, and parallel reduction jobs after applying or restoring.

## Combined patch: limit excessive influence from sparse background samples

[nirspec_flicker_patch.py](nirspec_flicker_patch.py) combines kernel centering
and saturated-group rejection with an inverse-density weight cap (default 4.0).
The cap prevents isolated usable background pixels from receiving excessive
weight in the noise fit. Its value remains experimental.

```python
from nirspec_ninja.pipeline.patches.nirspec_flicker_patch import install_patch
status = install_patch(max_weight=4.0)
```

Call `install_patch` in each cleaning worker. This patch applies only to the
running Python process; it does not edit installed JWST files. It verifies
the kernels, weight cap, and group-DQ exclusions. `save_patch_record` writes
settings and verification results beside the mask in a `*.flicker_patch.json` file.

LEGGOS also supplies calibrated pixel_dq to mask refinement to exclude
DO_NOT_USE, NONLINEAR, NO_LIN_CORR and NO_SAT_CHECK: unreliable pixels should
not drive the noise fit. That separate mask policy is not applied by
`install_patch` alone.

## Versions and validation

The Python patch tools accept only the tested JWST 3.0.0 code. Before upgrading,
check whether the official STScI JWST release includes these corrections. Stop
using each local correction once the official release includes it and its
regression tests pass.

The notebook Spec2 worker installs and verifies the combined patch before
cleaning, with the weight cap always enabled at 4.0.
Each exposure mask receives a companion *.flicker_patch.json record. Restart
notebook kernels after updating the worker code.

> Historical migration snapshot, 2026-11-19. See [current status](../../../MIGRATION_STATUS.md).
> Subsequent migrations moved notebook workers and added redirects in the supplied jrr.

# Mask migration — 2026-11-19

Authoritative implementation moved from Custom_masks_for_NIRSpec_CleanFlickerNoise:

- custom_nirspec_mask_2.py → src/nirspec_ninja/masks/custom_nirspec_mask_2.py
- aperture_extrapolation.py → src/nirspec_ninja/masks/aperture_extrapolation.py
- nirspec_order_geometry.py → src/nirspec_ninja/masks/nirspec_order_geometry.py

Only local imports and documentation paths changed in algorithm modules.
The original function name and signature are preserved.
The old modules alias the actual new module object, preserving legacy imports
and patches without a second implementation. The older custom_nirspec_mask.py
remains historical code and has not been migrated.

Validation in JWST 3.0.0: 11 regression tests passed with package imports;
the same 11 passed through legacy imports. Module identity and public alias
identity were checked. Tests use synthetic data/mocked WCS and no CRDS access;
full real-exposure scientific equivalence has not yet been tested.

No Scratch code was read or migrated; only its file modification timestamps
were inspected. Supplied jrmulti.py and notebook were newer than corresponding
versions found elsewhere. Supplied jrr is staged unchanged for a later phase.

Git initialization is outstanding: /usr/bin/git fails because Apple developer
tools are unavailable, and no Git executable exists in the usual Homebrew
locations. Once Git is available, run git init in this repository.

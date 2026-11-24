# Apply the saturation patch using tcsh

Activate the Python environment used for JWST. This installs the saturation-only fix described in [issue #11062](https://github.com/spacetelescope/jwst/issues/11062). Updated LEGGOS already applies its broader correction automatically within each parallel reduction job.

Locate the package and check that the patch applies:

```tcsh
set patch_file = "/Users/jrrigby1/AI_experiments/Coding/nirspec_ninja/src/nirspec_ninja/pipeline/patches/jwst_3.0.0_saturated_groups.patch"
set jwst_pkg_parent = "`python -c 'from pathlib import Path; import jwst; print(Path(jwst.__file__).resolve().parent.parent)'`"

patch --dry-run -p1 -d "$jwst_pkg_parent" < "$patch_file"
```

The --dry-run command checks whether the changes can be applied without editing
any files. Only if it succeeds, save a backup of the Python file and apply:

```tcsh
cp -n "$jwst_pkg_parent/jwst/clean_flicker_noise/clean_flicker_noise.py" "$jwst_pkg_parent/jwst/clean_flicker_noise/clean_flicker_noise.py.before_saturation_fix"
patch -p1 -d "$jwst_pkg_parent" < "$patch_file"
```

Verify by starting a new Python session with this command:

```tcsh
python -c 'import numpy as np; from jwst.clean_flicker_noise.clean_flicker_noise import _mask_unusable; mask = np.ones(4, dtype=bool); _mask_unusable(mask, np.array([0, 1, 2, 4], dtype=np.uint8)); assert mask.tolist() == [True, False, False, False]; print("PASS: saturated groups are excluded")'
```

Restart notebooks and parallel reduction jobs afterward. This change to the installed JWST files adds saturation rejection only; the configurable weight cap and additional calibrated pixel-DQ exclusions come from the LEGGOS integration.

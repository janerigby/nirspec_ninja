# CleanFlickerNoise saturation patch

This JWST 3.0.0 patch adds `SATURATED` to the group-DQ flags rejected by `clean_flicker_noise._mask_unusable`, alongside `DO_NOT_USE` and `JUMP_DET`. It applies even when a user supplies a background mask. Valid early groups remain usable if a pixel saturates later in the integration.

The patch is needed because a cosmic-ray event can saturate later groups while ramp fitting still recovers a faint, valid slope from earlier groups. Admitting the saturated groups to the noise fit can introduce vertical artifacts into otherwise usable spectra. See [JWST issue #11062](https://github.com/spacetelescope/jwst/issues/11062) for the diagnosis and example results.

The file describing the saturation-only code changes is `jwst_3.0.0_saturated_groups.patch`; it includes regression tests. Restart notebooks, Python sessions, and parallel reduction jobs after applying it to the installed package.

The updated LEGGOS workflow also applies a verified correction applied within each running Python session that retains both kernel-centering fixes, caps inverse-density weights, and excludes calibrated pixels flagged `NONLINEAR`, `NO_LIN_CORR`, or `NO_SAT_CHECK`. Its user-changeable `nsclean_max_weight` defaults to **4.0**. Each exposure's mask has a companion `*.flicker_patch.json` recording the cap and verification results. These additional changes are separate from the saturation-only patch.

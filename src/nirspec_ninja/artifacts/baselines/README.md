# Historical reference scores

`reference_scores_v02.csv` preserves scores from the earlier v0.2 artifact
method for MSA_bad, MSA_okay, IFU_bad, IFU_okay, and Mystery_set1 samples.
These reference-normalized column scores are incompatible with the current
noise-normalized scores. Current scoring and plotting do not load this file
by default. Rescore the exposures with the current method before comparing them
with new results.

`manifest.json` records checksums and group counts. The CSV is a copy of
`fresh_v2/many_masks_expanded_results_v02/scores.csv`; that original directory
contains the historical diagnostics. One exposure has no mask and no glow score.

Keep this file as a historical record and save new results in separate run
directories. See [current scoring instructions](../README.md).

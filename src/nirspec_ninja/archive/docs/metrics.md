# Comparing AWS and MAST

All Python download workflows write the same `attempts.csv` and `summary.csv`
under timestamped `transfer_metrics/` directories (beneath the output directory
for program/file-list downloads; beneath the working directory for
the science query). Existing science-query run metrics remain available.

The summary separates source and operation and reports attempts, successes,
failures, retries, skipped files, timeouts, timeout fraction, bytes, elapsed
wall time, and MB/s. Per-attempt rows include filename, duration, error, worker
count, and configured timeout when known. HEAD checks are separate from transfers.
Files reused locally do not contribute bytes or attempts to transfer rates.
Results are saved after each attempt, including failures before a run aborts.

Use `MAST --direct-mast` with the science-query tool to measure actual MAST
downloads; default MAST curl generation provides no transfer measurements.
The file-list tool accepts `--source AWS`, `--source MAST`, or `--source AUTO`
(AWS with MAST fallback). The program downloader already accepts either source.

Compare identical filenames in fresh, separate destinations with equal worker
counts and timeout policies. The science-query tool reuses complete files, so
running both sources against one output tree will skip the second transfer.
Use the file-list downloader for paired downloads; no separate benchmark is
needed. Run from the repository root:

```sh
python MAST_download_rawfilelist.py filenames.txt --source AWS --workers 3 -o downloads/aws_test
python MAST_download_rawfilelist.py filenames.txt --source MAST --workers 3 -o downloads/mast_test
```

Compare each destination's `transfer_metrics/*/summary.csv` download rows.
Use fresh destinations for each repetition and alternate which source runs
first. Repeat with enough files to estimate timeout frequency. File-list AUTO
fallback selects a harder subset for MAST and is not
an unbiased source comparison.

`timeout_fraction` is recognized timeout attempts divided by all explicit
attempts, not a prediction for future runs. Wrapped timeout exceptions are
recognized; Astroquery may hide internal retries or return failures in manifests,
so these counts cannot capture all internal timeouts. Blank timeout settings
mean library defaults, not unlimited time. Library-backed failed transfers have
unknown partial byte counts (recorded as zero); recovery transfers count bytes
actually written, including partial attempts and resumed data.

Wall throughput spans first attempt to last completion per source/operation,
including intervening retries and idle gaps. MB per attempt-second divides bytes
by summed attempt durations and excludes retry sleeps; it is not parallel wall
throughput. Archive metadata queries use MAST for both sources and are not AWS
transfer measurements. Metrics writing itself is outside measured attempt time.

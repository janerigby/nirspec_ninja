# Downloading known JWST filenames

`MAST_download_rawfilelist.py` replaces the two recovery scripts and the
missing-file shell downloader. It deduplicates the union of all supplied lists.
Use `--dry-run` to validate inputs and inspect destinations without transfers.

```sh
# Bare filenames or paths, downloaded into one flat directory
python MAST_download_rawfilelist.py filenames.txt -o downloads/recovered

# Restore date/visit paths and cross-check an old recursive inventory
python MAST_download_rawfilelist.py paths.txt --inventory inventory.txt \
  --layout date-visit --require-paths -o downloads/recovered_tree
```

Inputs may contain blank lines and `#` comments. Path lists use
`YYYY-MM-DD/VID###########/filename`, optionally preceded by an old root path.
Dates and visit IDs are validated; conflicting destinations are errors.
`--inventory` accepts recursive `ls -lR` or `ls -lhR` text and may be repeated.
Inventory paths are never trusted as destinations: files may have been moved
into administrative holding directories. Explicit path lists are authoritative.

In `date-visit` mode, filenames without explicit paths are placed using FITS
`DATE-OBS` after download. `--require-paths` instead rejects missing destinations
before downloading. `--visit-dates FILE` optionally supplies date/VID hints from
a recovery metadata file. Flat output ignores old directory paths.

Defaults: three workers, AWS then MAST, 30-second connection timeout,
120-second read timeout, five attempts per source with progressive delays,
and one final gap-fill pass. `--source AWS` or `--source MAST` disables fallback;
`--allow-missing-lists` explicitly permits missing input files.

Existing files are retained when their sizes match the archive. Transfers
resume `.part` files and replace final files atomically after size verification.
Adjacent `.download.json` checkpoints support final verification. A completed
file still awaiting date/visit placement remains a recovery gap.
`recovery_report.json` and `missing_files.txt` record the final result.
Size matching is a completeness check, not a checksum verification.

## Previous recovery sets

The former small recovery used 357 unique filenames; its local
`dead_hd_list.txt` input contained 303 of those names. It restored the
`RAW_prism_extra` date/visit tree. The large wrapper restored 3,238 filenames
into `RAW_prism_new`; explicit paths and the recursive inventory agreed.
These are historical inventory counts, not current recovery status.

Local lists now live in ignored `local_inputs/`. To retain either previous
recovery layout, supply its original inputs, `--layout date-visit`, and the
original destination explicitly. No dataset-specific or personal absolute
paths remain in the downloader.

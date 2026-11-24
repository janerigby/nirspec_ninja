# Archive queries and downloads

These tools find large NIRSpec data sets and download JWST files from AWS or
MAST. They support parallel queries and downloads, retries, and recovery from
filename lists.

## Choose a download tool

- `MAST_query_large_nirspec_datasets.py` selects large NIRSpec samples (for example, all the NIRSpec spectra taken with the prism as the disperser) and  downloads data or prepares MAST curl scripts.
- `MAST_download_PID.py` downloads raw data for specified JWST programs.
- `MAST_download_rawfilelist.py` recovers raw files from filename lists.  This is useful for large datasets, if you let the archive be your backup.

```sh
python -m nirspec_ninja.archive.MAST_download_PID 04125 AWS --instrument nirspec
python -m nirspec_ninja.archive.MAST_download_rawfilelist filenames.txt --source AUTO -o downloads/recovered
python -m nirspec_ninja.archive.MAST_query_large_nirspec_datasets AWS --query-workers 2 --workers 3
```

See [tool options and examples](TOOLS.md) for source selection, output layouts,
and restarting downloads. Retry and restart behavior differs between tools;
program-ID downloads stop after a failed transfer.

## Download MSA metadata companions

MSA exposures need the companion FITS file named by `MSAMETFL` in the UNCAL
header. This tool downloads it beside the UNCAL file, preserves existing files,
and validates newly downloaded files. Non-MSA exposures are skipped.

```sh
python -m nirspec_ninja.archive.download_msa_metadata /path/to/*_uncal.fits
```

```python
from nirspec_ninja.archive.download_msa_metadata import download_metadata
```

## Check your setup

See [live-service checks](../../../scripts/README.md) to test WISE retrieval
and downloads from both archives. Offline tests use simulated archive responses:

```sh
python -m nirspec_ninja.archive.tests.run_offline
```

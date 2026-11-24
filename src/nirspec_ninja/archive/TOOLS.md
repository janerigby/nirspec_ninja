# JWST MAST query and download tools

| Tool | Purpose |
| --- | --- |
| `MAST_query_large_nirspec_datasets.py` | Select public JWST NIRSpec observations using filters like grating/filter and exposure time; download UNCAL files or generate MAST curl scripts. |
| `MAST_download_PID.py` | Download raw data for a JWST program, organized by PID and target; supports LEGGOS and other programs. |
| `MAST_download_rawfilelist.py` | Download raw JWST data from a list of filenames. |

## Setup and use

Use the existing JWST environment with nirspec_ninja installed; see
[installation](../../../INSTALL.md). From any working directory:

```sh
python -m pip install --no-deps --no-build-isolation -e /Users/jrrigby1/AI_experiments/Coding/nirspec_ninja

python -m nirspec_ninja.archive.MAST_query_large_nirspec_datasets AWS --query-workers 2 --workers 3
python -m nirspec_ninja.archive.MAST_download_PID 04125 AWS --instrument nirspec
python -m nirspec_ninja.archive.MAST_download_PID 04125 03843 AWS --instrument nirspec
python -m nirspec_ninja.archive.MAST_download_rawfilelist filenames.txt --source AUTO -o downloads/recovered
```

Configure the large-query tool's scientific settings and `RUN_MODE` near its
header. Current selection: F070LP/G140M, exposure >2,000 seconds, release year
>2022.5, absolute Galactic latitude >20 degrees, and either fixed or moving
targets. Prism selection retains NRS1. `DOWNLOAD_ONLY` reuses the final product
manifest. Query workers and transfer workers are independent (defaults 2 and 3).
MAST mode creates curl scripts; add `--direct-mast` to transfer files directly.
Existing checkpoint and output paths are preserved for restart compatibility.

The PID tool accepts `--instrument miri`, `--instrument all`, and
`--output-root PATH`. Supply one or more PIDs; each gets its own PID/target
directory and transfer metrics. The default is `LEGGOS_PIDS = [4125, 3843]`
(confirmed interactively when no PIDs are supplied), and the default instrument
is NIRSpec. Use `--source AWS` or `--source MAST`, or a trailing source as above.
Repeated PIDs are downloaded once. A program with no matching data is skipped;
a download error stops the run, with completed files available for restart.

See [file-list recovery](docs/rawfilelist.md), [AWS/MAST metrics](docs/metrics.md),
and [moving-target metadata limitations](docs/known_archive_issues.md).
To compare sources, download the same list into fresh, separate directories:

```sh
python -m nirspec_ninja.archive.MAST_download_rawfilelist filenames.txt --source AWS -o downloads/aws_test
python -m nirspec_ninja.archive.MAST_download_rawfilelist filenames.txt --source MAST -o downloads/mast_test
```

Compare the `transfer_metrics/*/summary.csv` files. Keep worker counts equal;
use enough files and repeat runs before drawing conclusions about timeouts.

## Repository contents

Code, tests, docs, and small example inventories belong in Git. Local lists
live in ignored `local_inputs/`; science data, checkpoints, run metrics, and
transfer logs also remain local. Output paths are relative to the launch
directory unless explicitly supplied. Do not discard restart manifests when
resuming an existing download.

The former query and LEGGO scripts have been renamed to the first two tools
above. The recovery wrappers and missing-file shell script are consolidated
into `MAST_download_rawfilelist.py`; update old commands to use explicit inputs
and output paths. Reconstruction prompts and one-off reports have been retired.

## Offline checks

```sh
python -m nirspec_ninja.archive.tests.run_offline
```

The 44-test suite blocks outgoing network connections and uses temporary files
and mocked archive responses. The old CI workflow and requirements file remain
in Coding/MAST_queries; package CI is not yet established. See
[testing](docs/testing.md) for coverage and maintenance guidance.

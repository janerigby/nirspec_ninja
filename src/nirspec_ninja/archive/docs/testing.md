# Regression testing

Run the default suite from the repository root:

```sh
python -m nirspec_ninja.archive.tests.run_offline
```

The runner blocks outgoing socket connections. Tests use temporary files,
synthetic Astropy tables, small FITS headers, and archive response doubles.
No science downloads or live archive accounts are needed. Plain unittest
also works from the package root:
`python -m unittest discover -s src/nirspec_ninja/archive/tests -t src -v`.

## Protection now implemented

The archive suite contains 44 tests and passed at migration in the local
Python 3.14 environment. This count is separate from the package suite.

| Component | Approximate line coverage | Protected behavior |
| --- | ---: | --- |
| Common metrics | 99% | Concurrent persistence, timeout recognition, exact wall/attempt throughput, retry gaps, resumed bytes, cached-file exclusion, and separate source/operation denominators. |
| Large query | 78% | Scientific thresholds, fixed/moving criteria, prism detector selection, UNCAL/MIRROR filtering, checkpoint reuse, actual transfer retries, complete success/failure workflows, saved manifests, and DOWNLOAD_ONLY restart. |
| Program downloader | 89% | Single/multiple PID selection, confirmed LEGGOS defaults, instrument selection, calibration masks, target joins and missing metadata, cached files, complete workflows, manifest failures, and persisted timeout metrics. |
| File-list downloader | 86% | Inventory/path validation, dry-run, HEAD retries and size checks, interrupted/resumed transfers, ignored Range, exhaustion, fallback, existing-file retention, real FITS date recovery, final gap pass, reports, and exit status. |

Figures were measured with Python's standard-library `trace`, including
worker threads. They are approximate line coverage, include module initialization,
and do not measure branches or CLI subprocess coverage. Maintain observable
behavior and useful failure scenarios rather than pursuing 100% coverage.

The failed-transfer regression found and fixed a bug in the large-query tool:
size verification now occurs in staging before an existing file is replaced.

## Continuous integration

The historical workflow remains in the original
`Coding/MAST_queries/.github/workflows/tests.yml`, alongside its
`requirements.txt`. It invokes the old repository test path. Neither file was
migrated into nirspec_ninja, and package CI has not been established. The earlier
YAML was checked locally; its GitHub matrix was not verified running. Use the
installed JWST environment and the package test command above for current checks.

## Maintenance guidance

Add a small failing regression when fixing a bug, then preserve that scenario.
Mock archive/network boundaries and assert files, selected filenames, reports,
metrics, and failure status. Avoid assertions about log wording, thread order,
or a specific private helper call unless it protects a public contract.
Shared response fixtures live in `tests/helpers.py`.

Remaining gaps include interruption during atomic writes, unusual corrupted
checkpoint schemas, and archive/library changes. Live AWS/MAST smoke checks
should remain optional and separate from CI. Offline tests do not establish
current server compatibility, real performance, or timeout probability.

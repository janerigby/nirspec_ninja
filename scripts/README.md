# Live-service smoke tests

Check that your reduction environment can retrieve WISE images and download
public JWST files from both AWS and MAST. Run this script from the package root
in your installed JWST reduction Python environment:

```sh
python scripts/test_live_services.py /path/to/public_nirspec_exposure.fits
```

Supply a local public NIRSpec exposure with pointing metadata. The script:

1. Queries IRSA for WISE W1, downloads an image, and creates a finder PNG,
   DS9 regions and provenance record using the package finder tool.
2. Downloads the corresponding public RATE file once from AWS and once from
   MAST, with one worker each, using the filename-list downloader.
3. Checks file sizes and FITS structure, then compares SHA-256 hashes of both
   archive copies. Each test records pass/fail separately.

Outputs and logs go in a new nirspec_live_services_* directory under
AI_experiments/Temp. summary.json contains results. No existing reduction is
rerun or overwritten. The script does not test complete program-query workflows
or scientific accuracy of the finder footprints.

--archive-file selects a different public JWST FITS basename. By default an
UNCAL input selects its corresponding RATE to reduce download size; another
input uses its own basename. --max-download-mb defaults to 512 MB per archive
copy and checks advertised size before downloading. --timeout defaults to 900
seconds per chart/download command. --output must name a new directory.

The WISE cache is fresh for each test run so a prior image does not bypass
retrieval. Logs are retained when a service fails; the other tests still run.
Open wise_finder.png to check the chart visually after the script finishes.

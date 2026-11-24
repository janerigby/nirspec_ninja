# Known archive issues

## Moving-target metadata

Some JWST observations have FITS `TARGTYPE` indicating a moving target while
MAST CAOM reports `mtFlag=False`. A local check on October 23, 2026 found this
for the examples in [the original report](mtFlag_bug_report.txt); the current
archive status has not been rechecked.

Using `mtFlag=True` can therefore omit moving-target observations, and
`mtFlag=False` can include them. The large-dataset tool currently defaults to
`TARGET_TYPE="EITHER"`. Verify FITS headers when accurate target classification
matters; the scripts do not correct archive classifications automatically.

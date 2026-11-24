#!/usr/bin/env python3
"""Download MSAMETFL companions for NIRSpec MSA UNCAL files (requires astropy).
Prompt:  A script that, for a set of input JWST nirspec UNCAL files, finds all those with EXP_TYPE indicating
they are MSA observations, get the auxilary filename from MSAMETFL header keyword, and downloads
 it to the same directory as the UNCAL file

Usage: python download_msa_metadata.py /path/to/*_uncal.fits
"""

import argparse
from pathlib import Path
import shutil
import sys
from urllib.parse import urlencode
from urllib.request import urlopen

from astropy.io import fits


def download_metadata(uncal):
    header = fits.getheader(uncal, ext=0)
    exp_type = str(header.get("EXP_TYPE", "")).strip().upper()
    if not exp_type.startswith("NRS_") or "MSA" not in exp_type:
        print(f"Skip {uncal.name}: EXP_TYPE={exp_type}")
        return

    filename = str(header.get("MSAMETFL", "")).strip()
    if filename.upper() in ("", "N/A", "NONE", "NULL"):
        raise ValueError(f"{uncal}: missing MSAMETFL")
    if Path(filename).name != filename:
        raise ValueError(f"{uncal}: MSAMETFL must be a filename: {filename!r}")

    destination = uncal.parent / filename
    if destination.exists():
        print(f"Already exists: {destination}")
        return

    url = "https://mast.stsci.edu/api/v0.1/Download/file?" + urlencode(
        {"uri": f"mast:JWST/product/{filename}"}
    )
    temporary = destination.with_name(destination.name + ".part")
    try:
        with urlopen(url, timeout=120) as response, temporary.open("wb") as output:
            shutil.copyfileobj(response, output)
        # Reject error pages or other responses that are not FITS files.
        fits.getheader(temporary, ext=0)
        temporary.replace(destination)
    finally:
        temporary.unlink(missing_ok=True)
    print(f"Downloaded: {destination}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("uncal", nargs="+", type=Path, help="Input UNCAL FITS files")
    args = parser.parse_args()
    failures = 0
    for uncal in args.uncal:
        try:
            download_metadata(uncal)
        except Exception as error:
            print(f"Error: {uncal}: {error}", file=sys.stderr)
            failures += 1
    return int(failures > 0)


if __name__ == "__main__":
    sys.exit(main())

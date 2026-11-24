#!/usr/bin/env python3
"""Download JWST UNCAL FITS products for one or more programs from MAST or AWS.

Run as, for example:
    python MAST_download_PID.py 04125 AWS
    python MAST_download_PID.py 04125 MAST
    python MAST_download_PID.py 04125 AWS --instrument miri
    python MAST_download_PID.py 04125 AWS --instrument all
    python MAST_download_PID.py 04125 03843 AWS
    python MAST_download_PID.py --source MAST
"""

from nirspec_ninja.archive import common, transfers
import argparse
import re
import sys
import tempfile
import time
from pathlib import Path

import numpy as np
from astroquery.mast import Observations
from astropy.table import Table, unique, vstack
from nirspec_ninja.archive.download_metrics import DownloadMetrics


LEGGOS_PIDS = [4125, 3843]
DOWNLOAD_SOURCE = "AWS"
OBS_COLLECTION = "JWST"
DEFAULT_INSTRUMENT = "NIRSPEC"
PRODUCT_LIST_CHUNK_SIZE = 100
DOWNLOAD_CHUNK_SIZE = 100


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Query MAST and download every raw level-1 UNCAL FITS product "
            "for one or more JWST programs from MAST or public AWS S3."
        )
    )
    parser.add_argument("pids", nargs="*", metavar="PID",
                        help="program IDs, optionally followed by AWS or MAST; defaults to LEGGOS_PIDS")
    parser.add_argument(
        "--source",
        dest="download_source",
        type=str.upper,
        choices=("MAST", "AWS"),
        help=f"download source (default configured in script: {DOWNLOAD_SOURCE})",
    )
    parser.add_argument(
        "--instrument",
        default=DEFAULT_INSTRUMENT,
        help=(
            "JWST instrument to select (default: nirspec). Instrument names are "
            "case-insensitive; use 'all' to disable instrument filtering."
        ),
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=Path("Rawdata"),
        help="root for PID and target directories (default: Rawdata)",
    )
    args = parser.parse_intermixed_args(argv)
    if args.pids and args.pids[-1].upper() in {"MAST", "AWS"}:
        positional_source = args.pids.pop().upper()
        if args.download_source is not None and args.download_source != positional_source:
            parser.error("trailing download source conflicts with --source")
        args.download_source = positional_source
    try:
        args.pids = list(dict.fromkeys(int(value) for value in args.pids))
    except ValueError:
        parser.error("PIDs must be positive integers")
    if any(value <= 0 for value in args.pids):
        parser.error("PIDs must be positive integers")
    return args


def instrument_criterion(instrument):
    """Return a MAST instrument_name filter, or None for all instruments."""
    name = instrument.strip()
    if not name:
        raise ValueError("Instrument must be a name or 'all'")
    if name.casefold() == "all":
        return None

    name = name.upper()
    # MAST's JWST instrument_name values include observing modes after a slash
    # (for example NIRSPEC/MSA). A bare instrument therefore needs a wildcard.
    if "/" not in name and "*" not in name and "?" not in name:
        name += "*"
    return name


def observation_query_criteria(pid, instrument):
    """Build criteria for the MAST observation query."""
    criteria = {
        "obs_collection": OBS_COLLECTION,
        "proposal_id": str(pid),
    }
    instrument_name = instrument_criterion(instrument)
    if instrument_name is not None:
        criteria["instrument_name"] = instrument_name
    return criteria


def choose_pids(pids):
    """Return supplied IDs or interactively confirm the LEGGOS defaults."""
    if pids:
        if any(pid <= 0 for pid in pids):
            raise ValueError("PIDs must be positive integers")
        return list(dict.fromkeys(pids))

    answer = input(
        f"No PIDs were supplied. Proceed with LEGGOS PIDs {LEGGOS_PIDS}? [y/N] "
    ).strip().lower()
    if answer not in {"y", "yes"}:
        print("Cancelled; no query or download was started.")
        return None
    return list(LEGGOS_PIDS)


def table_chunks(table, size):
    yield from common.table_chunks(table, size)


def get_products(observations):
    """Retrieve product metadata in bounded requests and remove duplicates."""
    product_tables = []
    chunks = list(table_chunks(observations, PRODUCT_LIST_CHUNK_SIZE))
    for number, observation_chunk in enumerate(chunks, start=1):
        print(f"Retrieving product list {number}/{len(chunks)} ...", flush=True)
        product_tables.append(Observations.get_product_list(observation_chunk))

    if not product_tables:
        return Table()
    products = vstack(product_tables, metadata_conflicts="silent")
    key = "dataURI" if "dataURI" in products.colnames else "productFilename"
    return unique(products, keys=key)


def configure_download_source(source_override):
    """Select direct MAST downloads or the public AWS dataset."""
    source = source_override or DOWNLOAD_SOURCE
    if source == "MAST":
        common.configure_cloud_dataset(Observations, "MAST")
    elif source == "AWS":
        try:
            common.configure_cloud_dataset(Observations, "AWS", profile=None)
        except ModuleNotFoundError as exc:
            raise RuntimeError("AWS downloads require the boto3 package") from exc
    else:
        raise ValueError("Download source must be MAST or AWS")
    print(f"Download source: {source}", flush=True)
    return source


def is_level_one(value):
    if np.ma.is_masked(value):
        return False
    try:
        return float(value) == 1.0
    except (TypeError, ValueError):
        return False


def select_raw_fits(products):
    """Select calibration-level 1 products whose names end in uncal.fits."""
    required = {"calib_level", "productFilename"}
    missing = required.difference(products.colnames)
    if missing:
        raise RuntimeError(
            "MAST's product response lacks required column(s): "
            + ", ".join(sorted(missing))
        )
    mask = [
        is_level_one(level) and str(filename).lower().endswith("uncal.fits")
        for level, filename in zip(
            products["calib_level"], products["productFilename"]
        )
    ]
    selected = products[np.asarray(mask, dtype=bool)]
    if len(selected):
        selected.sort("productFilename")
    return selected


def normalize_id(value):
    """Normalize MAST integer-like identifiers for table joins."""
    if np.ma.is_masked(value):
        return ""
    text = str(value).strip()
    return text[:-2] if text.endswith(".0") else text


def add_target_names(products, observations):
    """Attach each product's observation target name via its parent obsid."""
    required_observation_columns = {"obsid", "target_name"}
    missing = required_observation_columns.difference(observations.colnames)
    if missing:
        raise RuntimeError(
            "MAST's observation response lacks required column(s): "
            + ", ".join(sorted(missing))
        )
    if "parent_obsid" not in products.colnames:
        raise RuntimeError("MAST's product response lacks parent_obsid")

    target_by_obsid = {}
    for obsid, target_name in zip(
        observations["obsid"], observations["target_name"]
    ):
        key = normalize_id(obsid)
        target = "" if np.ma.is_masked(target_name) else str(target_name).strip()
        if key and target:
            target_by_obsid[key] = target

    targets = []
    for parent_obsid in products["parent_obsid"]:
        key = normalize_id(parent_obsid)
        try:
            targets.append(target_by_obsid[key])
        except KeyError as exc:
            raise RuntimeError(
                f"Could not match product parent_obsid {key!r} to a target name"
            ) from exc
    products = products.copy()
    products["target_name"] = targets
    return products


def safe_target_directory_name(target_name):
    """Make a MAST target name safe as one portable directory component."""
    safe_name = re.sub(r"[^A-Za-z0-9._+-]+", "_", target_name.strip())
    safe_name = safe_name.strip("._")
    if not safe_name:
        raise ValueError(f"Target name {target_name!r} has no usable characters")
    return safe_name


def flatten_manifest(manifest, output_dir):
    """Move files named by an Astroquery manifest into one flat directory."""
    if manifest is None or "Local Path" not in manifest.colnames:
        return 0

    moved = 0
    for value in manifest["Local Path"]:
        if np.ma.is_masked(value):
            continue
        source = Path(str(value))
        if not source.is_file():
            continue
        destination = output_dir / source.name
        if destination.exists():
            if source.resolve() != destination.resolve():
                if source.stat().st_size != destination.stat().st_size:
                    raise FileExistsError(
                        f"Refusing to replace existing file {destination}"
                    )
                source.unlink()
            continue
        source.replace(destination)
        moved += 1

    return moved


def download_flat(products, pid_dir, source):
    """Download batches grouped by target and flatten within each target."""
    output_dir = pid_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    moved = 0
    metrics = DownloadMetrics(pid_dir / "transfer_metrics", "program", workers=1)
    target_names = sorted({str(value) for value in products["target_name"]})
    for target_name in target_names:
        target_products = products[
            np.asarray(
                [str(value) == target_name for value in products["target_name"]],
                dtype=bool,
            )
        ]
        target_dir = pid_dir / safe_target_directory_name(target_name)
        target_dir.mkdir(parents=True, exist_ok=True)
        chunks = list(table_chunks(target_products, 1))
        for number, product_chunk in enumerate(chunks, start=1):
            filename = str(product_chunk["productFilename"][0])
            destination = target_dir / filename
            expected = (product_chunk["size"][0]
                        if "size" in product_chunk.colnames else None)
            if (expected is not None and not np.ma.is_masked(expected)
                    and int(expected) > 0 and destination.is_file()
                    and destination.stat().st_size == int(expected)):
                metrics.record(source, filename, "download", time.monotonic(),
                               status="skipped")
                continue
            print(
                f"Downloading {target_name}: batch {number}/{len(chunks)} "
                f"({len(product_chunk)} files) ...",
                flush=True,
            )
            with tempfile.TemporaryDirectory(
                dir=target_dir, prefix=".download_staging_"
            ) as staging_dir:
                filename = str(product_chunk["productFilename"][0])
                started = time.monotonic()
                try:
                    manifest = transfers.download_products(Observations, product_chunk, staging_dir, source)
                    paths = [] if manifest is None else [
                        Path(str(value)) for value in manifest["Local Path"]
                        if not np.ma.is_masked(value) and Path(str(value)).is_file()
                    ]
                    if not paths:
                        raise RuntimeError(f"Download returned no local file: {manifest}")
                    size = sum(path.stat().st_size for path in paths)
                    if (expected is not None and not np.ma.is_masked(expected)
                            and int(expected) > 0 and size != int(expected)):
                        raise RuntimeError(f"Incomplete file: {filename}; bytes={size}, expected={expected}")
                    moved += flatten_manifest(manifest, target_dir)
                    metrics.record(source, filename, "download", started, byte_count=size)
                except Exception as exc:
                    metrics.record(source, filename, "download", started, error=exc)
                    raise
    return moved


def download_pid(pid, source, instrument, output_root):
    """Query and download one program into its own PID/target tree."""
    pid_dir = output_root / f"PID{pid:05d}"
    criteria = observation_query_criteria(pid, instrument)
    instrument_label = criteria.get("instrument_name", "all instruments")
    print(
        f"Querying MAST for {OBS_COLLECTION} program {pid}, "
        f"instrument {instrument_label} ...",
        flush=True,
    )
    observations = Observations.query_criteria(**criteria)
    print(f"Found {len(observations)} observations.", flush=True)
    if not len(observations):
        print(f"No observations matched PID {pid}; nothing to download.")
        return 0

    products = get_products(observations)
    raw_fits = select_raw_fits(products)
    print(
        f"Selected {len(raw_fits)} raw level-1 FITS files from "
        f"{len(products)} products.",
        flush=True,
    )
    if not len(raw_fits):
        print("No matching files were found; nothing to download.")
        return 0

    raw_fits = add_target_names(raw_fits, observations)
    moved = download_flat(raw_fits, pid_dir, source)
    print(
        f"Finished. {moved} new files were placed by target beneath "
        f"{pid_dir.resolve()}"
    )
    return 0


def main(argv=None):
    args = parse_args(argv)
    pids = choose_pids(args.pids)
    if pids is None:
        return 0
    source = configure_download_source(args.download_source)
    for pid in pids:
        download_pid(pid, source, args.instrument, args.output_root)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except (EOFError, KeyboardInterrupt):
        print("\nCancelled.", file=sys.stderr)
        raise SystemExit(130)

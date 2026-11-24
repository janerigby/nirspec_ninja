#!/usr/bin/env python3
"""Query large NIRSpec datasets with parallel, restart-safe transfers."""

# Original from Rosalia O'Brien. jrigby added bells and whistles w GPT-5.6-Sol

from nirspec_ninja.archive import common, transfers
import argparse
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from contextlib import contextmanager
import csv
from datetime import datetime
import hashlib
import re
import signal
import subprocess
import threading
import time
from pathlib import Path

from astroquery.mast import Observations
from astroquery.mast.observations import ObservationsClass
from astropy.coordinates import SkyCoord
from astropy.table import Table, unique, vstack
import astropy.units as u
from astropy.time import Time
import numpy as np
from nirspec_ninja.archive.download_metrics import DownloadMetrics

transfer_metrics = None
DIRECT_MAST_DOWNLOADS = False


# --------------------------- User settings ---------------------------

filter_choices = [
    "CLEAR;prism",
    "F070LP;G140M",
    "F100LP;G140M",
    "F070LP;G140H",
    "F100LP;G140H",
    "F170LP;G235M",
    "F170LP;G235H",
    "F290LP;G395M",
    "F290LP;G395H",
]

# Scientific selection settings.
myfilter_grating_combo = filter_choices[1]  # set to [0] for prism!!!
minexptime = 2000
earliest_release_time = 2022.5
mingal_lat = 20.0 
TARGET_TYPE = "EITHER"   # Choose "FIXED", "MOVING", or "EITHER".

# Run mode:
#   "QUERY" performs the MAST observation query, selects products, saves the
#   final selection, and then downloads the files in the same run.
#   "DOWNLOAD_ONLY" skips all MAST query and product-list calls and reuses the
#   saved final selection.
#RUN_MODE = "DOWNLOAD_ONLY"
RUN_MODE = 'QUERY'

# Choose "MAST" to create curl scripts or "AWS" to download directly from S3.
# With RUN_MODE="QUERY", AWS discovery and downloading happen in one run.
DOWNLOAD_SOURCE = "AWS"
AWS_PROFILE = None

# Move AWS files out of Astroquery's mastDownload/JWST/<observation>/ tree and
# place them directly in downloads/prism_observations/.
FLATTEN_AWS_DOWNLOADS = True

# Exclude products whose case-insensitive MAST `filters` metadata contains MIRROR.
EXCLUDE_MIRROR_GRATING = True

# The final filtered product table is written atomically and can be reused by
# DOWNLOAD_ONLY mode. Keep the same path when switching download sources.
FINAL_PRODUCT_MANIFEST = "checkpoints/final_products.ecsv"

# One summary CSV is written per run beneath a timestamped directory.
RUN_METRICS_DIR = "run_metrics"

# Final failures from an AWS run are written atomically here.  The complete
# product manifest remains the authoritative list used for restart scans.
FAILED_DOWNLOAD_LOG = "checkpoints/aws_failed_downloads.csv"

# Cisco Secure Client command used to record VPN status in run metrics.
VPN_STATS_COMMAND = "/opt/cisco/secureclient/bin/vpn"

# Batching and product-list checkpoint settings.
PRODUCT_LIST_CHUNK_SIZE = 100
DOWNLOAD_CHUNK_SIZE = 100
PRODUCT_LIST_CHECKPOINT_DIR = "checkpoints/product_lists"

# Query concurrency and AWS download concurrency are controlled separately.
DEFAULT_WORKERS = 3
PRODUCT_QUERY_WORKERS = 2
query_benchmark = {}
AWS_PROGRESS_EVERY_FILES = 25
AWS_DOWNLOAD_ATTEMPTS = 5
AWS_RETRY_DELAYS = (5, 10, 20, 40)

# Progress-monitoring settings. Use None to disable the timeout or 0 to disable
# heartbeat messages.
API_TIMEOUT_SECONDS = 15 * 60
HEARTBEAT_SECONDS = 60

# ---------------------------------------------------------------------


def log(message):
    """Print a timestamped status message immediately."""
    timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
    print(f"[{timestamp}] {message}", flush=True)


@contextmanager
def progress_monitor(stage, heartbeat_seconds=None, timeout_seconds=None):
    """Report that a blocking stage is alive and impose a hard time limit."""
    if heartbeat_seconds is None:
        heartbeat_seconds = HEARTBEAT_SECONDS
    if timeout_seconds is None:
        timeout_seconds = API_TIMEOUT_SECONDS
    start = time.monotonic()
    stopped = threading.Event()

    def heartbeat():
        while not stopped.wait(heartbeat_seconds):
            elapsed = time.monotonic() - start
            log(f"Still working: {stage} ({elapsed:.0f} s elapsed)")

    heartbeat_thread = None
    if heartbeat_seconds and heartbeat_seconds > 0:
        heartbeat_thread = threading.Thread(target=heartbeat, daemon=True)
        heartbeat_thread.start()

    previous_handler = None
    if timeout_seconds is not None:
        if threading.current_thread() is not threading.main_thread():
            raise RuntimeError("The progress timeout must run in the main thread")

        def timed_out(signum, frame):
            raise TimeoutError(
                f"{stage} exceeded the {timeout_seconds}-second timeout"
            )

        previous_handler = signal.signal(signal.SIGALRM, timed_out)
        signal.setitimer(signal.ITIMER_REAL, timeout_seconds)

    succeeded = False
    try:
        yield
        succeeded = True
    finally:
        if timeout_seconds is not None:
            signal.setitimer(signal.ITIMER_REAL, 0)
            signal.signal(signal.SIGALRM, previous_handler)
        stopped.set()
        if heartbeat_thread is not None:
            heartbeat_thread.join(timeout=1)
        elapsed = time.monotonic() - start
        outcome = "Finished" if succeeded else "Stopped with an error"
        log(f"{outcome}: {stage} ({elapsed:.1f} s)")


def log_count(stage, table):
    log(f"{stage}: {len(table):,} rows")


def chunks(table, chunk_size):
    """Yield consecutive table slices of at most chunk_size rows."""
    if chunk_size <= 0:
        raise ValueError("Chunk sizes must be positive integers")
    yield from common.table_chunks(table, chunk_size)


def configure_download_source(source_override=None):
    """Configure Astroquery for MAST curl scripts or direct AWS downloads."""
    source = (source_override or DOWNLOAD_SOURCE).upper()
    if source == "MAST":
        common.configure_cloud_dataset(Observations, "MAST")
        log("Download source=MAST (" +
            ("direct downloads" if DIRECT_MAST_DOWNLOADS else "creating curl scripts") + ")")
    elif source == "AWS":
        try:
            common.configure_cloud_dataset(Observations, "AWS", profile=AWS_PROFILE)
        except ModuleNotFoundError as exc:
            raise RuntimeError(
                "DOWNLOAD_SOURCE='AWS' requires the optional boto3 package"
            ) from exc
        log("Download source=AWS S3 (downloading files directly; no MAST fallback)")
    else:
        raise ValueError("DOWNLOAD_SOURCE must be either 'MAST' or 'AWS'")
    return source


def parse_arguments(argv=None):
    """Parse the source override and independent query/download worker counts."""
    parser = argparse.ArgumentParser(
        description="Query MAST and create curl scripts or download from AWS."
    )
    parser.add_argument(
        "download_source",
        nargs="?",
        type=str.upper,
        choices=("MAST", "AWS"),
        help="override the DOWNLOAD_SOURCE value configured in the script",
    )
    parser.add_argument(
        "-j", "--worker", "--workers",
        dest="workers",
        type=int,
        default=DEFAULT_WORKERS,
        help=("number of concurrent AWS download workers "
              f"(default: {DEFAULT_WORKERS})"),
    )
    parser.add_argument(
        "-q", "--query-workers", type=int, default=PRODUCT_QUERY_WORKERS,
        help=f"concurrent MAST product-list queries (default: {PRODUCT_QUERY_WORKERS})",
    )
    parser.add_argument("--direct-mast", action="store_true",
                        help="download MAST files directly instead of generating curl scripts")
    arguments = parser.parse_args(argv)
    if arguments.workers <= 0:
        parser.error("--worker must be a positive integer")
    if arguments.query_workers <= 0:
        parser.error("--query-workers must be a positive integer")
    return arguments


def product_checkpoint_path(observation_chunk, chunk_number):
    """Return a checkpoint path fingerprinted by the chunk's MAST obsids."""
    obsids = "\0".join(str(obsid) for obsid in observation_chunk["obsid"])
    fingerprint = hashlib.sha256(obsids.encode("utf-8")).hexdigest()[:16]
    checkpoint_dir = Path(PRODUCT_LIST_CHECKPOINT_DIR)
    return checkpoint_dir / f"chunk_{chunk_number:04d}_{fingerprint}.ecsv"


def write_product_checkpoint(products, checkpoint_path):
    """Atomically save a completed product-list response."""
    checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = checkpoint_path.with_suffix(".tmp.ecsv")
    products.write(temporary_path, format="ascii.ecsv", overwrite=True)
    temporary_path.replace(checkpoint_path)


def format_duration(seconds):
    """Return an elapsed time that remains readable for multi-hour queries."""
    if seconds >= 3600:
        return f"{seconds / 3600:.2f} hours"
    if seconds >= 60:
        return f"{seconds / 60:.1f} minutes"
    return f"{seconds:.1f} seconds"


def fetch_product_chunk(observation_chunk, checkpoint_path):
    """Fetch and atomically checkpoint one product-list chunk."""
    # Give each worker its own Astroquery client/session rather than sharing
    # the module-level Observations singleton between threads.
    client = ObservationsClass()
    products = client.get_product_list(observation_chunk)
    write_product_checkpoint(products, checkpoint_path)
    return products


def get_product_list_in_chunks(observations, workers=PRODUCT_QUERY_WORKERS):
    """Retrieve missing product chunks concurrently and combine in order."""
    if workers <= 0:
        raise ValueError("Product-query worker count must be positive")

    observation_chunks = list(chunks(
        observations, PRODUCT_LIST_CHUNK_SIZE
    ))
    total_chunks = len(observation_chunks)
    product_tables = [None] * total_chunks
    pending = []
    reused_chunks = 0

    query_benchmark.update({
        "queried_observations": len(observations),
        "product_chunk_count": total_chunks,
        "product_chunk_size": PRODUCT_LIST_CHUNK_SIZE,
        "product_chunk_sizes": [len(chunk) for chunk in observation_chunks],
        "query_workers": workers,
    })

    for index, observation_chunk in enumerate(observation_chunks):
        checkpoint_path = product_checkpoint_path(
            observation_chunk, index + 1
        )
        if checkpoint_path.is_file():
            product_tables[index] = Table.read(
                checkpoint_path, format="ascii.ecsv"
            )
            reused_chunks += 1
        else:
            pending.append((index, observation_chunk, checkpoint_path))

    log(
        f"Retrieving products in {total_chunks:,} chunks of at most "
        f"{PRODUCT_LIST_CHUNK_SIZE} observations with "
        f"{workers} concurrent MAST queries"
    )
    log(
        f"Product-list checkpoints: {reused_chunks:,} reused; "
        f"{len(pending):,} pending"
    )

    completed = 0
    executor = ThreadPoolExecutor(max_workers=workers)
    futures = {
        executor.submit(
            fetch_product_chunk, observation_chunk, checkpoint_path
        ): (index, len(observation_chunk))
        for index, observation_chunk, checkpoint_path in pending
    }
    try:
        while futures:
            done, _ = wait(
                futures,
                timeout=HEARTBEAT_SECONDS,
                return_when=FIRST_COMPLETED,
            )
            if not done:
                log(
                    f"Still retrieving product lists: {len(futures):,} "
                    "chunks in flight or queued"
                )
                continue
            for future in done:
                index, observation_count = futures.pop(future)
                # Propagate the first MAST or checkpoint error immediately.
                product_tables[index] = future.result()
                completed += 1
                log(
                    f"Completed product chunk {index + 1:,}/{total_chunks:,} "
                    f"({observation_count:,} observations); "
                    f"{completed:,}/{len(pending):,} pending chunks finished"
                )
    except BaseException:
        for future in futures:
            future.cancel()
        executor.shutdown(wait=False, cancel_futures=True)
        raise
    else:
        executor.shutdown(wait=True)

    if not product_tables:
        query_benchmark.update({
            "product_chunks_reused": reused_chunks,
            "product_chunks_fetched": completed,
            "candidate_product_rows": 0,
            "unique_candidate_products": 0,
        })
        return Table()
    products = vstack(product_tables, metadata_conflicts="silent")
    candidate_product_rows = len(products)
    unique_key = (
        "dataURI" if "dataURI" in products.colnames else "productFilename"
    )
    products = unique(products, keys=unique_key)
    query_benchmark.update({
        "product_chunks_reused": reused_chunks,
        "product_chunks_fetched": completed,
        "candidate_product_rows": candidate_product_rows,
        "unique_candidate_products": len(products),
    })
    return products


def query_and_select_products(workers=PRODUCT_QUERY_WORKERS):
    """Run scientific selection and report query-stage measurements."""
    query_benchmark.clear()
    start = time.monotonic()
    products = _query_and_select_products(workers)
    elapsed = time.monotonic() - start
    returned_files = 0 if products is None else len(products)
    query_benchmark.update({
        "query_elapsed_seconds": elapsed,
        "query_elapsed_human": format_duration(elapsed),
        "returned_files": returned_files,
    })

    observations = query_benchmark.get("queried_observations", 0)
    chunk_count = query_benchmark.get("product_chunk_count", 0)
    chunk_limit = query_benchmark.get(
        "product_chunk_size", PRODUCT_LIST_CHUNK_SIZE
    )
    raw_candidates = query_benchmark.get("candidate_product_rows", 0)
    unique_candidates = query_benchmark.get("unique_candidate_products", 0)
    reused = query_benchmark.get("product_chunks_reused", 0)
    fetched = query_benchmark.get("product_chunks_fetched", 0)
    log(
        f"Query benchmark: total time={format_duration(elapsed)}; "
        f"queried {observations:,} observations in {chunk_count:,} product-list "
        f"chunks of at most {chunk_limit:,} observations "
        f"({fetched:,} fetched, {reused:,} reused); MAST returned "
        f"{raw_candidates:,} candidate product rows "
        f"({unique_candidates:,} unique); final selection returned "
        f"{returned_files:,} files"
    )
    return products


def remove_empty_subdirectories(root):
    """Remove only empty directories below root, deepest first."""
    directories = sorted(
        (path for path in root.rglob("*") if path.is_dir()),
        key=lambda path: len(path.parts),
        reverse=True,
    )
    for directory in directories:
        try:
            directory.rmdir()
        except OSError:
            pass


def flatten_aws_manifest(manifest, destination_dir, overwrite=False):
    """Move downloaded files to one directory and update manifest paths."""
    if manifest is None or "Local Path" not in manifest.colnames:
        return 0
    destination_dir.mkdir(parents=True, exist_ok=True)
    moved = 0
    for row_number, value in enumerate(manifest["Local Path"]):
        if np.ma.is_masked(value):
            continue
        source_path = Path(str(value))
        if not source_path.is_file():
            continue
        destination_path = destination_dir / source_path.name
        if source_path.resolve() == destination_path.resolve():
            continue
        if destination_path.exists() and not overwrite:
            raise FileExistsError(
                f"Cannot flatten {source_path}: {destination_path} already exists"
            )
        source_path.replace(destination_path)
        manifest["Local Path"][row_number] = str(destination_path)
        moved += 1
    return moved


def write_failed_download_log(failures):
    """Atomically record the files that exhausted all AWS retry attempts."""
    log_path = Path(FAILED_DOWNLOAD_LOG)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = log_path.with_suffix(".tmp.csv")
    with temporary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream)
        writer.writerow(("productFilename", "exceptionType", "exception"))
        for filename, exc in failures:
            writer.writerow((filename, type(exc).__name__, str(exc)))
    temporary_path.replace(log_path)
    log(f"Saved AWS failure log: {log_path} ({len(failures):,} files)")


def product_filename(product):
    """Return the local basename for one MAST product row."""
    for column in ("productFilename", "dataURI"):
        if column in product.colnames:
            value = str(product[column][0])
            if value and value != "--":
                return Path(value).name
    raise RuntimeError("Product has neither productFilename nor dataURI")


def product_expected_size(product):
    """Return MAST's expected byte count, or None when it is unavailable."""
    if "size" not in product.colnames or np.ma.is_masked(product["size"][0]):
        return None
    try:
        size = int(product["size"][0])
    except (TypeError, ValueError, OverflowError):
        return None
    return size if size > 0 else None


def completed_aws_path(product, destination_dir):
    """Recognize a complete flattened file left by this or an earlier run."""
    path = Path(destination_dir) / product_filename(product)
    if not path.is_file():
        return None
    actual_size = path.stat().st_size
    expected_size = product_expected_size(product)
    if actual_size > 0 and (expected_size is None or actual_size == expected_size):
        return path
    log(f"Will retry incomplete AWS file {path.name}: local={actual_size:,} "
        f"bytes, expected={expected_size if expected_size is not None else 'unknown'}")
    return None


def download_one_aws_product(product, destination_dir, source="AWS"):
    """Download and flatten one product; the final file is its checkpoint.

    Every attempt creates a fresh Astroquery client and session.  Cleanup is
    restricted to this product's uniquely fingerprinted staging directory so
    one worker cannot remove a directory another worker is about to use.
    """
    filename = product_filename(product)
    identity = str(product["dataURI"][0]) if "dataURI" in product.colnames else filename
    fingerprint = hashlib.sha256(identity.encode("utf-8")).hexdigest()[:16]
    staging_dir = Path(destination_dir) / ".aws_staging" / fingerprint
    def operation():
        staging_dir.mkdir(parents=True, exist_ok=True)
        client = ObservationsClass()
        if source == "AWS":
            client.enable_cloud_dataset(provider="AWS", profile=AWS_PROFILE)
        manifest = transfers.download_products(client, product, staging_dir, source)
        # Verify the staged file before replacing an existing final file.
        expected = product_expected_size(product)
        staged_paths = [] if manifest is None or "Local Path" not in manifest.colnames else [
            Path(str(value)) for value in manifest["Local Path"]
            if not np.ma.is_masked(value) and Path(str(value)).name == filename
            and Path(str(value)).is_file()
        ]
        if (not staged_paths or any(
                path.stat().st_size <= 0 or
                (expected is not None and path.stat().st_size != expected)
                for path in staged_paths)):
            raise RuntimeError(f"{source} returned no complete staged file for {filename}: {manifest}")
        # Path.replace is atomic on the destination filesystem.  Thus an
        # old truncated file remains usable as a restart marker until its
        # replacement has downloaded completely into staging.
        flatten_aws_manifest(manifest, Path(destination_dir), overwrite=True)
        completed = completed_aws_path(product, destination_dir)
        if completed is None:
            raise RuntimeError(
                f"{source} returned no complete local file for {filename}: {manifest}"
            )
        remove_empty_subdirectories(staging_dir)
        try:
            staging_dir.rmdir()
        except OSError:
            pass
        return completed
    def success(completed, started, attempt):
        if transfer_metrics is not None:
            transfer_metrics.record(source, filename, "download", started,
                                    attempt, completed.stat().st_size)
    def failure(error, started, attempt):
        if transfer_metrics is not None:
            transfer_metrics.record(source, filename, "download", started,
                                    attempt, error=error)
    return transfers.retry_transfer(operation, attempts=AWS_DOWNLOAD_ATTEMPTS,
        exceptions=Exception, clock=time.monotonic, sleep=time.sleep,
        delay=lambda attempt: AWS_RETRY_DELAYS[attempt - 1],
        success=success, failure=failure,
        retry_message=lambda error, attempt, seconds: log(
            f"{source} download attempt {attempt}/{AWS_DOWNLOAD_ATTEMPTS} failed for {filename}: {type(error).__name__}: {error}; retrying in {seconds} s"),
        exhausted=lambda error: RuntimeError(
            f"{source} download failed after {AWS_DOWNLOAD_ATTEMPTS} attempts for {filename}: {type(error).__name__}: {error}"))


def download_products_from_aws_parallel(products, download_dir, workers, source="AWS"):
    """Download AWS products concurrently with file-level restart safety."""
    if workers <= 0:
        raise ValueError("workers must be a positive integer")
    destination_dir = Path(download_dir)
    destination_dir.mkdir(parents=True, exist_ok=True)

    names = [product_filename(products[index:index + 1])
             for index in range(len(products))]
    if len(names) != len(set(names)):
        duplicates = sorted({name for name in names if names.count(name) > 1})
        raise RuntimeError("Cannot flatten products with duplicate filenames: "
                           + ", ".join(duplicates[:5]))

    completed_paths = []
    pending = []
    for index in range(len(products)):
        product = products[index:index + 1]
        path = completed_aws_path(product, destination_dir)
        if path is None:
            pending.append(product)
        else:
            completed_paths.append(path)
            if transfer_metrics is not None:
                transfer_metrics.record(source, path.name, "download", time.monotonic(),
                                        status="skipped")

    total = len(products)
    log(f"{source} restart scan: {len(completed_paths):,}/{total:,} complete; "
        f"{len(pending):,} to download with {workers} workers")
    failures = []
    newly_downloaded_paths = []
    start = time.monotonic()
    last_report = start
    executor = ThreadPoolExecutor(max_workers=workers,
                                  thread_name_prefix="aws-download")
    futures = {
        executor.submit(download_one_aws_product, product, destination_dir, source): product
        for product in pending
    }
    try:
        unfinished = set(futures)
        while unfinished:
            done, unfinished = wait(unfinished, timeout=1, return_when=FIRST_COMPLETED)
            for future in done:
                product = futures[future]
                try:
                    path = future.result()
                    completed_paths.append(path)
                    newly_downloaded_paths.append(path)
                except Exception as exc:
                    failures.append((product_filename(product), exc))
            now = time.monotonic()
            newly_handled = len(completed_paths) + len(failures)
            should_report = (
                done and (newly_handled == total or
                          (AWS_PROGRESS_EVERY_FILES > 0 and
                           newly_handled % AWS_PROGRESS_EVERY_FILES == 0))
            ) or (HEARTBEAT_SECONDS and now - last_report >= HEARTBEAT_SECONDS)
            if should_report:
                elapsed = now - start
                log(f"{source} progress: {newly_handled:,}/{total:,} handled "
                    f"({len(completed_paths):,} complete, {len(failures):,} failed, "
                    f"{len(unfinished):,} active/queued; {elapsed:.0f} s)")
                last_report = now
    except BaseException:
        for future in futures:
            future.cancel()
        executor.shutdown(wait=False, cancel_futures=True)
        raise
    else:
        executor.shutdown(wait=True)

    # No worker is using staging now, so recursive empty-directory cleanup is
    # safe.  During transfers, workers clean only their own staging directory.
    staging_root = destination_dir / ".aws_staging"
    if staging_root.is_dir():
        remove_empty_subdirectories(staging_root)
        try:
            staging_root.rmdir()
        except OSError:
            pass

    write_failed_download_log(failures)
    if failures:
        examples = "; ".join(
            f"{name}: {type(exc).__name__}: {exc}"
            for name, exc in failures[:3]
        )
        raise RuntimeError(
            f"{len(failures):,} of {total:,} {source} downloads failed. "
            f"Completed files are restart-safe. First failures: {examples}"
        )
    newly_downloaded = {path.resolve() for path in newly_downloaded_paths}
    return Table({
        "Local Path": [str(path) for path in completed_paths],
        "Downloaded This Run": [
            path.resolve() in newly_downloaded for path in completed_paths
        ],
    })


def download_products_in_chunks(products, download_dir, source,
                                workers=DEFAULT_WORKERS):
    """Create MAST curl scripts or directly download AWS files in chunks."""
    if source == "AWS" or DIRECT_MAST_DOWNLOADS:
        return download_products_from_aws_parallel(products, download_dir, workers, source)

    product_chunks = list(chunks(products, DOWNLOAD_CHUNK_SIZE))
    manifests = []
    total_chunks = len(product_chunks)
    action = "curl manifests" if source == "MAST" else "AWS download batches"
    log(f"Processing {total_chunks:,} {action} for at most "
        f"{DOWNLOAD_CHUNK_SIZE} files each")

    base_download_dir = Path(download_dir)
    for chunk_number, product_chunk in enumerate(product_chunks, start=1):
        if source == "AWS" and FLATTEN_AWS_DOWNLOADS:
            chunk_dir = base_download_dir
        else:
            chunk_dir = base_download_dir / f"chunk_{chunk_number:04d}"
        chunk_dir.mkdir(parents=True, exist_ok=True)
        stage_action = ("create curl manifest" if source == "MAST"
                        else "download from AWS")
        stage = (f"{stage_action} chunk {chunk_number:,}/{total_chunks:,} "
                 f"({len(product_chunk):,} files)")
        with progress_monitor(stage):
            manifest = Observations.download_products(
                product_chunk,
                curl_flag=(source == "MAST"),
                cloud_only=(source == "AWS"),
                download_dir=str(chunk_dir),
            )
        if source == "AWS" and FLATTEN_AWS_DOWNLOADS:
            moved = flatten_aws_manifest(manifest, base_download_dir)
            log(f"Flattened AWS chunk {chunk_number:,}/{total_chunks:,}: "
                f"{moved:,} files moved to {base_download_dir}")
        if manifest is not None:
            manifests.append(manifest)

    if not manifests:
        return None
    return vstack(manifests, metadata_conflicts="silent")


def mjd_to_decimalyear(table, input_col, output_col="decimal_year"):
    """Add a decimal-year column converted from an MJD column."""
    converted_time = Time(table[input_col], format="mjd")
    table[output_col] = converted_time.decimalyear
    return table


def exclude_mirror_products(products):
    """Exclude products whose case-insensitive MAST filter contains MIRROR."""
    if not EXCLUDE_MIRROR_GRATING:
        return products
    if "filters" not in products.colnames:
        raise RuntimeError(
            "Cannot exclude MIRROR products: MAST product metadata has no "
            "'filters' column"
        )
    return products[
        ["MIRROR" not in str(value).upper() for value in products["filters"]]
    ]


def filter_level1_uncal_products(products):
    """Select calibration-level 1 UNCAL products without case sensitivity."""
    required_columns = {"calib_level", "productSubGroupDescription"}
    missing_columns = required_columns - set(products.colnames)
    if missing_columns:
        raise RuntimeError(
            "Cannot select level-1 UNCAL products; missing MAST columns: "
            + ", ".join(sorted(missing_columns))
        )
    mask = [
        str(level).strip() in {"1", "1.0"}
        and str(subgroup).strip().upper() == "UNCAL"
        for level, subgroup in zip(
            products["calib_level"],
            products["productSubGroupDescription"],
        )
    ]
    return products[mask]


def local_download_totals(manifest, current_session_only=False):
    """Return unique local file count and bytes reported by a manifest."""
    if manifest is None or "Local Path" not in manifest.colnames:
        return 0, 0
    if current_session_only and "Downloaded This Run" not in manifest.colnames:
        return 0, 0
    paths = set()
    for row_number, value in enumerate(manifest["Local Path"]):
        if (current_session_only and
                not bool(manifest["Downloaded This Run"][row_number])):
            continue
        if np.ma.is_masked(value):
            continue
        path = Path(str(value))
        if path.is_file():
            paths.add(path.resolve())
    return len(paths), sum(path.stat().st_size for path in paths)


def rate(numerator, seconds):
    """Return a rate, or an empty value when it is not meaningful."""
    return numerator / seconds if seconds > 0 else ""


def write_run_metrics(metrics, run_timestamp):
    """Atomically save one run's query and download performance metrics."""
    if metrics.get("run_mode") == "QUERY":
        metrics.update(query_benchmark)
        metrics["product_chunk_sizes"] = ",".join(
            str(size) for size in query_benchmark.get("product_chunk_sizes", [])
        )
    run_stamp = run_timestamp.strftime("%Y%m%dT%H%M%S%z")
    metrics_dir = Path(RUN_METRICS_DIR) / run_stamp
    metrics_dir.mkdir(parents=True, exist_ok=True)
    metrics_path = metrics_dir / "summary.csv"
    temporary_path = metrics_dir / "summary.tmp.csv"
    with temporary_path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics))
        writer.writeheader()
        writer.writerow(metrics)
    temporary_path.replace(metrics_path)
    log(f"Run metrics written to {metrics_path}")
    return metrics_path


def get_vpn_status():
    """Return active, inactive, unknown, or unavailable for Cisco VPN."""
    try:
        result = subprocess.run(
            [VPN_STATS_COMMAND, "stats"],
            capture_output=True,
            text=True,
            timeout=10,
            check=False,
        )
    except (FileNotFoundError, PermissionError, subprocess.TimeoutExpired):
        return "unavailable"

    output = f"{result.stdout}\n{result.stderr}"
    match = re.search(
        r"^\s*Connection State:\s*(.*?)\s*$",
        output,
        flags=re.IGNORECASE | re.MULTILINE,
    )
    state = match.group(1).strip().lower() if match else ""
    if state == "connected":
        return "active"
    if state == "disconnected":
        return "inactive"
    if state in {"not available", "unknown"}:
        return "unknown"

    state_matches = re.findall(
        r">>\s*state:\s*([^\r\n]+)", output, flags=re.IGNORECASE
    )
    for value in reversed(state_matches):
        normalized = value.strip().lower()
        if normalized == "connected":
            return "active"
        if normalized == "disconnected":
            return "inactive"
    return "unknown" if output.strip() else "unavailable"


def write_final_product_manifest(products):
    """Atomically save the exact filtered product selection for later reuse."""
    manifest_path = Path(FINAL_PRODUCT_MANIFEST)
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    temporary_path = manifest_path.with_suffix(".tmp.ecsv")
    products = products.copy()
    products.meta["filter_grating"] = myfilter_grating_combo
    products.meta["minimum_exposure_seconds"] = minexptime
    products.meta["earliest_release_decimal_year"] = earliest_release_time
    products.meta["target_type"] = TARGET_TYPE.upper()
    products.write(temporary_path, format="ascii.ecsv", overwrite=True)
    temporary_path.replace(manifest_path)
    log(f"Saved final product selection: {manifest_path} "
        f"({len(products):,} files)")


def read_final_product_manifest():
    """Load a previously saved final selection without contacting MAST."""
    manifest_path = Path(FINAL_PRODUCT_MANIFEST)
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"No final product manifest found at {manifest_path}. "
            "Run once with RUN_MODE='QUERY' before using DOWNLOAD_ONLY."
        )
    products = Table.read(manifest_path, format="ascii.ecsv")
    log(f"Loaded final product selection: {manifest_path} "
        f"({len(products):,} files)")
    if products.meta:
        log("Saved selection settings: "
            f"filter/grating={products.meta.get('filter_grating', 'unknown')}; "
            f"minimum exposure="
            f"{products.meta.get('minimum_exposure_seconds', 'unknown')} s; "
            f"release year>"
            f"{products.meta.get('earliest_release_decimal_year', 'unknown')}; "
            f"target type={products.meta.get('target_type', 'unknown')}")
    return products


def _query_and_select_products(workers):
    """Run the scientific selection and persist its final product table."""
    target_type = TARGET_TYPE.upper()
    if target_type not in {"FIXED", "MOVING", "EITHER"}:
        raise ValueError(
            "TARGET_TYPE must be 'FIXED', 'MOVING', or 'EITHER'"
        )
    criteria = {
        "instrument_name": ["NIRSPEC/MSA", "NIRSPEC/IFU"],
        "provenance_name": ["CALJWST"],
        "filters": [myfilter_grating_combo.upper()],
        "dataRights": "PUBLIC",
    }
    if target_type != "EITHER":
        criteria["mtFlag"] = target_type == "MOVING"

    with progress_monitor("query MAST observations"):
        observations = Observations.query_criteria(**criteria)
    log_count("MAST query result", observations)

    observations = observations[observations["t_exptime"] > minexptime]
    log_count("After exposure-time filter", observations)
    observations = mjd_to_decimalyear(
        observations, "t_obs_release", output_col="t_obs_release_decimal_year"
    )
    observations = observations[
        observations["t_obs_release_decimal_year"] > earliest_release_time
    ]
    log_count("After release-date filter", observations)

    coordinates = SkyCoord(
        ra=observations["s_ra"] * u.deg,
        dec=observations["s_dec"] * u.deg,
        obstime=Time(observations["t_min"], format="mjd"),
    )
    observations = observations[np.abs(coordinates.galactic.b) > mingal_lat * u.deg]
    log_count("After Galactic-latitude filter", observations)
    if not len(observations):
        log("No observations remain; no final product manifest created")
        return None

    products = get_product_list_in_chunks(observations, workers)
    log_count("Product-list result", products)
    if "prism" in myfilter_grating_combo.lower():
        target_detector = "nrs1"
        products = products[
            [target_detector in str(value).lower() for value in products["obs_id"]]
        ]
        log_count(f"After {target_detector.upper()} detector filter", products)
    else:
        log("Detector filter skipped for non-prism configuration")

    products = filter_level1_uncal_products(products)
    log_count("After product filters", products)
    if EXCLUDE_MIRROR_GRATING:
        products = exclude_mirror_products(products)
        log_count("After MIRROR grating exclusion", products)
    write_final_product_manifest(products)
    return products


def main(download_source_override=None, workers=DEFAULT_WORKERS,
         query_workers=PRODUCT_QUERY_WORKERS, direct_mast=False):
    global transfer_metrics, DIRECT_MAST_DOWNLOADS
    DIRECT_MAST_DOWNLOADS = direct_mast
    transfer_metrics = DownloadMetrics("transfer_metrics", "science_query", workers)
    run_timestamp = datetime.now().astimezone()
    start = time.monotonic()
    vpn_status_start = get_vpn_status()
    log(f"VPN status at start={vpn_status_start}")
    mode = RUN_MODE.upper()
    if mode not in {"QUERY", "DOWNLOAD_ONLY"}:
        raise ValueError("RUN_MODE must be either 'QUERY' or 'DOWNLOAD_ONLY'")

    log(f"Large NIRSpec dataset query started; mode={mode}")
    source = configure_download_source(download_source_override)
    if download_source_override:
        log(f"Runtime download-source override={source}")
    query_seconds = 0.0
    if mode == "QUERY":
        log(f"Filter/grating={myfilter_grating_combo}; "
            f"minimum exposure={minexptime} s; "
            f"release year>{earliest_release_time}; "
            f"target type={TARGET_TYPE.upper()}")
        query_start = time.monotonic()
        products = query_and_select_products(query_workers)
        query_seconds = time.monotonic() - query_start
        if products is None:
            vpn_status_end = get_vpn_status()
            total_seconds = time.monotonic() - start
            write_run_metrics({
                "timestamp": run_timestamp.isoformat(timespec="seconds"),
                "run_mode": mode,
                "download_source": source,
                "filter_grating": myfilter_grating_combo,
                "target_type": TARGET_TYPE.upper(),
                "vpn_status_start": vpn_status_start,
                "vpn_status_end": vpn_status_end,
                "query_seconds": f"{query_seconds:.3f}",
                "selected_files": 0,
                "query_files_per_second": "",
                "download_seconds": "0.000",
                "downloaded_files": 0,
                "downloaded_bytes": 0,
                "downloaded_megabytes": "0.000",
                "download_megabytes_per_second": "",
                "total_elapsed_seconds": f"{total_seconds:.3f}",
                "whole_run_megabytes_per_second": "",
            }, run_timestamp)
            return
    else:
        products = read_final_product_manifest()

    download_dir = ("curls/prism_observations" if source == "MAST" and not direct_mast
                    else "downloads/prism_observations")
    download_start = time.monotonic()
    result = download_products_in_chunks(products, download_dir, source, workers)
    download_seconds = time.monotonic() - download_start
    if result is None:
        log("Download result: no files to include")
    else:
        log_count("Files available after download/restart scan", result)
    if source == "AWS" or direct_mast:
        # Restart-skipped files are present in the combined result for final
        # verification, but performance metrics describe this process only.
        downloaded_files, downloaded_bytes = local_download_totals(
            result, current_session_only=True
        )
    else:
        # MAST mode only creates curl scripts; it does not transfer FITS data.
        downloaded_files, downloaded_bytes = 0, 0
    downloaded_megabytes = downloaded_bytes / 1_000_000
    vpn_status_end = get_vpn_status()
    log(f"VPN status at end={vpn_status_end}")
    total_seconds = time.monotonic() - start
    query_file_rate = rate(len(products), query_seconds) if mode == "QUERY" else ""
    download_rate = rate(downloaded_megabytes, download_seconds)
    whole_run_rate = rate(downloaded_megabytes, total_seconds)
    metrics = {
        "timestamp": run_timestamp.isoformat(timespec="seconds"),
        "run_mode": mode,
        "download_source": source,
        "filter_grating": myfilter_grating_combo,
        "target_type": TARGET_TYPE.upper(),
        "vpn_status_start": vpn_status_start,
        "vpn_status_end": vpn_status_end,
        "query_seconds": f"{query_seconds:.3f}",
        "selected_files": len(products),
        "query_files_per_second": (f"{query_file_rate:.6f}"
                                   if query_file_rate != "" else ""),
        "download_seconds": f"{download_seconds:.3f}",
        "downloaded_files": downloaded_files,
        "downloaded_bytes": downloaded_bytes,
        "downloaded_megabytes": f"{downloaded_megabytes:.3f}",
        "download_megabytes_per_second": (f"{download_rate:.6f}"
                                          if download_rate != "" else ""),
        "total_elapsed_seconds": f"{total_seconds:.3f}",
        "whole_run_megabytes_per_second": (f"{whole_run_rate:.6f}"
                                           if whole_run_rate != "" else ""),
    }
    write_run_metrics(metrics, run_timestamp)
    log(f"Performance: query={query_seconds:.1f} s for {len(products):,} files; "
        f"download={download_seconds:.1f} s, {downloaded_megabytes:.1f} MB, "
        f"{download_rate if download_rate != '' else 0:.2f} MB/s; "
        f"whole run={total_seconds:.1f} s, "
        f"{whole_run_rate if whole_run_rate != '' else 0:.2f} MB/s")
    if source == "MAST" and not direct_mast:
        log("MAST mode created curl scripts only; downloaded byte/rate metrics are zero")
    log(f"Large NIRSpec dataset query completed in {total_seconds:.1f} s")


if __name__ == "__main__":
    try:
        arguments = parse_arguments()
        main(arguments.download_source, arguments.workers, arguments.query_workers,
             arguments.direct_mast)
    except Exception as exc:
        log(f"FAILED: {type(exc).__name__}: {exc}")
        raise

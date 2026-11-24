#!/usr/bin/env python3
"""Restart-safe retrieval of named public JWST products from AWS or MAST.

Accept filename/path lists and recursive ls inventories. Output is flat by
default; --layout date-visit restores date/VID/filename paths. AWS is tried
first with MAST fallback unless --source selects one archive explicitly.
Transfers resume .part files and complete atomically after size verification.
"""

from __future__ import annotations

from nirspec_ninja.archive import common, transfers
import argparse
import concurrent.futures
import json
import os
import random
import re
import threading
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from urllib.parse import quote

import requests

from nirspec_ninja.archive.download_metrics import DownloadMetrics

transfer_metrics = None
SOURCE_POLICY = "AUTO"


DEFAULT_WORKERS = 3
ATTEMPTS_PER_SOURCE = 5
RETRY_DELAYS_SECONDS = (5, 10, 20, 40)
CONNECT_TIMEOUT_SECONDS = 30
READ_TIMEOUT_SECONDS = 120
CHUNK_BYTES = 8 * 1024 * 1024
FINAL_GAP_PASSES = 1
DEFAULT_OUTPUT_DIR = "downloads/rawfilelist"

JWST_NAME = re.compile(r"^jw[0-9a-z_+.-]+\.fits$", re.IGNORECASE)
_print_lock = threading.Lock()


def log(message: str) -> None:
    with _print_lock:
        stamp = datetime.now().astimezone().isoformat(timespec="seconds")
        print(f"[{stamp}] {message}", flush=True)


def atomic_json(path: Path, value: dict) -> None:
    common.atomic_json(path, value)


def load_names(list_paths: list[Path], allow_missing: bool, allow_empty=False) -> list[str]:
    names: set[str] = set()
    missing = []
    for path in list_paths:
        if not path.is_file():
            missing.append(str(path))
            continue
        for number, raw_line in enumerate(path.read_text().splitlines(), 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            name = Path(line).name
            if not JWST_NAME.fullmatch(name):
                raise ValueError(f"{path}:{number}: not a JWST FITS filename: {line!r}")
            names.add(name)
        log(f"Read {path}")
    if missing and not allow_missing:
        raise FileNotFoundError(
            "Required input list(s) are missing: " + ", ".join(missing)
            + ". Restore/mount them, or intentionally use --allow-missing-lists."
        )
    for path in missing:
        log(f"WARNING: skipping missing input list: {path}")
    if not names and not allow_empty:
        raise ValueError("No filenames were found in the available input lists")
    return sorted(names)


def validate_relative_path(path: Path, source: Path, context: str) -> Path:
    """Validate a date/visit destination, discarding any old root prefix."""
    if ".." in path.parts or len(path.parts) < 3:
        raise ValueError(f"{source}: invalid destination at {context}: {path}")
    relative = Path(*path.parts[-3:])
    date, visit, name = relative.parts
    if (not re.fullmatch(r"\d{4}-\d{2}-\d{2}", date)
            or not re.fullmatch(r"VID\d{11}", visit)
            or not JWST_NAME.fullmatch(name) or name[2:13] != visit[3:]):
        raise ValueError(f"{source}: invalid date/visit path at {context}: {path}")
    datetime.strptime(date, "%Y-%m-%d")
    return relative


def parse_ls_recursive(path: Path) -> dict[str, Path]:
    """Read ls -lR headings and FITS entries; administrative paths are evidence only."""
    result: dict[str, Path] = {}
    current: Path | None = None
    for number, raw in enumerate(path.read_text(errors="replace").splitlines(), 1):
        if raw.endswith(":") and not raw.startswith((" ", "\t")):
            current = Path(raw[:-1])
            continue
        match = re.search(r"(jw[0-9A-Za-z_+.-]+\.fits)\s*$", raw, re.IGNORECASE)
        if not match:
            continue
        if current is None:
            raise ValueError(f"{path}: FITS entry before heading at line {number}")
        name = match.group(1)
        if name in result:
            raise ValueError(f"{path}: duplicate filename: {name}")
        result[name] = current / name
    if not result:
        raise ValueError(f"{path}: no JWST FITS entries in recursive inventory")
    return result


def relative_hints(list_paths: list[Path], names: list[str],
                   visit_dates: Path | None = None) -> dict[str, Path]:
    """Read trusted explicit paths and optional date/visit recovery metadata."""
    hints: dict[str, Path] = {}
    for path in list_paths:
        if not path.is_file():
            continue
        for number, raw_line in enumerate(path.read_text().splitlines(), 1):
            line = raw_line.strip()
            if not line or line.startswith("#"):
                continue
            candidate = Path(line)
            if len(candidate.parts) < 3:
                continue
            relative = validate_relative_path(candidate, path, f"line {number}")
            previous = hints.get(relative.name)
            if previous is not None and previous != relative:
                raise ValueError(f"Conflicting destinations for {relative.name}")
            hints[relative.name] = relative
    if visit_dates is not None:
        pattern = re.compile(r"(?:^|\s)(\d{4}-\d{2}-\d{2})/VID(\d{11})(?:\s|$)")
        for date, visit in pattern.findall(visit_dates.read_text()):
            # The filename list supplies exposure/detector suffixes; this file
            # supplies only the shared destination directory.
            for name in [value for value in names if value[2:13] == visit]:
                relative = validate_relative_path(Path(date) / f"VID{visit}" / name,
                                                  visit_dates, "visit metadata")
                if name in hints and hints[name] != relative:
                    raise ValueError(f"Conflicting destinations for {name}")
                hints[name] = relative
    return hints


def build_targets(names: list[str], list_paths: list[Path],
                  output_dir: Path, layout="date-visit", visit_dates=None,
                  require_paths=False) -> list[Target]:
    if layout == "flat":
        return [Target(name, output_dir / name) for name in names]
    requested = set(names)
    existing: dict[str, Path] = {}
    if output_dir.is_dir():
        for path in output_dir.rglob("*.fits"):
            if ".recovery_staging" in path.relative_to(output_dir).parts:
                continue
            if path.name not in requested:
                continue
            if path.name in existing:
                raise RuntimeError(f"Duplicate recovered filename: {path.name}: "
                                   f"{existing[path.name]} and {path}")
            existing[path.name] = path
    hints = relative_hints(list_paths, names, visit_dates)
    if require_paths:
        missing = sorted(set(names) - hints.keys())
        if missing:
            raise ValueError(f"No trusted date/visit path supplied for {missing[0]}")
    targets = []
    staging = output_dir / ".recovery_staging"
    for name in names:
        if name in existing and name in hints and existing[name] != output_dir / hints[name]:
            raise ValueError(f"Existing file conflicts with supplied path: {existing[name]}")
        if name in existing:
            targets.append(Target(name, existing[name]))
        elif name in hints:
            targets.append(Target(name, output_dir / hints[name]))
        else:
            # DATE-OBS is read after transfer; no observation date is guessed.
            targets.append(Target(name, staging / name, dynamic_path=True))
    log(f"Recovered-tree inventory: {len(existing)}/{len(names)} requested files "
        f"present; {len(names) - len(existing)} absent")
    log(f"Destination paths: {len(hints)} supplied by recovery metadata; "
        f"{sum(target.dynamic_path for target in targets)} require DATE-OBS")
    return targets


def urls_for(name: str) -> tuple[tuple[str, str], ...]:
    """Return deterministic public archive URLs, AWS first."""
    dataset = name.split("_", 1)[0]
    if not re.fullmatch(r"jw\d{11}", dataset, re.IGNORECASE):
        raise ValueError(f"Cannot derive JWST dataset directory from {name!r}")
    program = dataset[:7]
    key = f"jwst/public/{program}/{dataset}/{name}"
    aws = f"https://stpubdata.s3.amazonaws.com/{key}"
    mast = ("https://mast.stsci.edu/api/v0.1/Download/file?uri="
            + quote(f"mast:JWST/product/{name}", safe=":/"))
    sources = (("AWS", aws), ("MAST", mast))
    return sources if SOURCE_POLICY == "AUTO" else tuple(
        pair for pair in sources if pair[0] == SOURCE_POLICY
    )


@dataclass(frozen=True)
class RemoteInfo:
    source: str
    url: str
    size: int
    etag: str | None


@dataclass
class Target:
    name: str
    destination: Path
    dynamic_path: bool = False


def retry_delay(attempt: int) -> float:
    return RETRY_DELAYS_SECONDS[attempt - 1] + random.uniform(0, 1)


def request_with_retries(session: requests.Session, method: str, url: str,
                         source: str, name: str, attempts: int = ATTEMPTS_PER_SOURCE,
                         **kwargs) -> requests.Response:
    def operation():
        response = session.request(method, url,
            timeout=(CONNECT_TIMEOUT_SECONDS, READ_TIMEOUT_SECONDS),
            allow_redirects=True, **kwargs)
        if response.status_code == 404:
            response.close()
            raise FileNotFoundError(f"not found at {source}")
        response.raise_for_status()
        return response
    def record(error, started, attempt):
        if method == "HEAD" and transfer_metrics is not None:
            transfer_metrics.record(source, name, "HEAD", started, attempt, error=error)
    return transfers.retry_transfer(operation, attempts=attempts,
        exceptions=(requests.RequestException, FileNotFoundError),
        clock=time.monotonic, sleep=time.sleep, delay=retry_delay,
        success=lambda result, started, attempt: record(None, started, attempt),
        failure=record,
        retry_message=lambda error, attempt, seconds: log(
            f"{name}: {source} {method} attempt {attempt}/{attempts} failed ({error}); retrying in {seconds:.1f}s"),
        exhausted=lambda error: RuntimeError(
            f"{source} {method} failed after {attempts} attempts: {error}"))


def inspect_remote(session: requests.Session, source: str, url: str,
                   name: str) -> RemoteInfo:
    response = request_with_retries(session, "HEAD", url, source, name)
    try:
        size_text = response.headers.get("Content-Length")
        if not size_text or int(size_text) <= 0:
            raise RuntimeError(f"{source} returned no valid Content-Length")
        return RemoteInfo(source, response.url, int(size_text),
                          response.headers.get("ETag"))
    finally:
        response.close()


def completion_valid(destination: Path, metadata_path: Path) -> bool:
    if not destination.is_file() or not metadata_path.is_file():
        return False
    try:
        metadata = json.loads(metadata_path.read_text())
        return (metadata.get("filename") == destination.name
                and destination.stat().st_size == int(metadata["size"])
                and int(metadata["size"]) > 0)
    except (OSError, ValueError, KeyError, TypeError, json.JSONDecodeError):
        return False


def transfer_once(session: requests.Session, remote: RemoteInfo,
                  destination: Path, metadata_path: Path, received=None) -> None:
    part = destination.with_suffix(destination.suffix + ".part")
    current = part.stat().st_size if part.exists() else 0
    if current > remote.size:
        log(f"{destination.name}: discarding oversized partial file")
        part.unlink()
        current = 0
    if current == remote.size:
        os.replace(part, destination)
    else:
        headers = {"Range": f"bytes={current}-"} if current else {}
        response = request_with_retries(
            session, "GET", remote.url, remote.source, destination.name,
            attempts=1, headers=headers, stream=True,
        )
        try:
            if current and response.status_code != 206:
                # The server ignored Range. Restart this file without appending.
                current = 0
                part.unlink(missing_ok=True)
            mode = "ab" if current else "wb"
            with part.open(mode) as stream:
                for chunk in response.iter_content(CHUNK_BYTES):
                    if chunk:
                        stream.write(chunk)
                        if received is not None:
                            received[0] += len(chunk)
                stream.flush()
                os.fsync(stream.fileno())
        finally:
            response.close()
        actual = part.stat().st_size
        if actual != remote.size:
            raise IOError(f"incomplete transfer: local={actual}, remote={remote.size}")
        os.replace(part, destination)
    atomic_json(metadata_path, {
        "filename": destination.name,
        "size": remote.size,
        "source": remote.source,
        "url": remote.url,
        "etag": remote.etag,
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    })


def transfer(session: requests.Session, remote: RemoteInfo, destination: Path,
             metadata_path: Path) -> None:
    """Transfer with five attempts; each retry resumes the growing .part file."""
    received = [0]
    def operation():
        received[0] = 0
        transfer_once(session, remote, destination, metadata_path, received)
    def record(error, started, attempt):
        if transfer_metrics is not None:
            transfer_metrics.record(remote.source, destination.name, "download",
                started, attempt, received[0], error=error)
    transfers.retry_transfer(operation, attempts=ATTEMPTS_PER_SOURCE,
        exceptions=(requests.RequestException, OSError, RuntimeError),
        clock=time.monotonic, sleep=time.sleep, delay=retry_delay,
        success=lambda result, started, attempt: record(None, started, attempt),
        failure=record,
        retry_message=lambda error, attempt, seconds: log(
            f"{destination.name}: {remote.source} transfer attempt {attempt}/{ATTEMPTS_PER_SOURCE} failed ({error}); resuming in {seconds:.1f}s"),
        exhausted=lambda error: RuntimeError(
            f"{remote.source} transfer failed after {ATTEMPTS_PER_SOURCE} attempts: {error}"))


def metadata_path(destination: Path) -> Path:
    return destination.with_name(destination.name + ".download.json")


def archive_match(destination: Path, remote: RemoteInfo) -> bool:
    return destination.is_file() and destination.stat().st_size == remote.size


def record_completion(destination: Path, remote: RemoteInfo) -> None:
    atomic_json(metadata_path(destination), {
        "filename": destination.name,
        "size": remote.size,
        "source": remote.source,
        "url": remote.url,
        "etag": remote.etag,
        "completed_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    })


def fits_observation_date(path: Path) -> str:
    """Read DATE-OBS without loading the large science arrays."""
    from astropy.io import fits
    value = str(fits.getheader(path, 0).get("DATE-OBS", ""))[:10]
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", value):
        raise RuntimeError(f"No valid DATE-OBS in downloaded FITS file {path}")
    return value


def finalize_dynamic_target(target: Target, root: Path) -> None:
    date = fits_observation_date(target.destination)
    visit = "VID" + target.name[2:13]
    final = root / date / visit / target.name
    final.parent.mkdir(parents=True, exist_ok=True)
    os.replace(target.destination, final)
    old_metadata = metadata_path(target.destination)
    if old_metadata.exists():
        os.replace(old_metadata, metadata_path(final))
    target.destination = final
    target.dynamic_path = False


def download_one(target: Target, output_dir: Path) -> tuple[str, str]:
    name = target.name
    destination = target.destination
    metadata = metadata_path(destination)
    errors = []
    with requests.Session() as session:
        session.headers["User-Agent"] = "jwst-filelist-recovery/1.0"
        for source, url in urls_for(name):
            try:
                remote = inspect_remote(session, source, url, name)
                if archive_match(destination, remote):
                    record_completion(destination, remote)
                    if transfer_metrics is not None:
                        transfer_metrics.record(source, name, "download", time.monotonic(),
                                                status="skipped")
                    if target.dynamic_path:
                        finalize_dynamic_target(target, output_dir)
                    return name, f"verified existing file against {source}"
                destination.parent.mkdir(parents=True, exist_ok=True)
                transfer(session, remote, destination, metadata)
                if target.dynamic_path:
                    finalize_dynamic_target(target, output_dir)
                return name, f"downloaded from {source}"
            except Exception as exc:
                errors.append(f"{source}: {type(exc).__name__}: {exc}")
                log(f"{name}: {source} exhausted; trying next source")
    raise RuntimeError("; ".join(errors))


def scan_gaps(targets: list[Target]) -> list[Target]:
    return [target for target in targets if target.dynamic_path or not completion_valid(
        target.destination, metadata_path(target.destination))]


def run_batch(targets: list[Target], output_dir: Path,
              workers: int) -> dict[str, str]:
    global transfer_metrics
    if transfer_metrics is None:
        transfer_metrics = DownloadMetrics(
            output_dir / "transfer_metrics", "filelist", workers,
            f"connect={CONNECT_TIMEOUT_SECONDS};read={READ_TIMEOUT_SECONDS}",
        )
    failures: dict[str, str] = {}
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(download_one, target, output_dir): target.name
                   for target in targets}
        completed = 0
        for future in concurrent.futures.as_completed(futures):
            name = futures[future]
            completed += 1
            try:
                _, result = future.result()
                log(f"[{completed}/{len(targets)}] {name}: {result}")
            except Exception as exc:
                failures[name] = f"{type(exc).__name__}: {exc}"
                log(f"[{completed}/{len(targets)}] {name}: FAILED: {exc}")
    return failures


def write_reports(targets: list[Target], output_dir: Path,
                  failures: dict[str, str]) -> None:
    gaps = scan_gaps(targets)
    missing_names = [target.name for target in gaps]
    atomic_json(output_dir / "recovery_report.json", {
        "requested": len(targets),
        "complete": len(targets) - len(gaps),
        "missing_or_invalid": missing_names,
        "failures": {name: failures.get(name, "verification failed")
                     for name in missing_names},
        "checked_at": datetime.now().astimezone().isoformat(timespec="seconds"),
    })
    (output_dir / "missing_files.txt").write_text(
        "".join(f"{name}\n" for name in missing_names)
    )


def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("lists", nargs="*", type=Path,
                        help="one or more filename or date/visit path lists")
    parser.add_argument("--inventory", action="append", type=Path, default=[],
                        help="recursive ls -lR inventory; may be repeated")
    parser.add_argument("--layout", choices=("flat", "date-visit"), default="flat")
    parser.add_argument("--visit-dates", type=Path,
                        help="optional file containing date/VID visit paths")
    parser.add_argument("--require-paths", action="store_true",
                        help="require explicit trusted date/visit paths for every filename")
    parser.add_argument("--dry-run", action="store_true",
                        help="validate inputs and show destinations without network calls or writes")
    parser.add_argument("-o", "--output-dir", type=Path,
                        default=Path(DEFAULT_OUTPUT_DIR))
    parser.add_argument("-j", "--workers", type=int, default=DEFAULT_WORKERS,
                        help=f"parallel transfers (default: {DEFAULT_WORKERS})")
    parser.add_argument("--allow-missing-lists", action="store_true",
                        help="continue when an input list is unavailable")
    parser.add_argument("--source", choices=("AUTO", "AWS", "MAST"), default="AUTO",
                        help="AUTO tries AWS then MAST; choose one for source comparisons")
    args = parser.parse_args(argv)
    if not args.lists and not args.inventory:
        parser.error("supply at least one filename list or --inventory")
    if args.layout == "flat" and (args.require_paths or args.visit_dates):
        parser.error("--require-paths and --visit-dates require --layout date-visit")
    return args


def main(argv=None) -> int:
    global SOURCE_POLICY, transfer_metrics
    args = parse_args(argv)
    SOURCE_POLICY = args.source
    transfer_metrics = None
    if args.workers < 1:
        raise SystemExit("--workers must be at least 1")
    paths = args.lists
    names = set(load_names(paths, args.allow_missing_lists, bool(args.inventory))) if paths else set()
    for inventory in args.inventory:
        if not inventory.is_file() and args.allow_missing_lists:
            log(f"WARNING: skipping missing inventory: {inventory}")
            continue
        recovered = parse_ls_recursive(inventory)
        log(f"Inventory {inventory}: {len(recovered)} unique filenames")
        names.update(recovered)
    if not names:
        raise ValueError("No filenames found")
    names = sorted(names)
    for name in names:
        urls_for(name)  # Validate dataset IDs before any network calls or writes.
    targets = build_targets(names, paths, args.output_dir, args.layout,
                            args.visit_dates, args.require_paths)
    log(f"Recovery set: {len(names)} unique files; output={args.output_dir}; "
        f"workers={args.workers}")
    if args.dry_run:
        for target in targets:
            suffix = " (date will come from FITS DATE-OBS)" if target.dynamic_path else ""
            print(f"{target.name} -> {target.destination}{suffix}")
        return 0
    args.output_dir.mkdir(parents=True, exist_ok=True)
    failures = run_batch(targets, args.output_dir, args.workers)
    for pass_number in range(1, FINAL_GAP_PASSES + 1):
        gaps = scan_gaps(targets)
        if not gaps:
            break
        log(f"Final verification found {len(gaps)} gap(s); fill pass "
            f"{pass_number}/{FINAL_GAP_PASSES}")
        failures.update(run_batch(gaps, args.output_dir, args.workers))
    write_reports(targets, args.output_dir, failures)
    gaps = scan_gaps(targets)
    if gaps:
        log(f"INCOMPLETE: {len(gaps)} file(s) remain; see missing_files.txt "
            "and recovery_report.json")
        return 1
    log(f"COMPLETE: verified all {len(names)} requested files")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

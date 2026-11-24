"""Common, durable attempt metrics for AWS and MAST archive operations."""

import csv
from datetime import datetime, timezone
from pathlib import Path
import threading
import time


def is_timeout(error):
    """Recognize typed timeouts, including chained transport exceptions."""
    seen = set()
    while error is not None and id(error) not in seen:
        seen.add(id(error))
        if isinstance(error, TimeoutError) or "timeout" in type(error).__name__.lower():
            return True
        error = error.__cause__ or error.__context__
    return False


class DownloadMetrics:
    """Write one row per explicit attempt; keep sources and HEAD separate."""

    def __init__(self, directory, workflow, workers=1, timeout_seconds=None):
        stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
        self.directory = Path(directory) / stamp
        self.workflow = workflow
        self.workers = workers
        self.timeout_seconds = timeout_seconds
        self.rows = []
        self.started = time.monotonic()
        self.lock = threading.Lock()

    def record(self, source, filename, operation, started, attempt=1,
               byte_count=0, error=None, status=None):
        row = dict(workflow=self.workflow, source=source, filename=filename,
                   operation=operation, attempt=attempt, workers=self.workers,
                   timeout_seconds=self.timeout_seconds,
                   start_seconds=round(started - self.started, 6),
                   seconds=round(time.monotonic() - started, 6),
                   bytes=byte_count, status=status or ("failed" if error else "ok"),
                   timeout=int(is_timeout(error)),
                   error="" if error is None else f"{type(error).__name__}: {error}")
        with self.lock:
            self.directory.mkdir(parents=True, exist_ok=True)
            self.rows.append(row)
            self._write("attempts.csv", self.rows)
            summaries = []
            for source_name, op in sorted({(r["source"], r["operation"]) for r in self.rows}):
                rows = [r for r in self.rows if (r["source"], r["operation"]) == (source_name, op)]
                attempts = [r for r in rows if r["status"] != "skipped"]
                successes = [r for r in attempts if r["status"] == "ok"]
                seconds = sum(r["seconds"] for r in attempts)
                size = sum(r["bytes"] for r in attempts)
                timeouts = sum(r["timeout"] for r in attempts)
                wall = (max(r["start_seconds"] + r["seconds"] for r in attempts)
                        - min(r["start_seconds"] for r in attempts)) if attempts else 0
                summaries.append(dict(
                    workflow=self.workflow, source=source_name, operation=op,
                    workers=self.workers, timeout_seconds=self.timeout_seconds,
                    attempts=len(attempts), successes=len(successes),
                    failures=len(attempts) - len(successes),
                    retries=sum(r["attempt"] > 1 for r in attempts),
                    skipped=sum(r["status"] == "skipped" for r in rows),
                    timeouts=timeouts,
                    timeout_fraction=timeouts / len(attempts) if attempts else "",
                    bytes=size, attempt_seconds=round(seconds, 6),
                    wall_seconds=round(wall, 6),
                    mb_per_wall_second=size / 1e6 / wall if wall else "",
                    mb_per_attempt_second=size / 1e6 / seconds if seconds else "",
                ))
            self._write("summary.csv", summaries)

    def _write(self, name, rows):
        path = self.directory / name
        temporary = path.with_suffix(".tmp")
        with temporary.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader()
            writer.writerows(rows)
        temporary.replace(path)

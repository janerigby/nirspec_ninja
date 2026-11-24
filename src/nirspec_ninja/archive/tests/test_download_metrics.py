"""Offline checks: python -m unittest discover -s tests."""

import csv
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import tempfile
import time
import unittest
from unittest.mock import patch

from nirspec_ninja.archive.download_metrics import DownloadMetrics, is_timeout


class MetricsTests(unittest.TestCase):
    def test_exact_rates_include_retry_gaps_and_separate_sources_and_head(self):
        with tempfile.TemporaryDirectory() as directory, patch("nirspec_ninja.archive.download_metrics.time.monotonic", return_value=0) as clock:
            metrics = DownloadMetrics(directory, "test", workers=2)
            clock.return_value = 2
            metrics.record("AWS", "a", "download", 0, byte_count=1_000_000, error=TimeoutError())
            clock.return_value = 5
            metrics.record("AWS", "b", "download", 1, byte_count=2_000_000)
            clock.return_value = 11
            metrics.record("AWS", "a", "download", 7, attempt=2, byte_count=3_000_000)
            clock.return_value = 20
            metrics.record("AWS", "cached", "download", 20, byte_count=999_000_000, status="skipped")
            clock.return_value = 100
            metrics.record("MAST", "a", "HEAD", 0, error=TimeoutError())
            clock.return_value = 4
            metrics.record("MAST", "a", "download", 1, byte_count=6_000_000)
            with (metrics.directory / "summary.csv").open() as stream:
                rows = {(row["source"], row["operation"]): row for row in csv.DictReader(stream)}
            aws = rows[("AWS", "download")]
            self.assertEqual(float(aws["attempt_seconds"]), 10)
            self.assertEqual(float(aws["wall_seconds"]), 11)
            self.assertEqual(float(aws["mb_per_attempt_second"]), 0.6)
            self.assertAlmostEqual(float(aws["mb_per_wall_second"]), 6 / 11)
            self.assertAlmostEqual(float(aws["timeout_fraction"]), 1 / 3)
            self.assertEqual(aws["bytes"], "6000000")
            mast = rows[("MAST", "download")]
            self.assertEqual(float(mast["mb_per_wall_second"]), 2)
            self.assertEqual(float(mast["timeout_fraction"]), 0)
            self.assertEqual(float(rows[("MAST", "HEAD")]["timeout_fraction"]), 1)

    def test_zero_duration_and_skip_only_have_no_invented_rates(self):
        with tempfile.TemporaryDirectory() as directory, patch("nirspec_ninja.archive.download_metrics.time.monotonic", return_value=0):
            metrics = DownloadMetrics(directory, "test")
            metrics.record("AWS", "cached", "download", 0, status="skipped")
            metrics.record("MAST", "empty", "download", 0)
            with (metrics.directory / "summary.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            self.assertEqual(rows[0]["timeout_fraction"], "")
            for row in rows:
                self.assertEqual(row["mb_per_wall_second"], "")
                self.assertEqual(row["mb_per_attempt_second"], "")

    def test_wrapped_timeout(self):
        error = RuntimeError("transfer failed")
        error.__cause__ = TimeoutError("read expired")
        self.assertTrue(is_timeout(error))
        self.assertFalse(is_timeout(FileNotFoundError("404")))

    def test_concurrent_failures_retries_and_skips(self):
        with tempfile.TemporaryDirectory() as directory:
            metrics = DownloadMetrics(directory, "test", workers=2)
            def record(index):
                metrics.record("AWS", f"file{index}", "download",
                               time.monotonic() - 1, byte_count=100)
            with ThreadPoolExecutor(max_workers=2) as pool:
                list(pool.map(record, range(10)))
            metrics.record("AWS", "bad", "download", time.monotonic() - 1,
                           attempt=2, error=TimeoutError())
            metrics.record("AWS", "cached", "download", time.monotonic(),
                           byte_count=999, status="skipped")
            metrics.record("MAST", "bad", "HEAD", time.monotonic(),
                           error=FileNotFoundError())
            with (metrics.directory / "summary.csv").open() as stream:
                rows = list(csv.DictReader(stream))
            aws = next(row for row in rows if row["source"] == "AWS")
            self.assertEqual(int(aws["attempts"]), 11)
            self.assertEqual(int(aws["successes"]), 10)
            self.assertEqual(int(aws["failures"]), 1)
            self.assertEqual(int(aws["skipped"]), 1)
            self.assertEqual(int(aws["retries"]), 1)
            self.assertEqual(int(aws["bytes"]), 1000)
            self.assertAlmostEqual(float(aws["timeout_fraction"]), 1 / 11)
            self.assertGreater(float(aws["mb_per_wall_second"]), 0)
            self.assertFalse(list(Path(metrics.directory).glob("*.tmp")))


if __name__ == "__main__":
    unittest.main()

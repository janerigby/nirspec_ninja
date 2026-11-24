"""Network-boundary doubles exercise real retry and resume loops."""

from pathlib import Path
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import requests
from astropy.table import Table

from nirspec_ninja.archive.download_metrics import DownloadMetrics
from nirspec_ninja.archive import MAST_download_rawfilelist as raw
from nirspec_ninja.archive import MAST_query_large_nirspec_datasets as query
from nirspec_ninja.archive.tests.helpers import NAME, Response, read_csv


class RecoveryTransferTests(unittest.TestCase):
    def test_head_timeout_then_success_records_retry(self):
        with tempfile.TemporaryDirectory() as directory:
            metrics = DownloadMetrics(directory, "test")
            session = MagicMock()
            response = Response()
            session.request.side_effect = [requests.ReadTimeout("slow"), response]
            with patch.object(raw, "transfer_metrics", metrics), patch.object(raw.time, "sleep") as sleep:
                remote = raw.inspect_remote(session, "AWS", "https://example.invalid", NAME)
            self.assertEqual(remote.size, 6)
            self.assertTrue(response.closed)
            sleep.assert_called_once()
            summary = read_csv(metrics.directory / "summary.csv")[0]
            self.assertEqual((summary["attempts"], summary["timeouts"], summary["retries"]), ("2", "1", "1"))
            self.assertEqual(session.request.call_args.kwargs["timeout"], (30, 120))

    def test_head_exhaustion_and_invalid_archive_sizes(self):
        session = MagicMock()
        session.request.side_effect = requests.ConnectTimeout("offline")
        with patch.object(raw.time, "sleep") as sleep, patch.object(raw, "transfer_metrics", None):
            with self.assertRaisesRegex(RuntimeError, "after 3 attempts") as error:
                raw.request_with_retries(session, "HEAD", "https://example.invalid", "MAST", NAME, attempts=3)
        self.assertIsInstance(error.exception.__cause__, requests.ConnectTimeout)
        self.assertEqual(session.request.call_count, 3)
        self.assertEqual(sleep.call_count, 2)
        for headers in [{}, {"Content-Length": "0"}, {"Content-Length": "-3"}]:
            response = Response(headers=headers)
            with self.subTest(headers=headers), patch.object(raw, "request_with_retries", return_value=response):
                with self.assertRaisesRegex(RuntimeError, "Content-Length"):
                    raw.inspect_remote(session, "MAST", "https://example.invalid", NAME)
                self.assertTrue(response.closed)

    def test_partial_timeout_then_resume_counts_each_attempt(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = root / NAME
            metrics = DownloadMetrics(root / "metrics", "test")
            first = Response(chunks=[b"abc", requests.ReadTimeout("interrupted")])
            second = Response(status=206, chunks=[b"def"])
            session = MagicMock()
            session.request.side_effect = [first, second]
            with patch.object(raw, "transfer_metrics", metrics), patch.object(raw.time, "sleep") as sleep:
                raw.transfer(session, raw.RemoteInfo("AWS", "https://example.invalid", 6, None),
                             target, raw.metadata_path(target))
            self.assertEqual(target.read_bytes(), b"abcdef")
            self.assertTrue(raw.completion_valid(target, raw.metadata_path(target)))
            self.assertTrue(first.closed and second.closed)
            self.assertEqual(session.request.call_args.kwargs["headers"], {"Range": "bytes=3-"})
            rows = read_csv(metrics.directory / "attempts.csv")
            self.assertEqual([row["bytes"] for row in rows], ["3", "3"])
            self.assertEqual([row["status"] for row in rows], ["failed", "ok"])
            self.assertEqual(rows[0]["timeout"], "1")
            sleep.assert_called_once()

    def test_ignored_range_restarts_instead_of_appending(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / NAME
            target.with_suffix(".fits.part").write_bytes(b"old")
            session = MagicMock()
            session.request.return_value = Response(b"abcdef", status=200)
            count = [0]
            raw.transfer_once(session, raw.RemoteInfo("MAST", "https://example.invalid", 6, None),
                              target, raw.metadata_path(target), count)
            self.assertEqual(target.read_bytes(), b"abcdef")
            self.assertEqual(count, [6])

    def test_short_response_preserves_old_final_and_partial(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / NAME
            target.write_bytes(b"old")
            session = MagicMock()
            session.request.return_value = Response(b"ab")
            with self.assertRaisesRegex(OSError, "incomplete transfer"):
                raw.transfer_once(session, raw.RemoteInfo("AWS", "https://example.invalid", 6, None),
                                  target, raw.metadata_path(target))
            self.assertEqual(target.read_bytes(), b"old")
            self.assertEqual(target.with_suffix(".fits.part").read_bytes(), b"ab")
            self.assertFalse(raw.metadata_path(target).exists())

    def test_transfer_exhaustion_records_all_timeouts(self):
        with tempfile.TemporaryDirectory() as directory:
            target = Path(directory) / NAME
            metrics = DownloadMetrics(Path(directory) / "metrics", "test")
            session = MagicMock()
            session.request.side_effect = requests.ReadTimeout("slow")
            with patch.object(raw, "transfer_metrics", metrics), patch.object(raw, "ATTEMPTS_PER_SOURCE", 3), \
                    patch.object(raw.time, "sleep") as sleep:
                with self.assertRaisesRegex(RuntimeError, "after 3 attempts"):
                    raw.transfer(session, raw.RemoteInfo("MAST", "https://example.invalid", 6, None),
                                 target, raw.metadata_path(target))
            summary = read_csv(metrics.directory / "summary.csv")[0]
            self.assertEqual((summary["attempts"], summary["failures"], summary["timeouts"]), ("3", "3", "3"))
            self.assertEqual(summary["successes"], "0")
            self.assertEqual(sleep.call_count, 2)
            self.assertFalse(target.exists())


class LargeQueryTransferTests(unittest.TestCase):
    def product(self):
        return Table({"productFilename": [NAME], "dataURI": ["mast:JWST/product/" + NAME], "size": [6]})

    def test_timeout_then_success_uses_new_clients(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            metrics = DownloadMetrics(root / "metrics", "test")
            first, second = MagicMock(), MagicMock()
            first.download_products.side_effect = requests.ReadTimeout("slow")
            def download(product, download_dir, **kwargs):
                target = Path(download_dir) / "nested" / NAME
                target.parent.mkdir()
                target.write_bytes(b"abcdef")
                return Table({"Local Path": [str(target)]})
            second.download_products.side_effect = download
            with patch.object(query, "ObservationsClass", side_effect=[first, second]), \
                    patch.object(query, "transfer_metrics", metrics), patch.object(query.time, "sleep") as sleep:
                result = query.download_one_aws_product(self.product(), root)
            self.assertEqual(result.read_bytes(), b"abcdef")
            self.assertTrue(second.download_products.call_args.kwargs["cloud_only"])
            sleep.assert_called_once_with(5)
            rows = read_csv(metrics.directory / "attempts.csv")
            self.assertEqual([row["timeout"] for row in rows], ["1", "0"])
            self.assertEqual([row["bytes"] for row in rows], ["0", "6"])

    def test_failed_transfer_does_not_replace_existing_file(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / NAME).write_bytes(b"old")
            client = MagicMock()
            def download(product, download_dir, **kwargs):
                target = Path(download_dir) / NAME
                target.write_bytes(b"xx")
                return Table({"Local Path": [str(target)]})
            client.download_products.side_effect = download
            with patch.object(query, "ObservationsClass", return_value=client), \
                    patch.object(query, "AWS_DOWNLOAD_ATTEMPTS", 2), \
                    patch.object(query, "transfer_metrics", None), patch.object(query.time, "sleep"):
                with self.assertRaises(RuntimeError):
                    query.download_one_aws_product(self.product(), root, source="MAST")
            self.assertEqual((root / NAME).read_bytes(), b"old")
            self.assertEqual(client.download_products.call_count, 2)
            self.assertFalse(client.download_products.call_args.kwargs["cloud_only"])
            client.enable_cloud_dataset.assert_not_called()

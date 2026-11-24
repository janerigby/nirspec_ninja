"""Offline regression tests for consolidated inventory and transfer behavior."""

import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch, MagicMock

from nirspec_ninja.archive import MAST_download_rawfilelist as raw


NAME = "jw01473016001_02101_00001_nrs1_uncal.fits"
RELATIVE = Path("2022-04-23/VID01473016001") / NAME


class RecoveryTests(unittest.TestCase):
    def test_inventory_administrative_paths_do_not_override_explicit(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            paths = root / "paths.txt"
            paths.write_text(f"RAW_example/{RELATIVE}\n")
            inventory = root / "inventory.txt"
            inventory.write_text(f"RAW_example/Bad_moving:\n-rw-r--r-- 1 me staff 8 Oct 1 2026 {NAME}\n")
            names = raw.parse_ls_recursive(inventory)
            self.assertEqual(set(names), {NAME})
            targets = raw.build_targets(sorted(names), [paths], root / "output",
                                        require_paths=True)
            self.assertEqual(targets[0].destination, root / "output" / RELATIVE)
            self.assertFalse(targets[0].dynamic_path)

    def test_conflicting_and_unsafe_paths_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            paths = Path(directory) / "paths.txt"
            paths.write_text(f"{RELATIVE}\n2022-04-24/VID01473016001/{NAME}\n")
            with self.assertRaises(ValueError):
                raw.relative_hints([paths], [NAME])
            with self.assertRaises(ValueError):
                raw.validate_relative_path(Path("../") / RELATIVE, paths, "test")
            with self.assertRaises(ValueError):
                raw.validate_relative_path(Path("2022-04-23/VID00000000000") / NAME,
                                           paths, "test")

    def test_dry_run_has_no_network_or_output_writes(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            names = root / "names.txt"
            names.write_text(NAME + "\n")
            with patch.object(raw, "run_batch") as batch, contextlib.redirect_stdout(io.StringIO()):
                result = raw.main([str(names), "--dry-run", "-o", str(root / "output")])
            self.assertEqual(result, 0)
            batch.assert_not_called()
            self.assertFalse((root / "output").exists())

    def test_partial_transfer_resumes_and_counts_only_new_bytes(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / NAME
            path.with_suffix(".fits.part").write_bytes(b"abc")
            response = MagicMock()
            response.status_code = 206
            response.iter_content.return_value = [b"def"]
            session = MagicMock()
            session.request.return_value = response
            received = [0]
            raw.transfer_once(session, raw.RemoteInfo("AWS", "https://example.invalid", 6, None),
                              path, raw.metadata_path(path), received)
            self.assertEqual(path.read_bytes(), b"abcdef")
            self.assertEqual(received, [3])
            self.assertEqual(session.request.call_args.kwargs["headers"], {"Range": "bytes=3-"})
            self.assertTrue(raw.completion_valid(path, raw.metadata_path(path)))

    def test_aws_failure_falls_back_to_mast(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            target = raw.Target(NAME, root / NAME)
            remote = raw.RemoteInfo("MAST", "https://example.invalid", 6, None)
            with patch.object(raw, "SOURCE_POLICY", "AUTO"), patch.object(raw, "transfer_metrics", None), \
                    patch.object(raw, "inspect_remote", side_effect=[TimeoutError(), remote]) as inspect, \
                    patch.object(raw, "transfer") as transfer:
                name, status = raw.download_one(target, root)
            self.assertEqual(name, NAME)
            self.assertEqual(status, "downloaded from MAST")
            self.assertEqual([call.args[1] for call in inspect.call_args_list], ["AWS", "MAST"])
            transfer.assert_called_once()

    def test_unplaced_staging_file_remains_a_gap(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            stage = root / ".recovery_staging" / NAME
            stage.parent.mkdir()
            stage.write_bytes(b"abcdef")
            remote = raw.RemoteInfo("AWS", "https://example.invalid", 6, None)
            raw.record_completion(stage, remote)
            targets = raw.build_targets([NAME], [], root)
            self.assertTrue(targets[0].dynamic_path)
            self.assertEqual(raw.scan_gaps(targets), targets)
            with patch.object(raw, "transfer_metrics", None), patch.object(raw, "inspect_remote", return_value=remote), \
                    patch.object(raw, "fits_observation_date", return_value="2022-04-23"):
                raw.download_one(targets[0], root)
            self.assertEqual(targets[0].destination, root / RELATIVE)
            self.assertEqual(raw.scan_gaps(targets), [])

    def test_source_policy(self):
        with patch.object(raw, "SOURCE_POLICY", "MAST"):
            self.assertEqual([source for source, _ in raw.urls_for(NAME)], ["MAST"])


if __name__ == "__main__":
    unittest.main()

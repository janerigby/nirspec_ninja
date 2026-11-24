"""Public command names and invalid arguments remain usable and predictable."""

import contextlib
import io
from pathlib import Path
import subprocess
import sys
import unittest
from unittest.mock import patch

from nirspec_ninja.archive import MAST_download_PID as pid
from nirspec_ninja.archive import MAST_download_rawfilelist as raw
from nirspec_ninja.archive import MAST_query_large_nirspec_datasets as query


class CLITests(unittest.TestCase):
    def test_multiple_pids_preserve_single_pid_source_syntax(self):
        for arguments, expected, source in [
                (["04125", "AWS"], [4125], "AWS"),
                (["04125", "03843", "AWS"], [4125, 3843], "AWS"),
                (["4125", "3843", "4125", "--source", "MAST"], [4125, 3843], "MAST"),
                (["4125", "--instrument", "miri", "3843", "MAST"], [4125, 3843], "MAST")]:
            with self.subTest(arguments=arguments):
                args = pid.parse_args(arguments)
                self.assertEqual(args.pids, expected)
                self.assertEqual(args.download_source, source)

    def test_pid_list_validation_and_default_confirmation(self):
        for arguments in [["0"], ["-3"], ["not-a-pid"], ["4125", "AWS", "--source", "MAST"]]:
            with self.subTest(arguments=arguments), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as exit:
                    pid.parse_args(arguments)
                self.assertEqual(exit.exception.code, 2)
        with patch("builtins.input", return_value="yes"):
            self.assertEqual(pid.choose_pids([]), [4125, 3843])
        with patch("builtins.input", return_value="no"), patch.object(pid.Observations, "query_criteria") as archive:
            self.assertEqual(pid.main(["--source", "MAST"]), 0)
            archive.assert_not_called()

    def test_all_three_commands_import_and_show_help(self):
        root = Path(__file__).resolve().parent.parent
        for name in ["MAST_query_large_nirspec_datasets.py", "MAST_download_PID.py", "MAST_download_rawfilelist.py"]:
            with self.subTest(command=name):
                result = subprocess.run([sys.executable, "-B", str(root / name), "--help"],
                                        cwd=root, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn("usage:", result.stdout)

    def test_invalid_concurrency_and_missing_inputs_are_rejected(self):
        for arguments in [["--workers", "0"], ["--query-workers", "-1"]]:
            with self.subTest(arguments=arguments), contextlib.redirect_stderr(io.StringIO()):
                with self.assertRaises(SystemExit) as exit:
                    query.parse_arguments(arguments)
                self.assertEqual(exit.exception.code, 2)
        with contextlib.redirect_stderr(io.StringIO()), self.assertRaises(SystemExit):
            raw.parse_args([])
        with self.assertRaises(ValueError):
            pid.choose_pids([0])

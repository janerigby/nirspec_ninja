"""Run real workflow orchestration, files, manifests, and reports offline."""

import contextlib
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import requests
from astropy.io import fits
from astropy.table import Table

from nirspec_ninja.archive import MAST_download_PID as pid
from nirspec_ninja.archive import MAST_download_rawfilelist as raw
from nirspec_ninja.archive import MAST_query_large_nirspec_datasets as query
from nirspec_ninja.archive.tests.helpers import NAME, Response, read_csv, working_directory


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.stack = contextlib.ExitStack()
        self.addCleanup(self.stack.close)
        self.root = Path(self.stack.enter_context(tempfile.TemporaryDirectory()))
        self.stack.enter_context(working_directory(self.root))
        self.stack.enter_context(patch.object(query, "get_vpn_status", return_value="unavailable"))
        for key, value in dict(minexptime=0, earliest_release_time=2020,
                               mingal_lat=0, myfilter_grating_combo="CLEAR;prism",
                               TARGET_TYPE="EITHER", RUN_MODE="QUERY",
                               FINAL_PRODUCT_MANIFEST="checkpoints/final_products.ecsv",
                               PRODUCT_LIST_CHECKPOINT_DIR="checkpoints/product_lists",
                               RUN_METRICS_DIR="run_metrics", FAILED_DOWNLOAD_LOG="checkpoints/failed.csv",
                               EXCLUDE_MIRROR_GRATING=True, transfer_metrics=None,
                               DIRECT_MAST_DOWNLOADS=False).items():
            self.stack.enter_context(patch.object(query, key, value))
        self.stack.enter_context(patch.object(raw, "SOURCE_POLICY", "AUTO"))
        self.stack.enter_context(patch.object(raw, "transfer_metrics", None))
        self.stack.enter_context(patch.object(query.time, "sleep"))

    def observations(self):
        return Table({"obsid": [123], "target_name": ["Alpha Centauri"],
                      "t_exptime": [2001], "t_obs_release": [60000],
                      "s_ra": [0.0], "s_dec": [60.0], "t_min": [60000]})

    def products(self):
        return Table({"productFilename": [NAME], "dataURI": ["mast:JWST/product/" + NAME],
                      "parent_obsid": [123], "calib_level": [1],
                      "productSubGroupDescription": ["UNCAL"],
                      "obs_id": ["test_nrs1"], "filters": ["CLEAR;PRISM"], "size": [6]})

    @staticmethod
    def download(product, download_dir, **kwargs):
        target = Path(download_dir) / "JWST" / NAME
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(b"abcdef")
        return Table({"Local Path": [str(target)]})

    def science_archive(self, fail=False):
        products = self.products()
        download = self.download
        class Client:
            def get_product_list(self, observations):
                return products.copy()
            def enable_cloud_dataset(self, **kwargs):
                pass
            def download_products(self, product, **kwargs):
                if fail:
                    raise requests.ReadTimeout("offline simulated timeout")
                return download(product, **kwargs)
        self.stack.enter_context(patch.object(query.Observations, "query_criteria", return_value=self.observations()))
        self.stack.enter_context(patch.object(query, "ObservationsClass", side_effect=Client))

    def test_science_success_writes_selection_files_and_both_metrics(self):
        self.science_archive()
        query.main("MAST", workers=2, direct_mast=True)
        self.assertEqual((self.root / "downloads/prism_observations" / NAME).read_bytes(), b"abcdef")
        selected = query.read_final_product_manifest()
        self.assertEqual(list(selected["productFilename"]), [NAME])
        run = read_csv(next(self.root.glob("run_metrics/*/summary.csv")))[0]
        self.assertEqual((run["selected_files"], run["downloaded_files"], run["downloaded_bytes"]), ("1", "1", "6"))
        transfer = read_csv(next(self.root.glob("transfer_metrics/*/summary.csv")))[0]
        self.assertEqual((transfer["source"], transfer["successes"], transfer["bytes"]), ("MAST", "1", "6"))
        self.assertEqual(read_csv(self.root / "checkpoints/failed.csv"), [])

    def test_download_only_reuses_manifest_and_files_without_archive_calls(self):
        query.write_final_product_manifest(self.products())
        target = self.root / "downloads/prism_observations" / NAME
        target.parent.mkdir(parents=True)
        target.write_bytes(b"abcdef")
        with patch.object(query, "RUN_MODE", "DOWNLOAD_ONLY"), \
                patch.object(query.Observations, "query_criteria", side_effect=AssertionError("query performed")), \
                patch.object(query, "ObservationsClass", side_effect=AssertionError("transfer performed")):
            query.main("MAST", direct_mast=True)
        run = read_csv(next(self.root.glob("run_metrics/*/summary.csv")))[0]
        self.assertEqual((run["run_mode"], run["downloaded_files"], run["downloaded_bytes"]), ("DOWNLOAD_ONLY", "0", "0"))
        transfer = read_csv(next(self.root.glob("transfer_metrics/*/summary.csv")))[0]
        self.assertEqual((transfer["skipped"], transfer["attempts"]), ("1", "0"))

    def test_science_failure_retains_manifest_and_failure_metrics(self):
        self.science_archive(fail=True)
        with patch.object(query, "AWS_DOWNLOAD_ATTEMPTS", 2):
            with self.assertRaisesRegex(RuntimeError, "1 of 1 MAST downloads failed"):
                query.main("MAST", direct_mast=True)
        self.assertTrue((self.root / "checkpoints/final_products.ecsv").exists())
        self.assertEqual(read_csv(self.root / "checkpoints/failed.csv")[0]["productFilename"], NAME)
        metrics = read_csv(next(self.root.glob("transfer_metrics/*/summary.csv")))[0]
        self.assertEqual((metrics["failures"], metrics["timeouts"], metrics["successes"]), ("2", "2", "0"))
        self.assertFalse((self.root / "downloads/prism_observations" / NAME).exists())

    def program_archive(self, failure=None):
        self.stack.enter_context(patch.object(pid.Observations, "query_criteria", return_value=self.observations()))
        self.stack.enter_context(patch.object(pid.Observations, "get_product_list", return_value=self.products()))
        return self.stack.enter_context(patch.object(pid.Observations, "download_products",
                                                     side_effect=failure or self.download))

    def test_pid_success_runs_selection_join_and_directory_layout(self):
        transfer = self.program_archive()
        result = pid.main(["4125", "MAST", "--output-root", "program"])
        self.assertEqual(result, 0)
        output = self.root / "program/PID04125"
        self.assertEqual((output / "Alpha_Centauri" / NAME).read_bytes(), b"abcdef")
        self.assertFalse(transfer.call_args.kwargs["cloud_only"])
        summary = read_csv(next(output.glob("transfer_metrics/*/summary.csv")))[0]
        self.assertEqual((summary["successes"], summary["bytes"]), ("1", "6"))

    def test_multiple_programs_have_separate_outputs_and_metrics(self):
        transfer = self.program_archive()
        with patch.object(pid.Observations, "query_criteria", return_value=self.observations()) as archive:
            self.assertEqual(pid.main(["4125", "3843", "4125", "MAST", "--output-root", "program"]), 0)
        self.assertEqual([call.kwargs["proposal_id"] for call in archive.call_args_list], ["4125", "3843"])
        self.assertEqual(transfer.call_count, 2)
        for number in [4125, 3843]:
            output = self.root / f"program/PID{number:05d}"
            self.assertEqual((output / "Alpha_Centauri" / NAME).read_bytes(), b"abcdef")
            summary = read_csv(next(output.glob("transfer_metrics/*/summary.csv")))[0]
            self.assertEqual((summary["successes"], summary["bytes"]), ("1", "6"))

    def test_default_list_runs_both_leggos_programs(self):
        self.program_archive()
        with patch("builtins.input", return_value="yes"), \
                patch.object(pid.Observations, "query_criteria", return_value=self.observations()) as archive:
            self.assertEqual(pid.main(["--source", "MAST", "--output-root", "program"]), 0)
        self.assertEqual([call.kwargs["proposal_id"] for call in archive.call_args_list], ["4125", "3843"])
        self.assertTrue((self.root / "program/PID04125/Alpha_Centauri" / NAME).exists())
        self.assertTrue((self.root / "program/PID03843/Alpha_Centauri" / NAME).exists())

    def test_empty_program_does_not_skip_remaining_programs(self):
        transfer = self.program_archive()
        with patch.object(pid.Observations, "query_criteria", side_effect=[self.observations()[:0], self.observations()]) as archive:
            self.assertEqual(pid.main(["4125", "3843", "MAST", "--output-root", "program"]), 0)
        self.assertEqual(archive.call_count, 2)
        self.assertEqual(transfer.call_count, 1)
        self.assertFalse((self.root / "program/PID04125").exists())
        self.assertTrue((self.root / "program/PID03843/Alpha_Centauri" / NAME).exists())

    def test_pid_failure_reports_manifest_error_as_failure(self):
        self.program_archive(lambda *args, **kwargs: Table({"Local Path": [""], "Status": ["ERROR"]}))
        with self.assertRaisesRegex(RuntimeError, "no local file"):
            pid.main(["4125", "MAST", "--output-root", "program"])
        summary = read_csv(next(self.root.glob("program/PID04125/transfer_metrics/*/summary.csv")))[0]
        self.assertEqual((summary["failures"], summary["successes"], summary["bytes"]), ("1", "0", "0"))

    def raw_archive(self, payload=None, failure=False):
        hdu = fits.PrimaryHDU()
        hdu.header["DATE-OBS"] = "2022-04-23"
        buffer = io.BytesIO()
        hdu.writeto(buffer)
        payload = buffer.getvalue() if payload is None else payload
        def request(session, method, url, **kwargs):
            if failure:
                raise requests.ReadTimeout("offline simulated timeout")
            return Response(payload)
        self.stack.enter_context(patch.object(requests.Session, "request", autospec=True, side_effect=request))
        names = self.root / "names.txt"
        names.write_text(NAME + "\n")
        return names, payload

    def test_filelist_success_reads_real_fits_date_and_writes_reports(self):
        names, payload = self.raw_archive()
        self.assertEqual(raw.main([str(names), "--source", "MAST", "--layout", "date-visit", "-o", "recovered"]), 0)
        target = self.root / "recovered/2022-04-23/VID01473016001" / NAME
        self.assertEqual(target.read_bytes(), payload)
        report = json.loads((self.root / "recovered/recovery_report.json").read_text())
        self.assertEqual(report["complete"], 1)
        self.assertEqual(report["missing_or_invalid"], [])
        self.assertEqual((self.root / "recovered/missing_files.txt").read_text(), "")
        rows = read_csv(next(self.root.glob("recovered/transfer_metrics/*/summary.csv")))
        transfer = next(row for row in rows if row["operation"] == "download")
        self.assertEqual(int(transfer["bytes"]), len(payload))
        self.assertEqual(transfer["successes"], "1")

    def test_filelist_exhaustion_and_gap_pass_produce_nonzero_status(self):
        names, _ = self.raw_archive(failure=True)
        self.assertEqual(raw.main([str(names), "--source", "MAST", "-o", "recovered"]), 1)
        report = json.loads((self.root / "recovered/recovery_report.json").read_text())
        self.assertEqual(report["complete"], 0)
        self.assertEqual(report["missing_or_invalid"], [NAME])
        self.assertIn(NAME, report["failures"])
        self.assertEqual((self.root / "recovered/missing_files.txt").read_text(), NAME + "\n")
        metrics = read_csv(next(self.root.glob("recovered/transfer_metrics/*/summary.csv")))[0]
        self.assertEqual((metrics["operation"], metrics["failures"], metrics["timeouts"]), ("HEAD", "10", "10"))

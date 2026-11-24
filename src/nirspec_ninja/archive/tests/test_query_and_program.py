"""Verify checkpoint reuse and common transfer metrics without archive access."""

import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from astropy.table import Table

from nirspec_ninja.archive import MAST_query_large_nirspec_datasets as query
from nirspec_ninja.archive import MAST_download_PID as program


NAME = "jw01473016001_02101_00001_nrs1_uncal.fits"


class QueryTests(unittest.TestCase):
    def test_parallel_checkpoints_are_reused_without_network_calls(self):
        observations = Table({"obsid": [101, 102, 103]})
        class Client:
            def get_product_list(self, chunk):
                number = int(chunk["obsid"][0])
                return Table({"dataURI": [str(number)], "productFilename": [f"{number}.fits"]})
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(query, "PRODUCT_LIST_CHECKPOINT_DIR", directory), \
                patch.object(query, "PRODUCT_LIST_CHUNK_SIZE", 1):
            with patch.object(query, "ObservationsClass", side_effect=Client) as clients:
                first = query.get_product_list_in_chunks(observations, workers=2)
            self.assertEqual(clients.call_count, 3)
            self.assertEqual(query.query_benchmark["product_chunks_fetched"], 3)
            with patch.object(query, "ObservationsClass", side_effect=AssertionError("network used")):
                resumed = query.get_product_list_in_chunks(observations, workers=2)
            self.assertEqual(list(first["dataURI"]), list(resumed["dataURI"]))
            self.assertEqual(query.query_benchmark["product_chunks_reused"], 3)
            self.assertEqual(query.query_benchmark["product_chunks_fetched"], 0)

    def test_direct_mast_routes_to_actual_transfers(self):
        with patch.object(query, "DIRECT_MAST_DOWNLOADS", True), \
                patch.object(query, "download_products_from_aws_parallel") as transfer:
            query.download_products_in_chunks(Table(), "output", "MAST", workers=2)
            self.assertEqual(transfer.call_args.args[2:], (2, "MAST"))


class ProgramTests(unittest.TestCase):
    def test_success_then_cached_file_is_excluded_from_rates(self):
        products = Table({"target_name": ["target"], "productFilename": [NAME], "size": [6]})
        def download(products, download_dir, **kwargs):
            path = Path(download_dir) / NAME
            path.write_bytes(b"abcdef")
            return Table({"Local Path": [str(path)]})
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch.object(program.Observations, "download_products", side_effect=download) as transfer:
                self.assertEqual(program.download_flat(products, output, "MAST"), 1)
                self.assertEqual(program.download_flat(products, output, "MAST"), 0)
            self.assertEqual(transfer.call_count, 1)
            summaries = []
            for path in output.glob("transfer_metrics/*/summary.csv"):
                with path.open() as stream:
                    summaries.extend(csv.DictReader(stream))
            self.assertEqual(sum(int(row["bytes"]) for row in summaries), 6)
            self.assertEqual(sum(int(row["skipped"]) for row in summaries), 1)

    def test_timeout_is_persisted_before_failure_propagates(self):
        products = Table({"target_name": ["target"], "productFilename": [NAME], "size": [6]})
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory)
            with patch.object(program.Observations, "download_products", side_effect=TimeoutError("expired")):
                with self.assertRaises(TimeoutError):
                    program.download_flat(products, output, "AWS")
            path = next(output.glob("transfer_metrics/*/summary.csv"))
            with path.open() as stream:
                row = next(csv.DictReader(stream))
            self.assertEqual(row["timeouts"], "1")
            self.assertEqual(row["failures"], "1")


if __name__ == "__main__":
    unittest.main()

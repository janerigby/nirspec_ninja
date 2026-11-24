"""Protect scientific selection and PID/target mapping using synthetic tables."""

import contextlib
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table, MaskedColumn
from astropy.time import Time
import astropy.units as u

from nirspec_ninja.archive import MAST_download_PID as pid
from nirspec_ninja.archive import MAST_query_large_nirspec_datasets as query
from nirspec_ninja.archive.tests.helpers import NAME


class SelectionTests(unittest.TestCase):
    def setUp(self):
        stack = contextlib.ExitStack()
        self.addCleanup(stack.close)
        for key, value in dict(minexptime=2000, earliest_release_time=2022.5,
                               TARGET_TYPE="EITHER", EXCLUDE_MIRROR_GRATING=True).items():
            stack.enter_context(patch.object(query, key, value))

    def observations(self):
        coords = SkyCoord(l=[0, 0, 0, 0, 0, 0] * u.deg,
                          b=[60, 60, 60, 60, 10, -60] * u.deg, frame="galactic").icrs
        return Table({
            "obsid": [1, 2, 3, 4, 5, 6],
            "t_exptime": [2001, 2000, 1999, 2001, 2001, 2001],
            "t_obs_release": Time([2023, 2023, 2023, 2022.5, 2023, 2023], format="decimalyear").mjd,
            "s_ra": coords.ra.deg, "s_dec": coords.dec.deg,
            "t_min": [60000] * 6,
        })

    def products(self):
        return Table({
            "productFilename": ["nrs1_uncal", "nrs2_uncal", "mirror", "cal", "other"],
            "obs_id": ["obs_nrs1", "obs_nrs2", "obs_nrs1", "obs_nrs1", "obs_nrs1"],
            "calib_level": ["1", "1.0", "1", "2", "1"],
            "productSubGroupDescription": ["uncal", " UNCAL ", "UNCAL", "UNCAL", "CAL"],
            "filters": ["CLEAR;PRISM", "CLEAR;PRISM", "clear;mirror", "CLEAR", "CLEAR"],
        })

    def test_scientific_boundaries_and_prism_detector(self):
        observations = self.observations()
        for grating, expected in [("CLEAR;prism", ["nrs1_uncal"]),
                                  ("F070LP;G140M", ["nrs1_uncal", "nrs2_uncal"])]:
            with self.subTest(grating=grating), tempfile.TemporaryDirectory() as directory, \
                    contextlib.ExitStack() as stack:
                for key, value in dict(myfilter_grating_combo=grating, minexptime=2000,
                                       earliest_release_time=2022.5, mingal_lat=20,
                                       TARGET_TYPE="EITHER", EXCLUDE_MIRROR_GRATING=True,
                                       FINAL_PRODUCT_MANIFEST=str(Path(directory) / "final.ecsv")).items():
                    stack.enter_context(patch.object(query, key, value))
                stack.enter_context(patch.object(query.Observations, "query_criteria", return_value=observations.copy()))
                products = stack.enter_context(patch.object(query, "get_product_list_in_chunks", return_value=self.products()))
                selected = query.query_and_select_products(workers=2)
                self.assertEqual(list(products.call_args.args[0]["obsid"]), [1, 6])
                self.assertEqual(list(selected["productFilename"]), expected)
                saved = query.read_final_product_manifest()
                self.assertEqual(list(saved["productFilename"]), expected)
                self.assertEqual(saved.meta["filter_grating"], grating)

    def test_latitude_boundary_is_strict(self):
        observations = self.observations()[:1]
        # Use the exact same computed coordinate to avoid floating-point threshold ambiguity.
        boundary = abs(SkyCoord(ra=observations["s_ra"] * u.deg,
                                dec=observations["s_dec"] * u.deg,
                                obstime=Time(observations["t_min"], format="mjd")).galactic.b[0].deg)
        with patch.object(query, "mingal_lat", boundary), \
                patch.object(query.Observations, "query_criteria", return_value=observations), \
                patch.object(query, "get_product_list_in_chunks") as products, \
                patch.object(query, "write_final_product_manifest") as manifest:
            self.assertIsNone(query.query_and_select_products())
            products.assert_not_called()
            manifest.assert_not_called()

    def test_target_type_criteria_and_empty_selection(self):
        empty = self.observations()[:0]
        for target, flag in [("FIXED", False), ("MOVING", True), ("EITHER", None)]:
            with self.subTest(target=target), patch.object(query, "TARGET_TYPE", target), \
                    patch.object(query.Observations, "query_criteria", return_value=empty.copy()) as archive:
                self.assertIsNone(query.query_and_select_products())
                criteria = archive.call_args.kwargs
                self.assertEqual(criteria["dataRights"], "PUBLIC")
                if flag is None:
                    self.assertNotIn("mtFlag", criteria)
                else:
                    self.assertIs(criteria["mtFlag"], flag)

    def test_missing_product_metadata_fails_clearly(self):
        with self.assertRaisesRegex(RuntimeError, "productSubGroupDescription"):
            query.filter_level1_uncal_products(Table({"calib_level": [1]}))
        with self.assertRaisesRegex(RuntimeError, "filters"):
            query.exclude_mirror_products(Table({"calib_level": [1]}))


class PIDSelectionTests(unittest.TestCase):
    def test_program_and_instrument_criteria(self):
        for instrument, expected in [("nirspec", "NIRSPEC*"), ("miri", "MIRI*"),
                                     ("NIRSPEC/IFU", "NIRSPEC/IFU"), ("all", None)]:
            with self.subTest(instrument=instrument):
                criteria = pid.observation_query_criteria(4125, instrument)
                self.assertEqual(criteria["proposal_id"], "4125")
                self.assertEqual(criteria["obs_collection"], "JWST")
                self.assertEqual(criteria.get("instrument_name"), expected)

    def test_raw_selection_handles_types_and_masks(self):
        table = Table({"productFilename": [NAME, "b_UNCAL.FITS", "c_cal.fits", "d_uncal.fits", "e_uncal.fits"]})
        table["calib_level"] = MaskedColumn(["1", "1.0", "1", "2", "1"], mask=[False, False, False, False, True])
        self.assertEqual(set(pid.select_raw_fits(table)["productFilename"]), {NAME, "b_UNCAL.FITS"})
        with self.assertRaisesRegex(RuntimeError, "calib_level"):
            pid.select_raw_fits(Table({"productFilename": [NAME]}))

    def test_target_join_normalizes_ids_and_rejects_missing_names(self):
        observations = Table({"obsid": [123, 124], "target_name": ["Alpha", "Beta"]})
        products = Table({"parent_obsid": ["123.0", "124"]})
        self.assertEqual(list(pid.add_target_names(products, observations)["target_name"]), ["Alpha", "Beta"])
        self.assertNotIn("target_name", products.colnames)
        for invalid in [Table({"parent_obsid": [999]}),
                        Table({"parent_obsid": MaskedColumn([123], mask=[True])})]:
            with self.subTest(products=invalid), self.assertRaisesRegex(RuntimeError, "Could not match"):
                pid.add_target_names(invalid, observations)
        observations["target_name"] = MaskedColumn(["Alpha", "Beta"], mask=[True, False])
        with self.assertRaisesRegex(RuntimeError, "Could not match"):
            pid.add_target_names(products, observations)
        with self.assertRaisesRegex(RuntimeError, "parent_obsid"):
            pid.add_target_names(Table({"other": [123]}), observations)

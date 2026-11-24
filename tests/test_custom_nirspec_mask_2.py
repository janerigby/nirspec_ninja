"""Regression checks for UNCAL-based illumination masking (no CRDS access)."""
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace as NS
from unittest.mock import patch

import numpy as np
from astropy.io import fits
from jwst import datamodels
from nirspec_ninja.masks import custom_nirspec_mask_2 as masks


class IlluminationTests(unittest.TestCase):
    @staticmethod
    def mock_trace(model, mask, region, limits):
        mask[70:76, 10:40] = False
        masks._add_trace_search_region(region, np.ones((6, 30), bool), 10, 70, *limits)
        return [((10., 39.), (70., 75.))]

    def test_curved_trace_limits_and_dark_shutter_gaps(self):
        footprint = np.zeros((35, 80), bool)
        for x in range(80):
            y = 4 + int(round(.004 * x*x))
            footprint[y:y+3, x] = True
        region = np.zeros((120, 200), bool)
        masks._add_trace_search_region(region, footprint, 50, 30, 30, 3)
        self.assertFalse(region[:, :20].any())
        self.assertFalse(region[:, 160:].any())
        self.assertTrue(region[65, 140])  # follows the red endpoint slope
        self.assertFalse(region[34, 140])  # horizontal bbox extension would admit this
        self.assertFalse(region[44, 60])  # distant across the trace
        gap_trace = np.zeros((20, 40), bool)
        gap_trace[2:5] = True
        gap_trace[15:18] = True
        region[:] = False
        masks._add_trace_search_region(region, gap_trace, 50, 30, 30, 3)
        self.assertTrue(region[33, 100])
        self.assertFalse(region[40, 100])  # preserve the unilluminated shutter gap

    def test_zero_limits_and_subarray_clipping(self):
        region = np.zeros((20, 25), bool)
        footprint = np.ones((3, 1), bool)
        masks._add_trace_search_region(region, footprint, 12, 10, 0, 0)
        self.assertEqual(region.sum(), 3)
        region[:] = False
        masks._add_trace_search_region(region, footprint, 12, 10, 5, 0)
        self.assertTrue(region[11, 7])
        self.assertTrue(region[11, 17])
        self.assertFalse(region[11, 6])
        self.assertFalse(region[11, 18])
        region[:] = False
        masks._add_trace_search_region(region, footprint, 0, 0, 30, 3)
        self.assertFalse(region[-1].any())
        self.assertTrue(region[0, -1])

    def test_ifu_cache_reuses_geometry_but_not_illumination(self):
        class LocalStep:
            def run(self, model):
                return NS(meta=NS(wcs=object()), close=lambda: None)

        def slits(model, mask, region, limits):
            for row in (10, 40, 70):
                mask[row:row+6, 10:40] = False
                masks._add_trace_search_region(region, np.ones((6, 30), bool), 10, row, *limits)
            return [((10., 39.), (10., 75.))]

        rng = np.random.default_rng(982)
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            cache = Path(directory) / 'geometry.npz'
            original_refine = masks._extend_mask_from_rate
            def refine_without_shared_lock(*args, **kwargs):
                # A separate descriptor must be able to take the same key lock
                # while exposure-specific illumination is being refined.
                import fcntl
                with open(str(cache) + '.lock', 'a') as handle:
                    fcntl.flock(handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
                    fcntl.flock(handle, fcntl.LOCK_UN)
                return original_refine(*args, **kwargs)
            with patch.object(masks, '_ifu_cache_path', return_value=cache), \
                    patch.object(masks, 'AssignWcsStep', LocalStep), \
                    patch.object(masks.nirspec, 'nrs_ifu_wcs', return_value=[]), \
                    patch.object(masks, '_extend_mask_from_rate', side_effect=refine_without_shared_lock), \
                    patch.object(masks, '_mask_slits_by_wcs', side_effect=slits) as project:
                outputs = []
                for i, row in enumerate((40, 70)):
                    data = rng.normal(0, .005, (100, 100)).astype(np.float32)
                    data[row:row+6, 40:60] += .3
                    with datamodels.IFUImageModel(data=data) as image:
                        image.meta.exposure.type = 'NRS_IFU'
                        image.meta.instrument.detector = 'NRS1'
                        image.meta.subarray.xsize = image.meta.subarray.ysize = 100
                        outputs.append(masks.make_custom_NIRSpec_mask(
                            image, Path(directory) / f'mask{i}.fits',
                            mask_FS_by_wcs=True, mask_FS_strip=False,
                            mask_failed_open=False, mask_by_illumination=True,
                            ifu_geometry_cache_dir=directory))
                self.assertEqual(project.call_count, 1)
                self.assertFalse(outputs[0][42, 50])
                self.assertTrue(outputs[1][42, 50])
                self.assertTrue(outputs[0][72, 50])
                self.assertFalse(outputs[1][72, 50])
                self.assertFalse(outputs[0][12, 20])
                self.assertFalse(outputs[1][12, 20])
                with np.load(cache) as saved:
                    self.assertTrue(saved['mask'][42, 50])
                    self.assertTrue(saved['mask'][72, 50])

    def test_ifu_rate_fits_illumination_and_grid_validation(self):
        class LocalStep:
            def run(self, model):
                return NS(meta=NS(wcs=object()), close=lambda: None)

        rng = np.random.default_rng(124)
        data = rng.normal(0, .005, (100, 100)).astype(np.float32)
        data[70:76, 40:60] += .3
        with datamodels.IFUImageModel(data=data) as image:
            image.meta.exposure.type = 'NRS_IFU'
            image.meta.instrument.detector = 'NRS1'
            image.meta.subarray.xsize = image.meta.subarray.ysize = 100
            with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory, \
                    patch.object(masks, 'AssignWcsStep', LocalStep), \
                    patch.object(masks, '_mask_slits_by_wcs', side_effect=self.mock_trace):
                rate = Path(directory) / 'rate.fits'
                output = Path(directory) / 'mask.fits'
                image.save(str(rate))
                with datamodels.open(rate) as reopened:
                    self.assertIsInstance(reopened, datamodels.IFUImageModel)
                    self.assertNotIsInstance(reopened, datamodels.ImageModel)
                result = masks.make_custom_NIRSpec_mask(
                    rate, output, mask_by_illumination=True,
                    mask_IFU=False, mask_FS_by_wcs=True, mask_FS_strip=False, mask_failed_open=False)
                self.assertFalse(result[72, 50])
                self.assertTrue(result[68, 75])
                wrong_grid = Path(directory) / 'wrong_grid.fits'
                with datamodels.IFUImageModel(data=np.zeros((99, 100), np.float32)) as wrong:
                    wrong.save(str(wrong_grid))
                with self.assertRaisesRegex(ValueError, r'IFUImageModel.*99, 100.*100, 100'):
                    masks.make_custom_NIRSpec_mask(
                        rate, Path(directory) / 'invalid.fits',
                        mask_illumination_reference=wrong_grid,
                        mask_IFU=False, mask_FS_by_wcs=True, mask_FS_strip=False, mask_failed_open=False)

    def test_unsigned_ramp_offsets_and_jump(self):
        shape = (100, 120)
        signal = np.zeros(shape, np.float32)
        signal[40:48, 50:80] = .2
        pedestal = 30000 + np.arange(120)[None, :] + np.zeros(shape)
        raw = np.empty((2, 9, *shape), np.uint16)
        for integration in range(2):
            for group in range(9):
                raw[integration, group] = np.round(
                    pedestal + group * 10 * (
                        signal + .03 * np.sin(np.arange(120) / 7)[None, :]
                        + .02 * np.cos(np.arange(100) / 9)[:, None]))
        raw[0, 3:, 20, 20] += 1000
        meta = NS(exposure=NS(group_time=10, readpatt='NRSRAPID'),
                  subarray=NS(slowaxis=1))
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
            path = Path(directory) / 'uncal.fits'
            fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(raw, name='SCI')]).writeto(path)
            result = masks._ramp_illumination_proxy(path, signal == 0, meta)
        self.assertAlmostEqual(float(np.median(result[40:48, 50:80])), .2, delta=.01)
        self.assertLess(abs(float(result[20, 20])), .03)
        self.assertLess(np.percentile(np.abs(result[signal == 0]), 99), .03)
        short = NS(data=np.zeros((1, 2, *shape), np.float32))
        with self.assertRaisesRegex(ValueError, 'three groups'):
            masks._ramp_illumination_proxy(short, signal == 0, meta)

    def test_local_illumination_and_distant_glow(self):
        rng = np.random.default_rng(124)
        proxy = rng.normal(0, .005, (100, 100))
        proxy[10:16, 40:65] += .3
        proxy[70:76, 60:90] += 10
        proxy[30, 30] = 100
        proxy[90, 90] = np.nan
        mask = np.ones(proxy.shape, bool)
        mask[10:16, 10:40] = False
        region = np.zeros(mask.shape, bool)
        masks._add_trace_search_region(region, np.ones((6, 30), bool), 10, 10, 30, 3)
        masks._extend_mask_from_rate(mask, proxy, np.isfinite(proxy), [], None, 5, region)
        self.assertFalse(mask[12, 50])
        self.assertTrue(mask[72, 75])
        self.assertTrue(mask[68, 75])
        self.assertTrue(mask[77, 75])
        self.assertTrue(mask[30, 30])
        self.assertTrue(mask[90, 90])  # invalid pixels cannot seed distant exclusions
        self.assertFalse(mask[12, 20])

    def test_rate_input_ownership_and_option_validation(self):
        rng = np.random.default_rng(124)
        data = rng.normal(0, .005, (100, 100)).astype(np.float32)
        data[70:76, 40:60] += .3
        image = datamodels.ImageModel(data=data.copy())
        image.meta.exposure.type = 'NRS_FIXEDSLIT'
        image.meta.instrument.detector = 'NRS1'
        image.meta.subarray.xsize = image.meta.subarray.ysize = 100

        class LocalStep:
            def run(self, model):
                return NS(meta=NS(wcs=object()), close=lambda: None)

        def geometry(model, mask, region, limits):
            mask[70:76, 10:40] = False
            masks._add_trace_search_region(region, np.ones((6, 30), bool), 10, 70, *limits)
            return [((10., 39.), (70., 75.))]

        try:
            with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as directory:
                output = Path(directory) / 'mask.fits'
                with patch.object(masks, 'AssignWcsStep', LocalStep), patch.object(
                        masks, '_mask_slits_by_wcs', side_effect=geometry):
                    result = masks.make_custom_NIRSpec_mask(
                        image, output, mask_FS_by_wcs=True,
                        mask_by_illumination=True, mask_failed_open=False)
                self.assertFalse(result[72, 50])
                self.assertTrue(result[68, 75])
                np.testing.assert_array_equal(image.data, data)
                with datamodels.open(output) as saved:
                    np.testing.assert_array_equal(saved.data, result.astype(np.float32))
                with self.assertRaisesRegex(ValueError, 'not both'):
                    masks.make_custom_NIRSpec_mask(
                        image, Path(directory) / 'invalid.fits', mask_by_illumination=True,
                        mask_illumination_reference='unused.fits')
        finally:
            image.close()


if __name__ == '__main__':
    unittest.main()

"""Mask-only MSA normalization; no CRDS downloads or cleaning."""
from contextlib import ExitStack
from pathlib import Path
import tempfile
import unittest
from types import SimpleNamespace as NS
from unittest.mock import patch

import numpy as np
from astropy.io import fits
from jwst import datamodels
from nirspec_ninja.masks import custom_nirspec_mask_2 as masks
from nirspec_ninja.masks.custom_nirspec_mask_2 import _mask_msa_metadata


class CrossQuadrantTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.path = Path(self.directory.name) / 'msa.fits'
        dtype = [('slitlet_id', 'i2'), ('msa_metadata_id', 'i2'),
                 ('dither_point_index', 'i2'), ('shutter_quadrant', 'i2'),
                 ('shutter_row', 'i2'), ('shutter_column', 'i2'),
                 ('primary_source', 'U1'), ('background', 'U1'),
                 ('fixed_slit', 'U8')]
        self.table = np.array([
            (10, 1, 1, 3, 116, 170, 'N', 'Y', 'NONE'),
            (10, 1, 1, 3, 116, 171, 'Y', 'N', 'NONE'),
            (10, 1, 1, 4, 116, 1, 'N', 'Y', 'NONE'),
            (11, 1, 1, 3, 200, 50, 'Y', 'N', 'NONE'),
            (10, 1, 2, 4, 116, 1, 'Y', 'N', 'NONE'),
            (10, 1, 2, 4, 116, 2, 'N', 'Y', 'NONE'),
        ], dtype=dtype)
        self.image = NS(meta=NS(instrument=NS(msa_metadata_file='msa.fits',
                            msa_metadata_id=1, grating='PRISM'),
                            dither=NS(position_number=1)))
        self.write()

    def write(self):
        fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU(self.table,
                      name='SHUTTER_INFO')]).writeto(self.path, overwrite=True)

    def test_split_is_private_and_preserves_other_records(self):
        before = self.path.read_bytes()
        with ExitStack() as stack:
            target = _mask_msa_metadata(self.image, self.directory.name, stack)
            table = fits.getdata(target, 'SHUTTER_INFO')
            np.testing.assert_array_equal(table[:2]['slitlet_id'], [10, 10])
            self.assertNotIn(int(table[2]['slitlet_id']), [10, 11])
            for name in self.table.dtype.names:
                np.testing.assert_array_equal(table[3:][name], self.table[3:][name])
            self.assertEqual(table[2]['shutter_quadrant'], 4)
            self.assertEqual(table[2]['shutter_column'], 1)
            self.assertEqual(self.path.read_bytes(), before)
        self.assertFalse(Path(target).exists())

    def test_unaffected_dither_bypasses_normalization(self):
        self.image.meta.dither.position_number = 2
        with ExitStack() as stack:
            self.assertIsNone(_mask_msa_metadata(self.image, self.directory.name, stack))

    def test_grating_warns_without_modifying_geometry(self):
        self.image.meta.instrument.grating = 'G235M'
        before = self.path.read_bytes()
        with ExitStack() as stack, self.assertLogs(masks.__name__, level='WARNING') as messages:
            self.assertIsNone(_mask_msa_metadata(self.image, self.directory.name, stack))
        self.assertIn('retaining existing grating', messages.output[0])
        self.assertEqual(self.path.read_bytes(), before)

    def test_ambiguous_primary_is_rejected(self):
        self.table['primary_source'][2] = 'Y'
        self.write()
        with ExitStack() as stack, self.assertRaisesRegex(ValueError, 'more than one primary'):
            _mask_msa_metadata(self.image, self.directory.name, stack)

    def test_background_only_crossing_is_supported(self):
        self.table['primary_source'][:3] = 'N'
        self.table['background'][:3] = 'Y'
        self.write()
        with ExitStack() as stack:
            target = _mask_msa_metadata(self.image, self.directory.name, stack)
            table = fits.getdata(target, 'SHUTTER_INFO')
            self.assertEqual(len(np.unique(table['slitlet_id'][:3])), 2)
            self.assertTrue((table['primary_source'][:3] == 'N').all())

    def test_ifu_does_not_read_or_normalize_msa_metadata(self):
        class LocalStep:
            def run(self, model):
                return NS(meta=NS(wcs=object()), close=lambda: None)

        with datamodels.IFUImageModel(data=np.zeros((20, 20), np.float32)) as image:
            image.meta.exposure.type = 'NRS_IFU'
            image.meta.instrument.detector = 'NRS1'
            image.meta.instrument.grating = 'PRISM'
            image.meta.instrument.msa_metadata_file = 'nonexistent.fits'
            with patch.object(masks, 'AssignWcsStep', LocalStep), \
                    patch.object(masks.nirspec, 'nrs_ifu_wcs', return_value=[]), \
                    patch.object(masks, '_mask_msa_metadata') as normalize:
                result = masks.make_custom_NIRSpec_mask(
                    image, Path(self.directory.name) / 'ifu_mask.fits',
                    mask_failed_open=False, mask_FS_strip=False)
                normalize.assert_not_called()
                self.assertTrue(result.all())


if __name__ == '__main__':
    unittest.main()

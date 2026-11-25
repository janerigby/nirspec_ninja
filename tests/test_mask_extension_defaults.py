"""Exercise the public builder's default continuation and explicit off switch."""
import json
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from jwst import datamodels
from nirspec_ninja.masks import custom_nirspec_mask_2 as masks

class ExtensionDefaultTests(unittest.TestCase):
    def test_full_aperture_continues_only_when_enabled(self):
        class SliceWcs:
            bounding_box = ((20, 70), (30, 37))
            def __call__(self, x, y):
                return x, y, np.ones_like(x)
        class LocalStep:
            def run(self, model):
                return SimpleNamespace(meta=SimpleNamespace(wcs=object()), close=lambda: None)
        cases = [
            ('G235M', 'F170LP', None, True),
            ('G235H', 'F170LP', None, True),
            ('G140H', 'F100LP', None, True),
            ('G235H', 'F170LP', False, False),
            ('G395H', 'F290LP', None, False),
            ('G395H', 'F290LP', True, True),
            ('PRISM', 'CLEAR', True, False),
        ]
        with tempfile.TemporaryDirectory() as directory, \
                patch.object(masks, 'AssignWcsStep', LocalStep), \
                patch.object(masks.nirspec, 'nrs_ifu_wcs', return_value=[SliceWcs()]):
            for i, (grating, filter_name, option, expected) in enumerate(cases):
                with self.subTest(grating=grating, option=option), \
                        datamodels.ImageModel(data=np.zeros((100, 100), np.float32)) as image:
                    image.meta.exposure.type = 'NRS_IFU'
                    image.meta.instrument.detector = 'NRS2'
                    image.meta.instrument.grating = grating
                    image.meta.instrument.filter = filter_name
                    image.meta.subarray.xsize = image.meta.subarray.ysize = 100
                    output = Path(directory) / f'mask{i}.fits'
                    kwargs = {} if option is None else {'mask_grating_extrapolate': option}
                    result = masks.make_custom_NIRSpec_mask(
                        image, output, mask_failed_open=False, mask_FS_strip=False, **kwargs)
                    self.assertFalse(result[30:38, 40].any())
                    self.assertEqual(bool((~result[30:38, 90]).all()), expected)
                    self.assertTrue(result[55:65, 90].all())
                    self.assertTrue(result[30:38, 10].all())
                    report = output.with_suffix('.extrapolation.json')
                    self.assertEqual(report.exists(), expected)
                    if expected:
                        self.assertTrue(json.loads(report.read_text())['fits'])

if __name__ == '__main__':
    unittest.main()

"""Offline finder chart generation with a known celestial-WCS image."""
import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from nirspec_ninja.finder import nirspec_finder as finder

class FinderTests(unittest.TestCase):
    def _check_chart(self, image_pixels, fszoom=False):
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp); exposure = root/'exposure.fits'
            hdu = fits.PrimaryHDU()
            hdu.header.update(INSTRUME='NIRSPEC', RA_V1=150., DEC_V1=20., PA_V3=30., EXP_TYPE='NRS_IFU')
            hdu.writeto(exposure)
            footprints, center, size, _, _, _ = finder.geometry(exposure)
            self.assertEqual(set(footprints), set(finder.NAMES))
            self.assertGreater(size, 0)
            wcs = WCS(naxis=2)
            wcs.wcs.ctype = ['RA---TAN', 'DEC--TAN']
            wcs.wcs.crval = [center.ra.deg, center.dec.deg]
            wcs.wcs.crpix = [(image_pixels + 1) / 2, (image_pixels + 1) / 2]
            wcs.wcs.cdelt = [-.001, .001]
            image = root/'image.fits'; output = root/'finder.png'
            fits.PrimaryHDU(np.random.default_rng(3).normal(size=(image_pixels,image_pixels)), wcs.to_header()).writeto(image)
            with patch('sys.argv', ['nirspec_finder', str(exposure), '--image', str(image), '--output', str(output), '--regions'] + (['--FSzoom'] if fszoom else [])), patch.object(finder, 'wise_image', side_effect=AssertionError('network forbidden')):
                finder.main()
            self.assertGreater(output.stat().st_size, 1000)
            self.assertEqual(output.with_suffix('.reg').read_text().count('polygon('), len(finder.NAMES))
            record = json.loads(output.with_suffix('.json').read_text())
            self.assertEqual(record['image'], str(image))
            self.assertEqual(record['fiducial']['aperture'], 'S200A1')
            if fszoom:
                self.assertTrue(record['FSzoom'])
                self.assertAlmostEqual(record['size_arcmin'], 100/60)
                self.assertAlmostEqual(record['height_arcmin'], 1)

    def test_local_chart_and_region_outputs(self):
        self._check_chart(400)

    def test_partial_image_still_shows_all_apertures(self):
        self._check_chart(40)

    def test_fixed_slit_zoom(self):
        self._check_chart(400, fszoom=True)

    def test_missing_attitude_is_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp)/'bad.fits'; hdu = fits.PrimaryHDU(); hdu.header['INSTRUME']='NIRSPEC'; hdu.writeto(path)
            with self.assertRaisesRegex(ValueError, 'Missing attitude'):
                finder.geometry(path)

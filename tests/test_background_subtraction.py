"""Subtracting a known supplied sky spectrum preserves the source signal."""
import tempfile
import unittest
from pathlib import Path
import numpy as np
from astropy.io import fits
from nirspec_ninja.subtract_background.subtraction import exp_bg_sub

class SuppliedBackgroundTests(unittest.TestCase):
    def test_known_cube_background_and_saved_result(self):
        with tempfile.TemporaryDirectory() as tmp:
            cube = Path(tmp) / 'cube.fits'; bg = Path(tmp) / 'sky.txt'; out = Path(tmp) / 'subtracted.fits'
            wave = np.arange(1., 5.); sky = 2 * wave
            data = np.broadcast_to((sky + 7)[:, None, None], (4, 3, 2)).copy()
            sci = fits.ImageHDU(data, name='SCI'); sci.header.update(CRVAL3=1., CDELT3=1.)
            fits.HDUList([fits.PrimaryHDU(), sci]).writeto(cube)
            np.savetxt(bg, np.column_stack([wave, sky]))
            np.testing.assert_allclose(exp_bg_sub(cube, bg), 7.)
            exp_bg_sub(cube, bg, out)
            np.testing.assert_allclose(fits.getdata(out, 'SCI'), 7.)
            np.testing.assert_array_equal(fits.getdata(cube, 'SCI'), data)

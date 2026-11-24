import csv
import re
import tempfile
import unittest
from pathlib import Path
import numpy as np
from astropy.io import fits
from nirspec_ninja.diagnostics.flicker_noise_diagnostics import measure_noise, write_noise_report


class NoiseTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.root = Path(self.tmp.name)
        # Large illuminated signal, invalid finite/DQ values, and DC offsets
        # must not affect the noise of the four valid background samples.
        a = np.array([[-3, -1, 1, 3], [900, np.nan, 999, 999]], dtype=float)
        b = np.array([[-1.5, -.5, .5, 1.5], [800, 99, 999, 999]]) + 10
        dq = np.array([[0, 0, 0, 0], [0, 0, 1, 512]], dtype=np.uint32)
        for name, data in [('before', a), ('example_rate', b)]:
            sci = fits.ImageHDU(data, name='SCI')
            sci.header['BUNIT'] = 'DN/s'
            fits.HDUList([fits.PrimaryHDU(), sci, fits.ImageHDU(dq, name='DQ')]).writeto(self.root / (name + '.fits'))
        fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(
            np.array([[1, 1, 1, 1], [0, 1, 1, 1]]), name='SCI')]).writeto(self.root / 'mask.fits')
        self.args = [self.root / n for n in ('before.fits', 'example_rate.fits', 'mask.fits')]

    def test_common_pixels_and_noise_reduction(self):
        row = measure_noise(*self.args)
        self.assertEqual(row['N_pixels'], 4)
        self.assertAlmostEqual(row['RMS_before'], np.sqrt(5))
        self.assertAlmostEqual(row['RMS_after'], np.sqrt(5) / 2)
        self.assertAlmostEqual(row['MAD_sigma_before'], 2 * 1.4826022185)
        self.assertAlmostEqual(row['RMS_reduction_percent'], 50)

    def test_report_upsert_keeps_other_exposures(self):
        write_noise_report(*self.args)
        other = self.root / 'other_rate.fits'
        other.write_bytes(self.args[1].read_bytes())
        write_noise_report(self.args[0], other, self.args[2])
        write_noise_report(*self.args)
        with (self.root / 'flicker_noise_rms.txt').open() as stream:
            reader = csv.DictReader(('\t'.join(re.split(r' {2,}', l.strip())) if '\t' not in l else l for l in stream if not l.startswith('#')), delimiter='\t')
            from nirspec_ninja.diagnostics.flicker_noise_diagnostics import HEADINGS
            reverse = {value: key for key, value in HEADINGS.items()}
            reader.fieldnames = [reverse.get(name.strip(), name.strip()) for name in reader.fieldnames]
            rows = list(reader)
        self.assertEqual([r['Rate_filename'].strip() for r in rows], ['example_rate.fits', 'other_rate.fits'])

    def test_report_metadata_precision_and_legacy_upgrade(self):
        with fits.open(self.args[1], mode='update') as hdus:
            hdus[0].header.update(EFFEXPTM=123.456789, READPATT='NRSIRS2',
                NGROUPS=15, NINTS=2, NFRAMES=1, GROUPGAP=0, TGROUP=14.5889,
                TFRAME=14.5889, DETECTOR='NRS1', SUBARRAY='FULL')
        report = self.root / 'flicker_noise_rms.txt'
        from nirspec_ninja.diagnostics.flicker_noise_diagnostics import LEGACY_COLUMNS
        original = measure_noise(*self.args)
        with report.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=LEGACY_COLUMNS, delimiter='\t')
            writer.writeheader()
            writer.writerow({key: original[key] for key in LEGACY_COLUMNS})
        row = write_noise_report(*self.args, fit_method='fft')
        self.assertEqual(row['Fit_method'], 'FFT')
        self.assertEqual(row['READPATT'], 'NRSIRS2')
        self.assertEqual(row['NGROUPS'], 15)
        self.assertAlmostEqual(row['MAD_reduction_percent'], 50)
        with report.open() as stream:
            reader = csv.DictReader(('\t'.join(re.split(r' {2,}', line.strip())) if '\t' not in line else line for line in stream if not line.startswith('#')), delimiter='\t')
            from nirspec_ninja.diagnostics.flicker_noise_diagnostics import HEADINGS
            reverse = {value: key for key, value in HEADINGS.items()}
            reader.fieldnames = [reverse.get(name.strip(), name.strip()) for name in reader.fieldnames]
            saved = {key: value.strip() for key, value in next(reader).items()}
        self.assertEqual(saved['RMS_before'], '2.236')
        self.assertEqual(saved['MAD_reduction_percent'], '50.00')
        self.assertEqual(saved['EFFEXPTM'], '123.457')

    def test_empty_mask_rejected(self):
        with fits.open(self.args[2], mode='update') as hdus:
            hdus['SCI'].data[:] = 0
        with self.assertRaisesRegex(ValueError, 'usable background'):
            measure_noise(*self.args)

    def test_unit_mismatch_rejected(self):
        with fits.open(self.args[1], mode='update') as hdus:
            hdus['SCI'].header['BUNIT'] = 'electrons/s'
        with self.assertRaisesRegex(ValueError, 'units differ'):
            measure_noise(*self.args)


if __name__ == '__main__':
    unittest.main()

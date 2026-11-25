"""Check extraction preserves a known spectrum and input arrays."""
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np
from astropy.io import fits
from nirspec_ninja.fixed_slit import extraction

class ExtractionTests(unittest.TestCase):
    def test_header_diagnostics_mask_a_detected_trace(self):
        wave = np.linspace(1, 5, 100)
        science = np.full((32, 100), 2.); science[12:18] = 20.
        with tempfile.TemporaryDirectory() as tmp:
            infile = Path(tmp) / 'synthetic_S200A1_s2d.fits'
            primary = fits.PrimaryHDU(); primary.header['EXPMID'] = 60000.
            sci = fits.ImageHDU(science, name='SCI'); sci.header['SLTNAME'] = 'S200A1'
            fits.HDUList([primary, sci, fits.ImageHDU(np.ones_like(science), name='ERR'),
                fits.ImageHDU(np.broadcast_to(wave, science.shape), name='WAVELENGTH')]).writeto(infile)
            result = dict(status='ok', significance=10., detected=True, bounds_flux_rows=[11., 18.])
            with patch.object(extraction, 'diagnostic_trace', return_value=result):
                frame = extraction.median_combine_level3_nirspecFS2(str(infile), 'S200A1', tmp + '/')
            self.assertEqual(frame.attrs['extraction_report']['spatial_flatness']['fractional_departure'], 0.)
            self.assertTrue(frame.attrs['extraction_report']['spatial_flatness']['trace_masked'])
            text = (Path(tmp) / 'synthetic_S200A1_s2d_extract.csv').read_text()
            self.assertIn('##trace_sig2 10.0', text)
            self.assertIn('##non-flat 0.0', text)
            np.testing.assert_array_equal(fits.getdata(infile, 'SCI'), science)

    def test_cal_diagnostics_use_the_combined_2D_profile(self):
        wave = np.linspace(1, 5, 100)
        image = np.full((32, 100), 2.)
        trace, flatness = extraction.measure_extraction_diagnostics('unused_cal.fits', 'S200A1', 'CAL',
            image, wave, 1, clip_nan=3, detect_traces=False)
        self.assertEqual(trace['status'], 'skipped')
        self.assertEqual(flatness['fractional_departure'], 0.)

    def test_synthetic_s2d_spectrum_and_new_detector(self):
        wave = np.linspace(1, 5, 40)
        science = 2 + 3 * wave[None, :] + (np.arange(32) - 15.5)[:, None] * .01
        with tempfile.TemporaryDirectory() as tmp:
            infile = Path(tmp) / 'synthetic_S200A1_s2d.fits'
            primary = fits.PrimaryHDU(); primary.header['EXPMID'] = 60000.
            sci = fits.ImageHDU(science, name='SCI'); sci.header['SLTNAME'] = 'S200A1'
            fits.HDUList([primary, sci, fits.ImageHDU(np.ones_like(science), name='ERR'),
                fits.ImageHDU(np.broadcast_to(wave, science.shape), name='WAVELENGTH')]).writeto(infile)
            with patch.object(extraction, 'diagnostic_trace', return_value={'status': 'ok', 'significance': 7.5, 'detected': False}) as detector:
                frame = extraction.median_combine_level3_nirspecFS2(str(infile), 'S200A1', tmp + '/')
            detector.assert_called_once()
            self.assertEqual(detector.call_args.kwargs['clip_nan'], 5)
            self.assertTrue(detector.call_args.kwargs['use_error'])
            np.testing.assert_allclose(frame['wave'], wave)
            np.testing.assert_allclose(frame['fnu'], 2 + 3 * wave, atol=1e-12)
            np.testing.assert_array_equal(fits.getdata(infile, 'SCI'), science)
            output = (Path(tmp) / 'synthetic_S200A1_s2d_extract.csv').read_text()
            self.assertIn('##trace_sig2 7.5', output)
            self.assertIn('##non-flat ', output)
            self.assertIn('##non-flat_status ok', output)
            self.assertAlmostEqual(100 * frame.attrs['extraction_report']['spatial_flatness']['fractional_departure'],
                                   100 * .085 / 11., places=8)
            self.assertNotIn('##trace_sig ', output)
            with patch.object(extraction, 'diagnostic_trace') as detector:
                skipped = extraction.median_combine_level3_nirspecFS2(str(infile), 'S200A1', tmp + '/', detect_traces=False)
            detector.assert_not_called()
            np.testing.assert_array_equal(skipped['fnu'], frame['fnu'])
    def test_trace_failure_still_writes_csv_and_reports_path(self):
        wave = np.linspace(1, 5, 40)
        science = np.broadcast_to(wave, (32, 40)).copy()
        with tempfile.TemporaryDirectory() as tmp:
            infile = Path(tmp) / 'synthetic_s2d.fits'
            primary = fits.PrimaryHDU(); primary.header['EXPMID'] = 60000.
            sci = fits.ImageHDU(science, name='SCI'); sci.header['SLTNAME'] = 'S200A1'
            fits.HDUList([primary, sci, fits.ImageHDU(np.ones_like(science), name='ERR'),
                fits.ImageHDU(science, name='WAVELENGTH')]).writeto(infile)
            failure = {'status': 'failed', 'reason': 'no usable rows', 'profile': np.full(32, np.nan)}
            with patch.object(extraction, 'diagnostic_trace', return_value=failure), patch('builtins.print') as output:
                reports = extraction.wrap_median_combine_level3_nirspecFS(tmp, tmp + '/', return_report=True)
            self.assertEqual(reports[0]['status'], 'ok')
            self.assertEqual(reports[0]['trace']['status'], 'failed')
            self.assertTrue(Path(reports[0]['output']).exists())
            self.assertNotIn('##trace_sig2 ', Path(reports[0]['output']).read_text())
            message = output.call_args.args[0]
            self.assertIn(str(infile), message)
            self.assertIn('CSV written', message)
            self.assertNotIn('array(', message)

    def test_reporting_continues_after_file_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            for name in ('a_s2d.fits', 'b_s2d.fits'):
                (Path(tmp) / name).touch()
            frame = __import__('pandas').DataFrame()
            frame.attrs['extraction_report'] = {'status': 'ok'}
            with patch.object(extraction, 'gethead', return_value='S200A1'), patch.object(
                    extraction, 'median_combine_level3_nirspecFS2', side_effect=[ValueError('bad input'), frame]), patch('builtins.print'):
                reports = extraction.wrap_median_combine_level3_nirspecFS(tmp, tmp + '/', return_report=True)
            self.assertEqual([r['status'] for r in reports], ['failed', 'ok'])
            self.assertIn('a_s2d.fits', reports[0]['input'])
            self.assertIn('bad input', reports[0]['reason'])

    def test_nan_rows_and_neighbors_are_excluded(self):
        science = np.ones((20, 10)); science[10, :5] = np.nan
        _, _, excluded = extraction._nirspecfs_nan_mask_rows(science, 2)
        self.assertTrue(excluded[8:13].all())
        self.assertTrue(excluded[:3].all())
        self.assertTrue(excluded[-3:].all())
        self.assertFalse(excluded[5])

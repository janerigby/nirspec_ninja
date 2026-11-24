"""Human slit selection and diagnostic-only S1600 support."""
import tempfile
import unittest
from pathlib import Path
from nirspec_ninja.diagnostics.fs_gradients import read_crowding_notes, GRADIENT_SLITS
from nirspec_ninja.subtract_background.fs_background import USABLE_SLITS
import numpy as np
from nirspec_ninja.diagnostics.fs_gradients import median_good_profiles
from nirspec_ninja.diagnostics.fs_gradients import combined_median_spectrum, spatial_median_spectrum


class CrowdingNotesTests(unittest.TestCase):
    def read(self, text, target='COSMIC.EYE'):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'notes.txt'
            path.write_text(text)
            return read_crowding_notes(path, target)

    def test_template_order_case_and_notes(self):
        marks, notes = self.read(
            '# 1 is use, 0 is do not use\n'
            'target S200A2 S200A1 S400A1 S1600A1 Notes\n'
            'Cosmic.eye 1 1 0 1 Galaxy in S400A1\n')
        self.assertEqual(marks, dict(S200A2=True, S200A1=True, S400A1=False, S1600A1=True))
        self.assertEqual(notes, 'Galaxy in S400A1')

    def test_optional_and_unlisted_targets(self):
        self.assertEqual(read_crowding_notes(None, 'anything'), ({}, ''))
        self.assertEqual(self.read('target S200A1 Notes\nOther 0 crowded\n'), ({}, ''))

    def test_invalid_mark_and_duplicate_target(self):
        for text in ('target S200A1\nCosmic.eye 2\n',
                     'target S200A1\nCosmic.eye 1\nCOSMICEYE 0\n'):
            with self.assertRaises(ValueError):
                self.read(text)

    def test_s1600_is_diagnostic_only(self):
        self.assertIn('S1600A1', GRADIENT_SLITS)
        self.assertNotIn('S1600A1', USABLE_SLITS)

    def test_candidate_median_excludes_bad_and_unmarked(self):
        profiles = [('S200A1', np.array([0., 1.]), np.array([1., 3.])),
                    ('S200A2', np.array([0., .5, 1.]), np.array([3., 4., 5.])),
                    ('S400A1', np.array([0., 1.]), np.array([100., 100.])),
                    ('S1600A1', np.array([0., 1.]), np.array([200., 200.]))]
        grid, median, names = median_good_profiles(profiles, dict(S200A1=True, S200A2=True, S400A1=False))
        np.testing.assert_allclose(median, 3.)
        self.assertEqual(names, ['S200A1', 'S200A2'])
        self.assertIsNone(median_good_profiles(profiles, {})[1])

    def test_finder_chart_order(self):
        self.assertEqual(GRADIENT_SLITS, ('S200A2', 'S400A1', 'S1600A1', 'S200A1'))

    def test_combined_spectrum_common_coverage_and_gaps(self):
        spectra = [('A', np.arange(6.), np.array([1., 2., np.nan, 4., 5., 6.])),
                   ('B', np.arange(1., 7.), np.arange(4., 10.))]
        wave, values = combined_median_spectrum(spectra)
        np.testing.assert_array_equal(wave, np.arange(1., 6.))
        np.testing.assert_allclose(values, [3., np.nan, 5., 6., 7.], equal_nan=True)

    def test_spatial_median_recovers_spectrum(self):
        from unittest.mock import patch
        from nirspec_ninja.subtract_background.fs_background import SlitPixels
        wave = np.broadcast_to(np.arange(10.), (6, 10)).copy()
        brightness = wave + np.arange(6.)[:, None]
        slit = SlitPixels('S200A1', brightness, wave, np.ones_like(wave), None)
        with patch('nirspec_ninja.diagnostics.fs_gradients.bg.spatial_pixel_mask',
                   return_value=(np.ones_like(wave, dtype=bool), np.zeros(6, dtype=bool))):
            x, values = spatial_median_spectrum(slit)
        np.testing.assert_allclose(values, x + 2.5)

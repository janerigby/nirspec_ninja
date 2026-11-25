"""Restrict candidate source centers without changing the default search."""
import unittest
import numpy as np
from nirspec_ninja.fixed_slit.source_detection import detect_source_s2d

class SourceCenterTests(unittest.TestCase):
    def test_requested_widths_cannot_exceed_slit_fwhm(self):
        y = np.arange(32)[:, None]
        data = np.broadcast_to(2 + np.exp(-.5 * ((y - 16) / 8)**2), (32, 80)).copy()
        rejected = detect_source_s2d(data, error_2d=np.full_like(data, .1), slit_length_pixels=16, widths=[8.])
        self.assertEqual(rejected['status'], 'failed')
        result = detect_source_s2d(data, error_2d=np.full_like(data, .1), slit_length_pixels=16, widths=[2., 6.])
        self.assertEqual(result['status'], 'ok')
        self.assertLessEqual(result['fwhm'], 16.)
    def test_allowed_centers_and_default_behavior(self):
        rng = np.random.default_rng(19)
        y = np.arange(32)[:, None]
        data = rng.normal(0, .1, (32, 80)) + 5 * np.exp(-.5 * ((y - 10) / 2)**2)
        original = detect_source_s2d(data)
        explicit = detect_source_s2d(data, allowed_source_center_rows=None)
        self.assertEqual(original['center'], explicit['center'])
        self.assertEqual(original['significance'], explicit['significance'])
        allowed = np.zeros(32, dtype=bool); allowed[16:26] = True
        result = detect_source_s2d(data, allowed_source_center_rows=allowed)
        self.assertEqual(result['status'], 'ok')
        self.assertTrue(allowed[int(result['center'])])
    def test_invalid_or_empty_masks_report_failure(self):
        data = np.ones((32, 80))
        for mask in (np.ones(31, bool), np.ones(32, int), np.zeros(32, bool)):
            result = detect_source_s2d(data, allowed_source_center_rows=mask)
            self.assertEqual(result['status'], 'failed')
            self.assertIn('allowed_source_center_rows', result['reason'])

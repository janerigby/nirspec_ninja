"""Human slit selection and diagnostic-only S1600 support."""
import tempfile
import unittest
from pathlib import Path
from nirspec_ninja.diagnostics.fs_gradients import read_crowding_notes, GRADIENT_SLITS
from nirspec_ninja.subtract_background.fs_background import USABLE_SLITS
import numpy as np
from nirspec_ninja.diagnostics.fs_gradients import median_good_profiles
from nirspec_ninja.diagnostics.fs_gradients import combined_median_spectrum, spatial_median_spectrum
from nirspec_ninja.diagnostics.fs_gradients import predicted_jbt_background
from nirspec_ninja.diagnostics.fs_gradients import spatial_flatness
from nirspec_ninja.diagnostics.fs_gradients import profile
from nirspec_ninja.diagnostics.fs_gradients import diagnostic_trace
from nirspec_ninja.diagnostics.fs_gradients import trace_row_mask
from nirspec_ninja.diagnostics.fs_gradients import interpolate_background_gaps
from nirspec_ninja.diagnostics.fs_gradients import smooth_background_spline, create_bkg_override


class CrowdingNotesTests(unittest.TestCase):
    def test_override_starter_never_overwrites_existing_choices(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'FS_bkg_override.txt'
            self.assertTrue(create_bkg_override(path, 'SGAS1110-IFU-POS1', ['G235H', 'G140H', 'G235H']))
            marks, _ = read_crowding_notes(path, 'SGAS1110-IFU-POS1', 'G235H')
            self.assertEqual(len(marks), 4)
            self.assertTrue(all(marks.values()))
            edited = path.read_text().replace('1 1 1 1', '0 -1 1 1')
            path.write_text(edited)
            self.assertFalse(create_bkg_override(path, 'other', ['PRISM']))
            self.assertEqual(path.read_text(), edited)

    def test_spline_suppresses_short_scale_structure_and_crosses_gap(self):
        wave = np.linspace(1., 3., 1000)
        truth = .2 + .05 * (wave - 1) + .02 * (wave - 2)**2
        noisy = truth + .03 * np.sin(2 * np.pi * wave / .015)
        noisy[::97] += .3
        noisy[450:550] = np.nan
        noisy[:20] = np.nan; noisy[-20:] = np.nan
        fitted, gaps = smooth_background_spline(wave, noisy, .1)
        self.assertTrue(np.isfinite(fitted[450:550]).all())
        self.assertTrue(gaps[450:550].all())
        self.assertTrue(np.isnan(fitted[:20]).all())
        self.assertLess(np.nanmedian(np.abs(fitted - truth)), .01)
        self.assertLess(np.max(np.abs(np.diff(fitted[440:560]))), .002)

    def read(self, text, target='COSMIC.EYE', grating=None):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'notes.txt'
            path.write_text(text)
            return read_crowding_notes(path, target, grating=grating)

    def test_template_order_case_and_notes(self):
        marks, notes = self.read(
            '# 1 is use, 0 is do not use\n'
            'target S200A2 S200A1 S400A1 S1600A1 Notes\n'
            'Cosmic.eye 1 1 0 1 Galaxy in S400A1\n')
        self.assertEqual(marks, dict(S200A2=True, S200A1=True, S400A1=False, S1600A1=True))
        self.assertEqual(notes, 'Galaxy in S400A1')

    def test_grating_specific_evaluations(self):
        table = ('target grating S200A2 S200A1 S400A1 S1600A1\n'
                 'Cosmic.eye G235M 1 1 0 1\n# second grating\n'
                 'Cosmic.eye prism 1 1 1 0\n')
        medium, notes = self.read(table, grating='g235m')
        prism, _ = self.read(table, grating='PRISM')
        self.assertFalse(medium['S400A1']); self.assertTrue(medium['S1600A1'])
        self.assertTrue(prism['S400A1']); self.assertFalse(prism['S1600A1'])
        self.assertEqual(notes, '')
        self.assertEqual(self.read(table, grating='G140M'), ({}, ''))
        with self.assertRaisesRegex(ValueError, 'duplicate'):
            self.read(table + 'COSMICEYE Prism 0 0 0 0\n', grating='PRISM')

    def test_position_names_use_base_target_with_exact_override(self):
        table = 'target grating S200A1\nSGAS1110 G235M 1\n'
        self.assertTrue(self.read(table, target='SGAS1110-IFU-POS1', grating='G235M')[0]['S200A1'])
        table += 'SGAS1110-IFU-POS1 G235M 0\n'
        self.assertFalse(self.read(table, target='SGAS1110-IFU-POS1', grating='G235M')[0]['S200A1'])
        self.assertTrue(self.read(table, target='SGAS1110-IFU-POS2', grating='G235M')[0]['S200A1'])
        self.assertEqual(self.read(table, target='SGAS11100-IFU-POS1', grating='G235M'), ({}, ''))

    def test_directory_selection_when_product_name_is_base_target(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'notes.txt'
            path.write_text('target grating S200A1\nSGAS1110-IFU-POS1 G140H 1\n'
                            'SGAS1110-IFU-POS2 G140H 0\nSGAS1110 G235M 1\n')
            first = read_crowding_notes(path, 'SGAS1110-IFU-POS1', 'G140H', fallback_target='SGAS1110')
            second = read_crowding_notes(path, 'SGAS1110-IFU-POS2', 'G140H', fallback_target='SGAS1110')
            self.assertTrue(first[0]['S200A1'])
            self.assertFalse(second[0]['S200A1'])
            self.assertTrue(read_crowding_notes(path, 'SGAS1110-IFU-POS1', 'G235M',
                                               fallback_target='SGAS1110')[0]['S200A1'])

    def test_whole_slit_override_is_a_good_mark(self):
        marks, _ = self.read('target grating S200A1\nCosmic.eye PRISM -1\n', grating='PRISM')
        self.assertEqual(marks['S200A1'], -1)
        _, median, names = median_good_profiles([('S200A1', np.array([0., 1.]), np.array([2., 2.]))], marks)
        np.testing.assert_allclose(median, 2.)
        self.assertEqual(names, ['S200A1'])

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

    def test_jbt_date_selection_and_cache(self):
        from types import SimpleNamespace
        from astropy.time import Time
        header = dict(TARG_RA=10., TARG_DEC=20., EXPMID=Time('2024-06-24T12:00:00').mjd)
        model = SimpleNamespace(cache_version='test', bkg_data=dict(
            calendar=np.array([174, 175]), wave_array=np.array([1., 2., 3.]),
            total_bg=np.array([[100., 100., 100.], [1., 2., 3.]])))
        with tempfile.TemporaryDirectory() as tmp:
            wave, values, provenance = predicted_jbt_background(header, tmp, {(10., 20.): model})
            np.testing.assert_array_equal(values, [1., 2., 3.])
            self.assertEqual(provenance['day_zero_based'], 175)
            # Cached predictions work without a model or network retrieval.
            _, cached, _ = predicted_jbt_background(header, tmp, {})
            np.testing.assert_array_equal(cached, values)
            header['EXPMID'] += 1
            with self.assertRaisesRegex(ValueError, 'unavailable'):
                predicted_jbt_background(header, tmp, {(10., 20.): model})

    def test_flatness_detects_edge_excess_not_single_pixel(self):
        x = np.linspace(0, 1, 20)
        flat = np.full(20, 2.)
        self.assertEqual(spatial_flatness(x, flat)['fractional_departure'], 0.)
        isolated = flat.copy(); isolated[8] = 20.
        self.assertEqual(spatial_flatness(x, isolated)['fractional_departure'], 0.)
        edge = flat.copy(); edge[:3] = 3.
        result = spatial_flatness(x, edge)
        self.assertEqual(result['fractional_departure'], .5)
        self.assertLess(result['peak_fraction'], .2)
        self.assertIsNone(spatial_flatness(x, np.zeros(20))['fractional_departure'])

    def test_spatial_gaps_break_profiles_and_candidate(self):
        from unittest.mock import patch
        from nirspec_ninja.subtract_background.fs_background import SlitPixels
        wave = np.broadcast_to(np.arange(10.), (10, 10)).copy()
        brightness = np.ones_like(wave)
        excluded = np.zeros(10, dtype=bool); excluded[4:6] = True
        good = np.broadcast_to(~excluded[:, None], wave.shape).copy()
        slit = SlitPixels('S200A1', brightness, wave, np.ones_like(wave), None)
        with patch('nirspec_ninja.diagnostics.fs_gradients.bg.spatial_pixel_mask',
                   return_value=(good, excluded)):
            x, full, matched, _ = profile(slit)
        self.assertEqual(len(x), 10)
        self.assertTrue(np.isnan(full[4:6]).all())
        self.assertTrue(np.isnan(matched[4:6]).all())
        grid, median, _ = median_good_profiles([('S200A1', x, matched)], {'S200A1': True})
        self.assertTrue(np.isnan(median[(grid > x[3]) & (grid < x[6])]).all())
        self.assertEqual(spatial_flatness(x, matched)['fractional_departure'], 0.)

    def test_diagnostic_trace_finds_synthetic_continuum(self):
        from nirspec_ninja.subtract_background.fs_background import SlitPixels
        wave = np.broadcast_to(np.linspace(1., 5., 120), (32, 120)).copy()
        rows = np.arange(32.)[:, None]
        brightness = np.broadcast_to(2 + 3 * np.exp(-.5 * ((rows - 12) / 1.5)**2), wave.shape).copy()
        slit = SlitPixels('S200A1', brightness, wave, np.full_like(wave, .1), None)
        result = diagnostic_trace(slit)
        self.assertTrue(result['detected'])
        self.assertAlmostEqual(result['center'], 12., delta=1.)
        self.assertGreater(result['fwhm'], 0)
        self.assertNotIn('accepted', result)
        half_width = (result['bounds_flux_rows'][1] - result['center'])
        sigma = result['fwhm'] / np.sqrt(8 * np.log(2))
        from scipy.special import erf
        self.assertAlmostEqual(erf(half_width / (sigma * np.sqrt(2))), .99)
        mask = trace_row_mask(result, 32)
        self.assertTrue(mask[12])
        self.assertFalse(mask[0])

    def test_spectrum_median_excludes_trace_rows(self):
        from unittest.mock import patch
        from nirspec_ninja.subtract_background.fs_background import SlitPixels
        wave = np.broadcast_to(np.arange(10.), (6, 10)).copy()
        brightness = np.full_like(wave, 2.); brightness[:3] = 20.
        slit = SlitPixels('S200A1', brightness, wave, np.ones_like(wave), None)
        extra = np.array([True, True, True, False, False, False])
        with patch('nirspec_ninja.diagnostics.fs_gradients.bg.spatial_pixel_mask',
                   side_effect=lambda _: (np.ones_like(wave, dtype=bool), np.zeros(6, dtype=bool))):
            _, raw = spatial_median_spectrum(slit)
            _, masked = spatial_median_spectrum(slit, extra_excluded_rows=extra)
        np.testing.assert_allclose(raw, 11.)
        np.testing.assert_allclose(masked, 2.)
        np.testing.assert_array_equal(slit.brightness, brightness)

    def test_flatness_uses_trace_masked_rows(self):
        x = np.linspace(0, 1, 20)
        values = np.full(20, 2.); values[6:10] = 4.
        self.assertEqual(spatial_flatness(x, values)['fractional_departure'], 1.)
        mask = trace_row_mask(dict(detected=True, bounds_flux_rows=[5., 10.]), 20)
        masked = np.where(mask, np.nan, values)
        self.assertEqual(spatial_flatness(x, masked)['fractional_departure'], 0.)

    def test_single_remaining_window_is_not_perfect_flatness(self):
        x = np.linspace(0, 1, 10)
        values = np.full(10, np.nan); values[3:6] = [1., 2., 3.]
        self.assertEqual(spatial_flatness(x, values)['status'], 'unavailable')

    def test_background_display_fills_only_interior_gaps(self):
        wave = np.arange(8.)
        values = np.array([np.nan, 1., 2., np.nan, np.nan, 5., 6., np.nan])
        filled, flags = interpolate_background_gaps(wave, values)
        np.testing.assert_allclose(filled, [np.nan, 1., 2., 3., 4., 5., 6., np.nan], equal_nan=True)
        np.testing.assert_array_equal(flags, [False, False, False, True, True, False, False, False])
        self.assertTrue(np.isnan(values[3:5]).all())

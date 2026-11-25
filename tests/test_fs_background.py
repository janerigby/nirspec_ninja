"""Numerical and real-step checks without raw data or CRDS downloads."""
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).parent / '.mplconfig'))
import numpy as np
from astropy import units as u
from jwst import datamodels
from nirspec_ninja.subtract_background import fs_background as bg


def write_s2d(path, slit_levels, start=60000.0, end=60000.001, with_dq=True):
    model = datamodels.MultiSlitModel()
    model.meta.instrument.name = 'NIRSPEC'
    model.meta.instrument.detector = 'NRS1'
    model.meta.instrument.grating = 'G140M'
    model.meta.instrument.filter = 'F100LP'
    model.meta.exposure.type = 'NRS_FIXEDSLIT'
    model.meta.exposure.start_time = start
    model.meta.exposure.end_time = end
    wave = np.linspace(1.0, 2.0, 80)
    for name, level, scatter in slit_levels:
        data = level + (np.arange(32) - 15.5)[:, None] * scatter + np.zeros((32, 80))
        slit = datamodels.SlitModel(data=data, err=np.full_like(data, .1),
                                   wavelength=np.broadcast_to(wave, data.shape).copy(),
                                   dq=np.zeros(data.shape, dtype=np.uint32) if with_dq else None)
        slit.name = name
        slit.meta.bunit_data = 'MJy/sr'
        slit.meta.exposure.type = 'NRS_FIXEDSLIT'
        slit.source_type = 'EXTENDED'
        model.slits.append(slit)
        slit.close()
    model.save(path)
    model.close()
    return wave


class ResamplingTests(unittest.TestCase):
    def test_load_saved_diagnostic_curve_and_missing_configuration(self):
        import json
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            path = tmp / 'G235H_F170LP_smoothed_bkg.csv'
            wave = np.linspace(1., 2., 100)
            np.savetxt(path, np.column_stack((wave, 2 + wave, np.zeros(len(wave)))),
                       delimiter=',', header='wavelength_um,median_bkg_MJy_sr,interpolated', comments='')
            path.with_suffix('.json').write_text(json.dumps(dict(grating='G235H', filter='F170LP',
                units='MJy/sr', contributing_slits=['S200A1'])))
            estimate = bg.load_diagnostic_background(tmp, 'G235H', 'F170LP', np.array([.5, 1.2, 1.8, 2.5]))
            np.testing.assert_allclose(estimate['brightness'], [np.nan, 3.2, 3.8, np.nan], equal_nan=True, atol=.01)
            self.assertTrue(np.isnan(estimate['uncertainty']).all())
            self.assertEqual(estimate['used'][0]['source'], 'diagnostic_smoothed_FS')
            self.assertEqual(len(estimate['used'][0]['sha256']), 64)
            with self.assertRaises(FileNotFoundError):
                bg.load_diagnostic_background(tmp, 'G140H', 'F100LP', wave)

    def test_constant_density_preserved_on_nonuniform_grid(self):
        native = np.array([1, 1.1, 1.3, 1.7, 2.2, 2.8], dtype=float)
        output = np.array([1.2, 1.5, 1.9, 2.3])
        values = bg.resample_surface_brightness(native, np.full(len(native), 7.0), output)
        np.testing.assert_allclose(values, 7.0)

    def test_frequency_integral_conserved(self):
        native = np.array([1, 1.2, 1.5, 1.9, 2.4])
        values = np.array([1, 4, 8, 2, 5])
        output = np.array([1.2, 1.6, 2.0])
        resampled = bg.resample_surface_brightness(native, values, output)
        native_edges = 1 / bg.bin_edges(native)
        out_edges = 1 / bg.bin_edges(output)
        expected_integral = sum(value * max(0, min(a, out_edges[0]) - max(b, out_edges[-1]))
                                for value, a, b in zip(values, native_edges[:-1], native_edges[1:]))
        self.assertAlmostEqual(np.sum(resampled * -np.diff(out_edges)), expected_integral)
        np.testing.assert_allclose(bg.resample_surface_brightness(native, values, native), values)

    def test_masked_pixels_and_missing_ranges_are_not_bridged(self):
        native = np.array([1, 1.2, 1.5, 1.9, 2.4])
        values = np.array([1, 4, np.nan, 2, 5])
        result = bg.resample_surface_brightness(native, values, native)
        self.assertTrue(np.isnan(result[2]))
        np.testing.assert_allclose(result[[0, 1, 3, 4]], values[[0, 1, 3, 4]])
        outside = bg.resample_surface_brightness(native, values, [3.0, 3.5])
        self.assertTrue(np.all(np.isnan(outside)))


class BackgroundTests(unittest.TestCase):
    def test_disabled_screening_uses_all_eligible_slits_without_detection(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / 'jw04125001001_02101_00001_nrs1_s2d.fits'
            wave = write_s2d(path, [('S200A1', 10., .001), ('S200A2', 20., .001),
                                    ('S400A1', 30., .001), ('S200B1', 999., .001)])
            with patch('nirspec_ninja.fixed_slit.source_detection.detect_source_s2d', side_effect=AssertionError('must not detect')):
                estimate = bg.estimate_background([path], None, wave)
            np.testing.assert_allclose(estimate['brightness'], 20.)
            self.assertEqual(len(estimate['used']), 3)
            self.assertTrue(all(not r['screening_enabled'] for r in estimate['used']))

    def test_l3_decisions_apply_by_configuration_to_contemporaneous_spectra(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            files = []
            for index, brightness in enumerate((10., 30.)):
                path = Path(tmp) / f'jw04125001001_02101_{index + 1:05d}_nrs1_s2d.fits'
                wave = write_s2d(path, [('S200A1', brightness, .001), ('S200A2', 999., .001)],
                                 start=60000 + index * 90, end=60000.001 + index * 90)
                files.append(path)
            records = [dict(file='Level3_TARGET_FS_g140m-f100lp_s200a1_s2d.fits',
                            screening_scope='configuration', grating='G140M', filter='F100LP',
                            slit='S200A1', accepted=True),
                       dict(file='Level3_TARGET_FS_g140m-f100lp_s200a2_s2d.fits',
                            screening_scope='configuration', grating='G140M', filter='F100LP',
                            slit='S200A2', accepted=False)]
            for path, level in zip(files, (10., 30.)):
                estimate = bg.estimate_background([path], records, wave)
                np.testing.assert_allclose(estimate['brightness'], level)
                self.assertEqual([r['slit'] for r in estimate['used']], ['S200A1'])
                self.assertEqual(estimate['used'][0]['screening_files'], [records[0]['file']])
            records[0]['filter'] = 'F170LP'
            with self.assertRaisesRegex(ValueError, 'Missing trace decision'):
                bg.estimate_background([files[0]], records, wave)

    def test_l3_screening_accepts_product_names_without_exposure_ids(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / 'Level3_TARGET_FS_g140m-f100lp_s200a1_s2d.fits'
            write_s2d(path, [('S200A1', 2., .001)])
            records = bg.scan_fixed_slit_traces([path], Path(tmp) / 'reports', 'TARGET',
                                              make_plots=False, screening_scope='configuration')
            self.assertEqual(records[0]['exposure'], 'combined L3')
            self.assertEqual(records[0]['screening_scope'], 'configuration')
            self.assertIn('combined L3', (Path(tmp) / 'reports/trace_report.txt').read_text())

    def test_trace_fit_masks_edges_and_never_centers_on_missing_rows(self):
        from nirspec_ninja.fixed_slit import source_detection as source_detect2
        rng = np.random.default_rng(8)
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / 'jw04125001001_02101_00001_nrs1_s2d.fits'
            write_s2d(path, [('S200A2', 2., .001)])
            with datamodels.open(path) as model:
                slit = model.slits[0]
                slit.data += rng.normal(0, .1, slit.data.shape)
                slit.data[:6] += 100  # Edge artifact must not enter the fit.
                slit.data[16] = np.nan  # Nor may a centroid land in an interior gap.
                model.save(path)
            records = bg.scan_fixed_slit_traces([path], Path(tmp) / 'reports', 'TARGET', make_plots=False)
            record = records[0]
            self.assertNotIn(record['center'], record['excluded_rows'])
            self.assertGreaterEqual(record['center'], 6)
            self.assertLess(record['center'], 26)
            pixels = list(bg.read_slit_pixels(path))[0]
            good, excluded = bg.spatial_pixel_mask(pixels)
            result = source_detect2.detect_source_s2d(
                np.where(good & ~excluded[:, None], pixels.brightness, np.nan),
                error_2d=pixels.error, allowed_source_center_rows=~excluded)
            self.assertTrue(np.all(np.isnan(result['profile'][excluded])))
            failure = source_detect2.detect_source_s2d(
                pixels.brightness, error_2d=pixels.error,
                allowed_source_center_rows=np.zeros(pixels.brightness.shape[0], bool))
            self.assertEqual(failure['status'], 'failed')

    def test_approved_estimator_errors_and_ignored_slits(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / 'jw04125001001_02101_00001_nrs1_s2d.fits'
            wave = write_s2d(path, [('S200A1', 10., .01), ('S200A2', 20., .02),
                                    ('S400A1', 30., .03), ('S1600A1', 1000., .1),
                                    ('S200B1', 1000., .1)])
            records = [dict(file=str(path.resolve()), slit=name, accepted=True)
                       for name in bg.USABLE_SLITS]
            estimate = bg.estimate_background([path], records, wave)
            np.testing.assert_allclose(estimate['brightness'], 20.)
            np.testing.assert_array_equal(estimate['n_slits'], 3)
            np.testing.assert_array_equal(estimate['n_rows'], 60)
            np.testing.assert_allclose(estimate['uncertainty'], 1.253 * 1.4826 * .15 / np.sqrt(20), rtol=1e-5)
            records[1]['accepted'] = False
            records[2]['accepted'] = False
            estimate = bg.estimate_background([path], records, wave)
            np.testing.assert_allclose(estimate['brightness'], 10.)
            np.testing.assert_array_equal(estimate['n_slits'], 1)
            saved = bg.save_background_spectrum(estimate, Path(tmp) / 'background')
            with datamodels.open(saved) as model:
                np.testing.assert_array_equal(model.spec[0].spec_table['WAVELENGTH'], wave)
                np.testing.assert_allclose(model.spec[0].spec_table['SB_ERROR'], estimate['uncertainty'])

    def test_s2d_without_dq_uses_weights_for_detection_and_background(self):
        rng = np.random.default_rng(46)
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / 'jw04125001001_02101_00001_nrs1_s2d.fits'
            wave = write_s2d(path, [('S200A1', 2., .001)], with_dq=False)
            with datamodels.open(path) as model:
                slit = model.slits[0]
                self.assertIsNone(slit.dq)
                slit.data += rng.normal(0, .1, slit.data.shape)
                slit.wht = np.ones(slit.data.shape, dtype=np.float32)
                slit.wht[:, :3] = 0
                slit.data[:, :3] = 1e8  # Finite but uncovered pixels must be excluded.
                model.save(path)
            pixels = list(bg.read_slit_pixels(path))[0]
            self.assertIsNone(pixels.dq)
            self.assertTrue(np.all(np.isnan(pixels.brightness[:, :3])))
            self.assertTrue(np.all(np.isnan(pixels.error[:, :3])))
            self.assertTrue(bg.usable_spatial_rows(pixels))
            records = bg.scan_fixed_slit_traces([path], Path(tmp) / 'reports', 'TARGET', make_plots=False)
            self.assertEqual(records[0]['status'], 'ok')
            self.assertFalse(records[0]['dq_available'])
            self.assertTrue(records[0]['accepted'])
            # Use the stored native grid for this mask regression; the S2D
            # schema stores wavelengths at float32 precision. Other tests
            # verify changes of wavelength grid and partial-bin coverage.
            estimate = bg.estimate_background([path], records, pixels.wavelength[16])
            self.assertTrue(np.all(np.isnan(estimate['brightness'][:3])))
            self.assertTrue(np.all(np.isfinite(estimate['brightness'][3:])))
            self.assertLess(np.nanmax(estimate['brightness']), 3)

    def test_exact_exposure_matching_rejects_other_seasons(self):
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / 'jw04125001001_02101_00001_nrs1_s2d.fits'
            write_s2d(path, [('S200A1', 2., .01)])
            ifu = Path(tmp) / 'jw04125001001_02101_00001_nrs1_cal.fits'
            from astropy.io import fits
            fits.PrimaryHDU(header=fits.getheader(path, 0)).writeto(ifu)
            other = Path(tmp) / 'jw04125002001_02101_00001_nrs1_s2d.fits'
            write_s2d(other, [('S200A1', 200., .01)], start=60100., end=60100.001)
            self.assertEqual(bg.contemporaneous_s2d(ifu, [path, other]), [str(path.resolve())])
            with self.assertRaises(FileNotFoundError):
                bg.contemporaneous_s2d(ifu, [other])
            write_s2d(path, [('S200A1', 200., .01)], start=60100., end=60100.001)
            with self.assertRaisesRegex(ValueError, 'do not overlap'):
                bg.contemporaneous_s2d(ifu, [path])

    def test_real_continuum_detector_and_human_report(self):
        rng = np.random.default_rng(46)
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / 'jw04125001001_02101_00001_nrs1_s2d.fits'
            write_s2d(path, [('S200A1', 2., .001), ('S200A2', 2., .001), ('S1600A1', 2., .001)])
            with datamodels.open(path) as model:
                for slit in model.slits:
                    slit.data += rng.normal(0, .1, slit.data.shape)
                    if slit.name == 'S200A2':
                        slit.data += 3 * np.exp(-.5 * ((np.arange(32) - 16) / 1.5)**2)[:, None]
                model.save(path)
            records = bg.scan_fixed_slit_traces([path], Path(tmp) / 'reports', 'TARGET', make_plots=True)
            decisions = {r['slit']: r for r in records}
            self.assertGreater(decisions['S200A2']['corrected_sigma'], 2.5)
            self.assertFalse(decisions['S200A2']['accepted'])
            self.assertTrue(decisions['S200A1']['accepted'])
            self.assertFalse(decisions['S1600A1']['accepted'])
            report = (Path(tmp) / 'reports' / 'trace_report.txt').read_text()
            self.assertIn('TARGET', report)
            self.assertIn('S200A2', report)
            self.assertIn('continuum trace', report)
            self.assertEqual(len(list((Path(tmp) / 'reports').glob('*_trace.png'))), 3)

    def test_nonlinear_cube_wavelength_grid_round_trip(self):
        from astropy.modeling import models
        from gwcs import WCS, coordinate_frames as cf
        wave = np.array([1., 1.12, 1.4, 2.0])
        spatial = cf.Frame2D(name='spatial', axes_order=(0, 1), unit=(u.deg, u.deg))
        spectral = cf.SpectralFrame(name='spectral', axes_order=(2,), unit=u.um)
        output = cf.CompositeFrame([spatial, spectral], name='world')
        detector = cf.CoordinateFrame(3, ('SPATIAL', 'SPATIAL', 'SPECTRAL'),
                                      axes_order=(0, 1, 2), unit=(u.pix, u.pix, u.pix), name='detector')
        transform = (models.Linear1D(slope=.1 * u.deg / u.pix, intercept=0 * u.deg)
                     & models.Linear1D(slope=.1 * u.deg / u.pix, intercept=0 * u.deg)
                     & models.Tabular1D(points=np.arange(4) * u.pix, lookup_table=wave * u.um))
        model = datamodels.IFUCubeModel(data=np.ones((4, 2, 2), dtype=np.float32))
        model.meta.wcs = WCS([(detector, transform), (output, None)])
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            path = Path(tmp) / 'cube_s3d.fits'
            model.save(path)
            np.testing.assert_array_equal(bg.cube_wavelengths(path), wave)
        model.close()

    def test_actual_master_background_subtraction_and_gap_flags(self):
        from jwst.assign_wcs import nirspec
        from types import SimpleNamespace
        wave = np.array([1., 1.5, 2.])
        def world(x, y):
            return x * 0, y * 0, wave[np.asarray(x).astype(int)]
        fake_wcs = SimpleNamespace(bounding_box=((0, 2), (0, 1)))
        class SliceWcs:
            bounding_box = fake_wcs.bounding_box
            __call__ = staticmethod(world)
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            source = Path(tmp) / 'ifu_cal.fits'
            model = datamodels.IFUImageModel(data=np.tile([11., 12., 13.], (2, 1)),
                                            err=np.full((2, 3), .2),
                                            dq=np.zeros((2, 3), dtype=np.uint32))
            model.meta.instrument.name = 'NIRSPEC'
            model.meta.exposure.type = 'NRS_IFU'
            model.meta.target.source_type = 'EXTENDED'
            model.save(source)
            model.close()
            estimate = dict(wavelength=wave, brightness=np.array([1., np.nan, 3.]),
                            uncertainty=np.ones(3), n_slits=np.array([1, 0, 1]), n_rows=np.array([20, 0, 20]))
            master = bg.save_background_spectrum(estimate, Path(tmp) / 'master')
            dest = Path(tmp) / 'subtracted_cal.fits'
            with patch.object(nirspec, 'nrs_ifu_wcs', return_value=[SliceWcs()]):
                flagged = bg.subtract_ifu_background(source, master, dest)
            self.assertEqual(flagged, 2)
            with datamodels.open(dest) as result:
                np.testing.assert_allclose(result.data[:, [0, 2]], 10.)
                self.assertTrue(np.all(np.isnan(result.data[:, 1])))
                self.assertTrue(np.all(result.dq[:, 1] & 1))
                np.testing.assert_allclose(result.err, .2)
                self.assertEqual(result.meta.cal_step.master_background, 'COMPLETE')
            with datamodels.open(source) as original:
                np.testing.assert_allclose(original.data, [[11, 12, 13]] * 2)


class CubeWorkflowTests(unittest.TestCase):
    def test_two_seasons_subtracted_independently_before_combination(self):
        import json
        from jwst.assign_wcs import nirspec
        from jwst.pipeline import Spec3Pipeline
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            tmp = Path(tmp)
            files, s2d_files, records = [], [], []
            wave = np.linspace(1., 2., 80)
            for index, (date, sky) in enumerate(((60000., 10.), (60100., 200.)), 1):
                prefix = f'jw0412500{index}001_02101_00001_nrs1'
                s2d = tmp / (prefix + '_s2d.fits')
                write_s2d(s2d, [('S200A1', sky, .01)], start=date, end=date+.001)
                s2d_files.append(str(s2d))
                records.append(dict(file=str(s2d.resolve()), slit='S200A1', accepted=True))
                source = tmp / (prefix + '_cal.fits')
                model = datamodels.IFUImageModel(data=np.full((2, 80), sky+5),
                                                err=np.full((2, 80), .2),
                                                dq=np.zeros((2, 80), np.uint32))
                model.meta.instrument.name = 'NIRSPEC'
                model.meta.instrument.grating = 'G140M'
                model.meta.instrument.filter = 'F100LP'
                model.meta.exposure.type = 'NRS_IFU'
                model.meta.exposure.start_time = date
                model.meta.exposure.end_time = date + .001
                model.meta.target.source_type = 'EXTENDED'
                model.save(source)
                model.close()
                files.append(source)
            reference = tmp / 'TARGET_NOBG_g140m-f100lp_s3d.fits'
            model = datamodels.IFUCubeModel(data=np.ones((80, 2, 2)))
            model.meta.instrument.grating = 'G140M'
            model.meta.instrument.filter = 'F100LP'
            model.save(reference)
            model.close()
            association = tmp / 'input_asn.json'
            association.write_text(json.dumps({'products': [{'name': 'TARGET_NOBG', 'members': [
                dict(expname=str(filename), exptype='science') for filename in files]}]}))
            observed = []
            class SliceWcs:
                bounding_box = ((0, 79), (0, 1))
                def __call__(self, x, y):
                    return x*0, y*0, wave[np.asarray(x).astype(int)]
            class CubeBuilder:
                def __init__(self, **parameters):
                    self.parameters = parameters
                    self.assertions = []
                def run(self, filename):
                    members = json.loads(Path(filename).read_text())['products'][0]['members']
                    for member in members:
                        with datamodels.open(member['expname']) as corrected:
                            observed.append((Path(member['expname']).name, corrected.data.copy(),
                                             corrected.meta.cal_step.master_background))
                    cube = datamodels.IFUCubeModel(data=np.full((80, 2, 2), 5.))
                    cube.meta.instrument.grating = 'G140M'
                    cube.meta.instrument.filter = 'F100LP'
                    cube.save(tmp / 'TARGET_FSBGSUB_g140m-f100lp_s3d.fits')
                    cube.close()
            with patch.object(bg, 'cube_wavelengths', return_value=wave), \
                 patch.object(nirspec, 'nrs_ifu_wcs', return_value=[SliceWcs()]), \
                 patch('jwst.pipeline.Spec3Pipeline', CubeBuilder):
                result = bg.build_subtracted_cube(
                    reference, association, {'steps': {'master_background': {'skip': False}}},
                    s2d_files, records, tmp / 'reports', 'TARGET')
            self.assertTrue(Path(result).is_file())
            self.assertEqual(len(observed), 2)
            for filename, values, status in observed:
                self.assertIn(filename, [p.name for p in files])
                np.testing.assert_allclose(values, 5., atol=1e-5)
                self.assertEqual(status, 'COMPLETE')
            manifest = json.loads((tmp / 'reports' / 'background_report.json').read_text())
            self.assertEqual([r['start_mjd'] for r in manifest], [60000., 60100.])
            for r in manifest:
                self.assertEqual(len(r['fs_files']), 1)
                self.assertEqual(bg.exposure_key(r['fs_files'][0]), r['exposure'])
                self.assertEqual(r['status'], 'subtracted')
            # Saved diagnostic mode must subtract the chosen CSV, not recompute
            # a background from the very different contemporaneous FS levels.
            diagnostics = tmp / 'diagnostics'; diagnostics.mkdir()
            saved = diagnostics / 'G140M_F100LP_smoothed_bkg.csv'
            np.savetxt(saved, np.column_stack((wave, np.full(len(wave), 2.), np.zeros(len(wave)))),
                       delimiter=',', header='wavelength_um,median_bkg_MJy_sr,interpolated', comments='')
            saved.with_suffix('.json').write_text(json.dumps(dict(grating='G140M', filter='F100LP',
                units='MJy/sr', contributing_slits=['S200A1'])))
            observed.clear()
            with patch.object(bg, 'cube_wavelengths', return_value=wave), \
                 patch.object(nirspec, 'nrs_ifu_wcs', return_value=[SliceWcs()]), \
                 patch('jwst.pipeline.Spec3Pipeline', CubeBuilder), \
                 patch.object(bg, 'estimate_background', side_effect=AssertionError('raw background used')):
                bg.build_subtracted_cube(reference, association, {'steps': {'master_background': {'skip': False}}},
                                         [], None, tmp / 'reports', 'TARGET',
                                         diagnostic_spectra_dir=diagnostics)
            for (_, values, status), sky in zip(observed, [10., 200.]):
                np.testing.assert_allclose(values, sky + 3., atol=1e-5)
                self.assertEqual(status, 'COMPLETE')
            manifest = json.loads((tmp / 'reports' / 'background_report.json').read_text())
            self.assertTrue(all(r['background_source'] == 'diagnostic_smoothed_FS' for r in manifest))
            self.assertFalse(list(tmp.glob('fsbg_*')))
            with datamodels.open(files[0]) as original:
                np.testing.assert_allclose(original.data, 15.)

    def test_offset_subset_retains_alignment_and_exposure_order(self):
        import asdf
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            source = Path(tmp) / 'offsets.asdf'
            with asdf.AsdfFile(dict(filename=['a_cal.fits', 'b_cal.fits', 'c_cal.fits'],
                                   raoffset=[1., 2., 3.], decoffset=[4., 5., 6.], units='arcsec')) as f:
                f.write_to(source)
            output = bg.subset_cube_offsets(source, ['/temp/c_cal.fits', '/temp/a_cal.fits'],
                                            Path(tmp) / 'subset.asdf')
            with asdf.open(output) as f:
                self.assertEqual(f['filename'], ['c_cal.fits', 'a_cal.fits'])
                self.assertEqual(f['raoffset'], [3., 1.])
                self.assertEqual(f['decoffset'], [6., 4.])
                self.assertEqual(f['units'], 'arcsec')
            with self.assertRaisesRegex(ValueError, 'Expected one cube offset'):
                bg.subset_cube_offsets(source, ['missing_cal.fits'], Path(tmp) / 'missing.asdf')

    def test_real_spec3_settings_can_be_cloned(self):
        from copy import deepcopy
        from jwst.pipeline import Spec3Pipeline
        original = Spec3Pipeline()
        original.cube_build.linear_wave = False
        original.cube_build.scalexy = .05
        original.extract_1d.skip = True
        original.cube_build.override_cubepar = '/test/jwst_nirspec_cubepar_0009.fits'
        parameters = bg.clone_spec3_parameters(original)
        parameters['steps']['master_background']['skip'] = True
        parameters['output_file'] = 'TARGET_FSBGSUB'
        cloned = Spec3Pipeline(**parameters)
        self.assertFalse(cloned.cube_build.linear_wave)
        self.assertEqual(cloned.cube_build.scalexy, .05)
        self.assertTrue(cloned.extract_1d.skip)
        self.assertTrue(cloned.master_background.skip)
        self.assertEqual(cloned.cube_build.get_ref_override('cubepar'), '/test/jwst_nirspec_cubepar_0009.fits')

    def test_withhold_policy_reports_insufficient_coverage(self):
        import json
        from astropy.io import fits
        with tempfile.TemporaryDirectory(dir=Path(__file__).parent) as tmp:
            tmp = Path(tmp)
            source = tmp / 'jw04125001001_02101_00001_nrs1_cal.fits'
            fs = tmp / 'jw04125001001_02101_00001_nrs1_s2d.fits'
            wave = write_s2d(fs, [('S200A1', 2., .01)])
            fits.PrimaryHDU(header=fits.getheader(fs, 0)).writeto(source)
            reference = tmp / 'NOBG_s3d.fits'
            fits.PrimaryHDU(header=fits.getheader(fs, 0)).writeto(reference)
            association = tmp / 'asn.json'
            association.write_text(json.dumps({'products': [{'name': 'NOBG', 'members': [
                dict(expname=str(source), exptype='science')]}]}))
            records = [dict(file=str(fs.resolve()), slit='S200A1', accepted=False)]
            with patch.object(bg, 'cube_wavelengths', return_value=wave), \
                 patch('jwst.pipeline.Spec3Pipeline') as pipeline:
                result = bg.build_subtracted_cube(reference, association, {}, [fs], records,
                                                  tmp / 'reports', 'TARGET', gap_policy='withhold')
            self.assertIsNone(result)
            pipeline.assert_not_called()
            self.assertIn('withheld', (tmp / 'reports' / 'background_report.txt').read_text())


if __name__ == '__main__':
    unittest.main()

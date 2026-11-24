"""Check migrated notebook workers and shared process mechanics offline."""
import multiprocessing
import tempfile
import unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch

from nirspec_ninja.pipeline import notebook_workers as workers
from nirspec_ninja.pipeline import multiprocessing as logistics


class PipelineMigrationTests(unittest.TestCase):
    def test_workers_are_importable_by_spawned_processes(self):
        # Imports in a fresh interpreter matter for macOS and notebook workers.
        with multiprocessing.get_context('spawn').Pool(1) as pool:
            self.assertEqual(pool.starmap(workers.simple_task_to_multiprocess,
                                         [(3, None, None)]), [4])

    def test_det1_still_skips_cleaning_without_mutating_parameters(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp) / 'example_uncal.fits'
            out = Path(tmp) / 'example_rate.fits'
            params = {'jump': {'expand_large_events': True}}
            pipeline = SimpleNamespace(save_model=lambda *a, **k: None,
                                       run=lambda _: out.touch())
            with patch('jwst.pipeline.calwebb_detector1.Detector1Pipeline.build_config',
                       return_value=({}, None)) as config, \
                 patch('jwst.pipeline.calwebb_detector1.Detector1Pipeline.from_config_section',
                       return_value=pipeline):
                self.assertEqual(workers.run_jwst_det1(source, tmp, params), str(out.resolve()))
            self.assertTrue(config.call_args.kwargs['steps']['clean_flicker_noise']['skip'])
            self.assertFalse(config.call_args.kwargs['steps']['ramp_fit']['skip'])
            self.assertTrue(config.call_args.kwargs['save_results'])
            self.assertEqual(params, {'jump': {'expand_large_events': True}})

    def test_det1_rejects_conflicting_stage_requirements_before_running(self):
        cases = [({'clean_flicker_noise': {'skip': False}}, 'clean_flicker_noise.skip=True'),
                 ({'ramp_fit': {'skip': True}}, 'ramp_fit')]
        for params, message in cases:
            with self.subTest(params=params), \
                 tempfile.TemporaryDirectory() as tmp, \
                 patch('jwst.pipeline.calwebb_detector1.Detector1Pipeline.build_config') as config:
                self.assertIsNone(workers.run_jwst_det1(
                    Path(tmp) / 'example_uncal.fits', tmp, params))
                config.assert_not_called()
                self.assertIn(message, (Path(tmp) / 'example_uncal_pipecrash.txt').read_text())

    def test_spec2_retains_fft_and_original_mode_source(self):
        import numpy as np
        from jwst import datamodels
        with tempfile.TemporaryDirectory() as tmp, datamodels.ImageModel(
                data=np.ones((8, 8), dtype=np.float32)) as model:
            source = Path(tmp) / 'example_rate.fits'
            model.meta.exposure.type = 'NRS_MSASPEC'
            model.meta.instrument.opmode = 'MSASPEC'
            model.meta.cal_step.clean_flicker_noise = 'SKIPPED'
            model.save(source)
            original = source.read_bytes()
            mask = Path(tmp) / 'custom_masks/example_mask.fits'
            mask.parent.mkdir()
            mask.touch()
            params = {'clean_flicker_noise': {'fit_method': 'fft', 'skip': False,
                      'save_results': True, 'save_noise': True, 'save_background': True}}
            seen = []
            def run(rate):
                seen.append((rate.meta.exposure.type, rate.meta.instrument.opmode))
                return None
            pipeline = SimpleNamespace(clean_flicker_noise=SimpleNamespace(run=run), run=run)
            with patch('jwst.pipeline.Spec2Pipeline.build_config', return_value=({}, None)) as config, \
                 patch('jwst.pipeline.Spec2Pipeline.from_config_section', return_value=pipeline), \
                 patch('nirspec_ninja.pipeline.patches.nirspec_flicker_patch.install_patch', return_value={'max_weight': 4.0, 'runtime_verified': True}) as install, \
                 patch('nirspec_ninja.diagnostics.flicker_noise_diagnostics.write_noise_report') as report:
                self.assertTrue(workers.run_jwst_spec2(source, tmp, params, save_results=False))
            report.assert_called_once_with(source.resolve(), source.with_name('example_rate_after_cleaning.fits').resolve(), mask.resolve(), report_file=source.parent.resolve() / 'flicker_noise_rms.txt', fit_method='fft')
            install.assert_called_once_with(max_weight=4.0)
            self.assertTrue(pipeline.clean_flicker_noise.flicker_patch_status['runtime_verified'])
            self.assertTrue(mask.with_suffix('.flicker_patch.json').exists())
            cleaning = config.call_args.kwargs['steps']['clean_flicker_noise']
            self.assertEqual(cleaning['fit_method'], 'fft')
            self.assertFalse(cleaning['skip'])
            self.assertEqual(cleaning['user_mask'], str(mask.resolve()))
            self.assertTrue(cleaning['save_results'])
            self.assertTrue(cleaning['save_noise'])
            self.assertTrue(cleaning['save_background'])
            self.assertFalse(config.call_args.kwargs['save_results'])
            self.assertEqual(seen, [('NRS_FIXEDSLIT', 'FIXEDSLIT')])
            self.assertEqual(source.read_bytes(), original)
            self.assertNotIn('user_mask', params['clean_flicker_noise'])

    def test_spec2_missing_opmode_preserved_in_cleaned_diagnostic(self):
        import numpy as np
        from astropy.io import fits
        from jwst import datamodels
        with tempfile.TemporaryDirectory() as tmp, datamodels.ImageModel(
                data=np.ones((8, 8), dtype=np.float32)) as model:
            source = Path(tmp) / 'example_rate.fits'
            model.meta.exposure.type = 'NRS_MSASPEC'
            model.meta.cal_step.clean_flicker_noise = 'SKIPPED'
            model.save(source)
            original = source.read_bytes()
            with datamodels.open(source) as saved:
                self.assertNotIn('opmode', saved.meta.instrument.instance)
                with self.assertRaises(AttributeError):
                    _ = saved.meta.instrument.opmode
            mask = Path(tmp) / 'custom_masks/example_mask.fits'
            mask.parent.mkdir()
            fits.HDUList([fits.PrimaryHDU(), fits.ImageHDU(np.ones((8, 8)), name='SCI')]).writeto(mask)
            seen = []

            def clean(rate):
                self.assertEqual(rate.meta.exposure.type, 'NRS_FIXEDSLIT')
                self.assertEqual(rate.meta.instrument.opmode, 'FIXEDSLIT')
                result = rate.copy()
                result.data *= 2
                result.meta.cal_step.clean_flicker_noise = 'COMPLETE'
                return result

            step = SimpleNamespace(run=clean)

            def run(rate):
                result = step.run(rate)
                seen.append((result.meta.exposure.type, result.meta.instrument.opmode))
                return result

            pipeline = SimpleNamespace(clean_flicker_noise=step, run=run)
            with patch('jwst.pipeline.Spec2Pipeline.build_config', return_value=({}, None)) as config, \
                 patch('jwst.pipeline.Spec2Pipeline.from_config_section', return_value=pipeline), \
                 patch('nirspec_ninja.pipeline.patches.nirspec_flicker_patch.install_patch', return_value={'max_weight': 4.0, 'runtime_verified': True}) as install:
                self.assertTrue(workers.run_jwst_spec2(source, tmp,
                    {'clean_flicker_noise': {'fit_method': 'fft'}}))
            cleaning = config.call_args.kwargs['steps']['clean_flicker_noise']
            for name in ('skip', 'save_results', 'save_noise', 'save_background'):
                self.assertFalse(cleaning[name])
            self.assertTrue(config.call_args.kwargs['save_results'])
            report = (Path(tmp) / 'flicker_noise_rms.txt').read_text()
            self.assertIn('example_rate_after_cleaning.fits', report)
            self.assertIn('FFT', report)
            self.assertIn('MAD_%', report)
            # A contradictory choice must fail explicitly, rather than be overwritten.
            with patch('jwst.pipeline.Spec2Pipeline.build_config') as rejected:
                self.assertFalse(workers.run_jwst_spec2(source, tmp,
                    {'clean_flicker_noise': {'skip': True}}))
                rejected.assert_not_called()
            self.assertIn('clean_flicker_noise.skip=False',
                          (Path(tmp) / 'example_rate_spec2_pipecrash.txt').read_text())
            self.assertEqual(source.read_bytes(), original)
            self.assertEqual(seen, [('NRS_FIXEDSLIT', 'FIXEDSLIT')])
            diagnostic = Path(tmp) / 'example_rate_after_cleaning.fits'
            with datamodels.open(diagnostic) as saved:
                self.assertEqual(saved.meta.exposure.type, 'NRS_MSASPEC')
                self.assertNotIn('opmode', saved.meta.instrument.instance)
                self.assertTrue(np.all(saved.data == 2))
            with fits.open(diagnostic) as hdus:
                self.assertNotIn('OPMODE', hdus[0].header)

    def test_crds_preflight_finishes_before_workers_and_reuses_context(self):
        with patch.object(logistics, '_crds_preflight_context', None), \
             patch.dict(logistics.os.environ, {}, clear=False), \
             patch('crds.get_default_context', return_value='jwst_test.pmap') as context, \
             patch('crds.client.api.dump_mappings') as mappings, \
             patch.object(logistics.subprocess, 'run') as command:
            self.assertEqual(logistics.preflight_crds(['example_rate.fits'], 'Spec2'),
                             'jwst_test.pmap')
            self.assertEqual(logistics.os.environ['CRDS_CONTEXT'], 'jwst_test.pmap')
            logistics.preflight_crds([], 'Detector1')
            context.assert_called_once_with('jwst')
            mappings.assert_called_once_with('jwst_test.pmap')
            self.assertTrue(command.call_args.kwargs['check'])
            self.assertIn('--sync-references', command.call_args.args[0])

    def test_scheduler_preserves_sequential_results(self):
        self.assertEqual(logistics.run_pipeline_tasks(
            workers.simple_task_to_multiprocess, [(1, None, None), (5, None, None)], 1), [2, 6])
        self.assertEqual(logistics.run_pipeline_tasks(
            workers.simple_task_to_multiprocess, [(1, None, None)], 2), [2])


if __name__ == '__main__':
    unittest.main()

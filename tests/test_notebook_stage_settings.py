"""Offline checks for default behavior, overrides, and notebook task arguments."""
import inspect
import json
import tempfile
import unittest
from copy import deepcopy
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import MagicMock, Mock, patch

import numpy as np
from jwst import datamodels
from nirspec_ninja.pipeline import notebook_workers as workers

NOTEBOOK = Path(__file__).resolve().parents[1] / 'Reduce-and-extract_fixedslit-multicores_v2.ipynb'


def stage_source(marker):
    notebook = json.loads(NOTEBOOK.read_text())
    return next(''.join(c['source']) for c in notebook['cells']
                if marker in ''.join(c['source']))


def stage_settings(marker):
    namespace = dict(RUN_MASKS=True, RUN_LEVEL3=True, RUN_EXTRACT=True,
                     dirs_to_handle=[], Pool=MagicMock(),
                     run_jwst_spec3=workers.run_jwst_spec3, run_extract=workers.run_extract, summarize_extract=Mock())
    exec(stage_source(marker), namespace)
    return namespace


class NotebookStageSettingsTests(unittest.TestCase):
    def test_mask_defaults_and_overrides_preserve_input_and_output_handling(self):
        explicit = stage_settings('if RUN_MASKS')['mask_param_nested_dict']
        cases = [None, {}, explicit,
                 {'mask_FS_by_wcs': False, 'mask_illumination_nsigma': 7.0}]
        with tempfile.TemporaryDirectory() as tmp, datamodels.ImageModel(
                data=np.ones((8, 8), dtype=np.float32)) as model:
            source = Path(tmp) / 'example_rate.fits'
            model.meta.cal_step.clean_flicker_noise = 'SKIPPED'
            model.save(source)
            original = source.read_bytes()
            default_kwargs = None
            for params in cases:
                with self.subTest(params=params):
                    before = deepcopy(params)
                    with patch('nirspec_ninja.masks.custom_nirspec_mask_2.make_custom_NIRSpec_mask',
                               side_effect=lambda src, dest, **kw: Path(dest).touch()) as builder:
                        if params is None:
                            result = workers.make_rate_custom_mask(source, tmp)
                        else:
                            result = workers.make_rate_custom_mask(source, tmp, params)
                    self.assertEqual(result, str((Path(tmp) / 'custom_masks/example_mask.fits').resolve()))
                    kwargs = builder.call_args.kwargs
                    self.assertTrue(kwargs['overwrite'])
                    self.assertEqual(kwargs['input_dir'], str(Path(tmp).resolve()))
                    self.assertEqual(source.read_bytes(), original)
                    self.assertEqual(params, before)
                    if params is None:
                        default_kwargs = kwargs
                        self.assertEqual(kwargs['mask_illumination_nsigma'], 5.0)
                    elif not params or params == explicit:
                        self.assertEqual(kwargs, default_kwargs)
                    else:
                        self.assertFalse(kwargs['mask_FS_by_wcs'])
                        self.assertEqual(kwargs['mask_illumination_nsigma'], 7.0)
                        self.assertTrue(kwargs['mask_grating_extrapolate'])

    def test_spec3_defaults_and_overrides_still_save_and_close(self):
        explicit = stage_settings('if RUN_LEVEL3')['spec3_param_nested_dict']
        for params in (None, {}, explicit, {'outlier_detection': {'skip': True}}):
            with self.subTest(params=params):
                before = deepcopy(params)
                result = SimpleNamespace(close=Mock())
                pipeline = SimpleNamespace(outlier_detection=SimpleNamespace(skip=None),
                                           run=Mock(return_value=result))
                with patch('jwst.pipeline.Spec3Pipeline', return_value=pipeline):
                    if params is None:
                        returned = workers.run_jwst_spec3('association.json', 'output')
                    else:
                        returned = workers.run_jwst_spec3('association.json', 'output', params)
                self.assertIsNone(returned)
                self.assertEqual(pipeline.outlier_detection.skip,
                                 bool(params and params.get('outlier_detection', {}).get('skip')))
                self.assertTrue(pipeline.save_results)
                self.assertEqual(pipeline.output_dir, 'output')
                pipeline.run.assert_called_once_with('association.json')
                result.close.assert_called_once_with()
                self.assertEqual(params, before)

    def test_extraction_clip_choice_and_existing_failure_handling(self):
        explicit = stage_settings('if RUN_EXTRACT')['clip_nan']
        with patch('nirspec_ninja.fixed_slit.extraction.wrap_median_combine_level3_nirspecFS') as extract:
            extract.return_value = [{'status': 'ok'}]
            self.assertEqual(workers.run_extract('input', 'output')['status'], 'ok')
            extract.assert_called_with('input', 'output', clip_nan=5, return_report=True, detect_traces=True)
            self.assertEqual(workers.run_extract('input', 'output', explicit)['status'], 'ok')
            extract.assert_called_with('input', 'output', clip_nan=5, return_report=True, detect_traces=True)
            self.assertEqual(workers.run_extract('input', 'output', 3)['status'], 'ok')
            extract.assert_called_with('input', 'output', clip_nan=3, return_report=True, detect_traces=True)
            extract.side_effect = RuntimeError('example extraction failure')
            with patch('builtins.print') as report:
                self.assertEqual(workers.run_extract('input', 'output', 3)['status'], 'failed')
            self.assertIn('example extraction failure', str(report.call_args_list))

    def test_notebook_pool_jobs_pass_visible_settings_to_correct_workers(self):
        captured = []

        class Pool:
            def __init__(self, *args, **kwargs):
                pass

            def __enter__(self):
                return self

            def __exit__(self, *args):
                pass

            def starmap(self, worker, tasks):
                tasks = list(tasks)
                for task in tasks:
                    inspect.signature(worker).bind(*task)
                captured.append((worker, tasks))
                return ['success' for _ in tasks]

        def glob(pattern):
            if pattern.endswith('_rate.fits'):
                return ['/out/L2a/example_rate.fits']
            if pattern.endswith('cal.fits'):
                return ['/out/L2b/example_cal.fits']
            return []

        namespace = dict(RUN_MASKS=True, RUN_LEVEL3=True, RUN_EXTRACT=True,
                         dirs_to_handle=['/input/'], leggos_out='/out/', pipev='3.0.0',
                         sort_out_paths=lambda *a, **k: ('epoch', 'target', '/input/', '/out/'),
                         folder_L2a='L2a/', folder_L2b='L2b/', folder_L3='L3/', folder_extracted='extracted/',
                         glob=SimpleNamespace(glob=glob), Path=Path,
                         os=SimpleNamespace(path=SimpleNamespace(exists=lambda _: True)),
                         writel3asn=Mock(), download_metadata=Mock(), Pool=Pool,
                         make_rate_custom_mask=workers.make_rate_custom_mask,
                         run_jwst_spec3=workers.run_jwst_spec3, run_extract=workers.run_extract, summarize_extract=Mock())
        for marker in ('if RUN_MASKS', 'if RUN_LEVEL3', 'if RUN_EXTRACT'):
            exec(stage_source(marker), namespace)
        self.assertEqual([w for w, _ in captured],
                         [workers.make_rate_custom_mask, workers.run_jwst_spec3, workers.run_extract])
        self.assertEqual(captured[0][1][0][2], namespace['mask_param_nested_dict'])
        self.assertEqual(captured[1][1][0][2], namespace['spec3_param_nested_dict'])
        self.assertEqual(captured[2][1][0][2], namespace['clip_nan'])


if __name__ == '__main__':
    unittest.main()

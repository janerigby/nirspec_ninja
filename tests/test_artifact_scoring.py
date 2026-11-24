"""Column scoring uses measured amplitudes without any clean reference."""
import csv, json, tempfile, unittest
from pathlib import Path
from types import SimpleNamespace
from unittest.mock import patch
import numpy as np
from nirspec_ninja.artifacts import detect

class ArtifactScoringTests(unittest.TestCase):
    def test_cleaned_rate_discovery_and_mask_identity(self):
        exposure = 'jw02127002001_03101_00005_nrs1'
        with tempfile.TemporaryDirectory() as tmp:
            nested = Path(tmp)/'nested'; nested.mkdir()
            products = [nested/(exposure+suffix) for suffix in
                        ('_rate.fits', '_rate_after_cleaning.fits', '_clean_flicker_noise.fits')]
            for product in products:
                product.touch()
                self.assertEqual(detect.exposure_stem(product), exposure)
                self.assertEqual(detect.mask_name(product), exposure+'_mask.fits')
            (nested/(exposure+'_mask.fits')).touch()
            self.assertEqual(detect.inputs([tmp]), sorted(products))
            self.assertEqual(detect.inputs([str(products[1])]), [products[1]])

    def test_scores_without_reference_or_baseline(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'synthetic_rate.fits'; source.touch()
            dest = Path(tmp)/'results'
            args = SimpleNamespace(inputs=[str(source)], output=str(dest), mask_dir=None,
                spectrum_percentile=90, min_glow_rows=30, glow_threshold=None,
                column_threshold=None, PNG=False, save_maps=False, detailed_csv=True,
                candidate_level=5, file_list=None)
            columns = np.zeros(2048); columns[100] = 120
            with patch.object(detect, 'metadata', return_value=('IFU','COMPLETE','DN/s')), patch.object(detect, 'extract', return_value=('IFU',2.,{'columns':columns})), patch.object(detect, 'measure', return_value=(np.zeros((256,256)),{'glow_status':'ok'},None)), patch.object(detect.np, 'load', side_effect=AssertionError('must not load a reference')), patch.object(detect.subprocess, 'run', side_effect=AssertionError('must not compare obsolete scores')):
                self.assertEqual(detect.run(args), 0)
            row = next(csv.DictReader((dest/'scores.csv').open()))
            self.assertEqual(float(row['column_score']), 120.)
            record = json.loads((dest/'run.json').read_text())
            self.assertNotIn('reference',record)
            self.assertEqual(record['version'],'0.3')
            self.assertFalse((dest/'score_comparison.png').exists())

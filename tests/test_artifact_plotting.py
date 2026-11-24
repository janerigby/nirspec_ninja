"""Plotting must not silently mix obsolete reference scores with Input."""
import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch
from nirspec_ninja.artifacts import plot_score_comparison as plotting

class PlottingTests(unittest.TestCase):
    def test_cli_does_not_add_reference_samples(self):
        with patch('sys.argv', ['plot_score_comparison', 'catscores.csv', '--show-histograms', '--hide-mystery']), patch.object(plotting, 'plot_scores') as plot:
            plotting.main()
        self.assertEqual(plot.call_args.args[3], [])
        self.assertFalse(plot.call_args.args[4])
        self.assertTrue(plot.call_args.args[6])
    def test_only_input_is_plotted_by_default(self):
        with tempfile.TemporaryDirectory() as tmp:
            source = Path(tmp)/'scores.csv'
            with source.open('w',newline='') as handle:
                writer = csv.DictWriter(handle,fieldnames=['file','mode','status','glow_status','column_score','glow_score'])
                writer.writeheader();writer.writerow(dict(file='/data/IFU_bad/exposure.fits',mode='IFU',status='ok',glow_status='ok',column_score=120,glow_score=1))
            from matplotlib.axes import Axes
            original = Axes.scatter
            labels=[]
            def capture(self,*args,**kwargs):
                labels.append(kwargs.get('label',''))
                return original(self,*args,**kwargs)
            with patch.object(Axes,'scatter',capture):
                output=plotting.plot_scores(source,show_histograms=True,show_mystery=False)
            self.assertTrue(output.exists())
            self.assertTrue(any('Input' in label for label in labels))
            self.assertFalse(any('Visual' in label for label in labels))

"""Check the shared writer's science/background association behavior."""
import json
import tempfile
import unittest
from pathlib import Path

from nirspec_ninja.pipeline.associations import writel3asn


class AssociationTests(unittest.TestCase):
    def test_science_order_product_name_and_empty_backgrounds(self):
        for backgrounds in (None, []):
            with self.subTest(backgrounds=backgrounds), tempfile.TemporaryDirectory() as tmp:
                path = Path(tmp) / 'association.json'
                self.assertIsNone(writel3asn(['second_cal.fits', 'first_cal.fits'],
                                            backgrounds, str(path), 'target_FS_g395m-f290lp'))
                result = json.loads(path.read_text())
                self.assertEqual(result['products'], [{
                    'name': 'target_FS_g395m-f290lp',
                    'members': [{'expname': 'second_cal.fits', 'exptype': 'science'},
                                {'expname': 'first_cal.fits', 'exptype': 'science'}]}])
                self.assertEqual(result['asn_type'], 'None')

    def test_background_order_and_existing_output_replacement(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / 'association.json'
            path.write_text('previous association')
            writel3asn(['science_cal.fits'], ['bg2_x1d.fits', 'bg1_x1d.fits'],
                       str(path), 'target')
            result = json.loads(path.read_text())
            self.assertEqual(result['products'][0]['members'], [
                {'expname': 'science_cal.fits', 'exptype': 'science'},
                {'expname': 'bg2_x1d.fits', 'exptype': 'background'},
                {'expname': 'bg1_x1d.fits', 'exptype': 'background'}])


if __name__ == '__main__':
    unittest.main()

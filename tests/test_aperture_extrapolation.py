"""Coverage and failure behavior of detector-edge aperture continuation."""
import unittest
import warnings
import numpy as np
from nirspec_ninja.masks.aperture_extrapolation import extend_aperture

class ExtrapolationTests(unittest.TestCase):
    def test_curvature_full_aperture_and_dark_gaps(self):
        footprint = np.zeros((70, 160),bool)
        # Wide aperture includes all rows, not just a central bright spectrum.
        for x in range(160):
            y = int(round(8+.0008*x*x))
            footprint[y:y+8,x] = True
        mask = np.ones((150,300),bool)
        mask[140,0] = False  # independent exclusion must survive
        region = np.zeros_like(mask)
        report = extend_aperture(mask,region,footprint,40,20)
        self.assertLess(report['withheld_error_pixels'],2)
        for col in range(200,300):
            y=int(round(28+.0008*(col-40)**2))
            self.assertFalse(mask[y:y+8,col].any())
        self.assertTrue(mask[120:135,250].all())
        self.assertFalse(mask[140,0])
        self.assertTrue(region[:,299].any())

    def test_short_fit_warns_and_widens(self):
        footprint=np.ones((6,12),bool)
        mask=np.ones((60,100),bool)
        with warnings.catch_warnings(record=True) as captured:
            warnings.simplefilter('always')
            report=extend_aperture(mask,None,footprint,30,20)
        self.assertTrue(report['warning'])
        self.assertTrue(captured)
        self.assertFalse(mask[18:28,99].all())
        self.assertTrue(mask[40:,99].all())

    def test_empty_has_no_effect(self):
        mask=np.ones((20,30),bool)
        self.assertIsNone(extend_aperture(mask,None,np.zeros((3,5),bool),2,3))
        self.assertTrue(mask.all())

    def test_wcs_endpoint_islands_do_not_drive_extrapolation(self):
        footprint=np.zeros((60,160),bool)
        for x in range(160):
            y=int(round(12-.03*x))
            footprint[y:y+8,x]=True
        footprint[:5,:12]=True
        footprint[40:45,-12:]=True
        mask=np.ones((90,300),bool)
        report=extend_aperture(mask,None,footprint,40,10)
        self.assertGreater(report['rejected_width_columns'],0)
        self.assertFalse(report['warning'])
        self.assertTrue(mask[45:,299].all())

if __name__=='__main__': unittest.main()

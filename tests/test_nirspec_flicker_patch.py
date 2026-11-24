"""Behavioral regression checks for process-local flicker patch."""
import sys,unittest,tempfile
from pathlib import Path
import numpy as np
from astropy.io import fits
from stdatamodels.jwst.datamodels import dqflags
from jwst.clean_flicker_noise import nsclean,clean_flicker_noise as cfn
from nirspec_ninja.pipeline.patches import nirspec_flicker_patch as patch
from nirspec_ninja.masks.custom_nirspec_mask_2 import refine_mask_from_rate,refine_mask_file_from_rate

class FlickerPatchTests(unittest.TestCase):
    def setUp(self):
        self.bindings=nsclean.NSClean,cfn.NSClean,cfn._mask_unusable
    def tearDown(self):
        nsclean.NSClean,cfn.NSClean,cfn._mask_unusable=self.bindings
    def test_caps_and_both_kernel_orientations(self):
        for cap in (2,4,8):
            status=patch.install_patch(cap)
            self.assertTrue(status['runtime_verified'])
            self.assertEqual(status['max_weight'],cap)
    def test_invalid_caps_fail_before_binding_changes(self):
        for cap in (0,.5,-1,np.nan,np.inf,True):
            with self.assertRaises(ValueError):patch.install_patch(cap)
            self.assertIs(cfn.NSClean,self.bindings[1])
    def test_saturated_samples_have_no_influence_but_early_samples_survive(self):
        patch.install_patch(4)
        mask=np.ones((128,128),bool)
        dq=np.zeros((128,128),np.uint8);dq[30:60,60]=dqflags.group['SATURATED']
        early=mask.copy();cfn._mask_unusable(early,np.zeros_like(dq))
        late=mask.copy();cfn._mask_unusable(late,dq)
        self.assertTrue(early[40,60]);self.assertFalse(late[40,60])
        rng=np.random.default_rng(22);a=rng.normal(size=(128,128)).astype('float32')
        b=a.copy();b[30:60,60]=50000
        za=cfn.NSClean('NRS1',late.copy()).clean(a.copy(),buff=True)
        zb=cfn.NSClean('NRS1',late.copy()).clean(b.copy(),buff=True)
        self.assertTrue(np.allclose((a-za)[late],(b-zb)[late],atol=1e-5,rtol=0))
    def test_calibrated_pixel_flags_excluded_and_other_flags_retained(self):
        rng=np.random.default_rng(33);data=rng.normal(1,.01,size=(32,32))
        dq=np.zeros(data.shape,np.uint32);mask=np.ones(data.shape,bool);mask[10,10]=False
        for col,name in enumerate(patch.PIXEL_EXCLUSIONS):dq[4,col]=dqflags.pixel[name]
        dq[4,4]=dqflags.pixel['JUMP_DET'];dq[4,5]=dqflags.pixel['SATURATED']
        result=refine_mask_from_rate(mask,data,nsigma=20,pixel_dq=dq)
        self.assertFalse(result[4,:4].any());self.assertTrue(result[4,4:6].all())
        self.assertFalse(result[10,10]);self.assertTrue(mask[4,:4].all())
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'mask.fits'
            fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(mask.astype('float32'),name='SCI')]).writeto(path)
            saved=refine_mask_file_from_rate(path,data,nsigma=20,pixel_dq=dq)
            self.assertTrue(np.array_equal(saved,fits.getdata(path,'SCI').astype(bool)))
    def test_bad_dq_shape_rejected(self):
        with self.assertRaises(ValueError):
            refine_mask_from_rate(np.ones((8,8),bool),np.ones((8,8)),pixel_dq=np.zeros((9,9),np.uint32))

if __name__=='__main__':unittest.main()

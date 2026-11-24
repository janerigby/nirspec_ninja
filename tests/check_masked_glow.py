"""Checks physical scope requirements, not detection accuracy on real data."""
from pathlib import Path
from tempfile import TemporaryDirectory
import numpy as np
from astropy.io import fits
from nirspec_ninja.artifacts.masked_glow import measure

with TemporaryDirectory(dir=Path(__file__).parent) as td:
 root=Path(td); rate=root/'synthetic_rate.fits'
 mask=np.ones((2048,2048),dtype=np.float32);mask[:,400:512]=0
 fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(mask,name='SCI')]).writeto(root/'synthetic_mask.fits')
 rng=np.random.default_rng(42)
 base=rng.normal(0,1,mask.shape).astype(np.float32);base[:,400:512]=10
 cases=[('no_glow',0,0,0),('29_rows',29,25,0),('30_rows',30,25,1),('64_rows',64,25,1),('fainter_than_spectra',64,5,0),('spectra_only',64,25,0)]
 for name,height,amplitude,expected in cases:
  a=base.copy()
  x=440 if name=='spectra_only' else 1000
  a[64:64+height,x:x+64]+=amplitude
  fits.HDUList([fits.PrimaryHDU(),fits.ImageHDU(a,name='SCI'),fits.ImageHDU(np.zeros(mask.shape,dtype=np.uint32),name='DQ')]).writeto(rate,overwrite=True)
  scores,info,_=measure(rate,root,1)
  assert (scores.max()>0)==bool(expected),(name,info)
  if name=='30_rows': assert info['glow_component_rows']==30,info
  print(name,'passed',flush=True)
print('All masked glow behavior checks passed.')

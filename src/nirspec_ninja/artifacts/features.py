import numpy as np
from scipy.ndimage import median_filter, gaussian_filter
from astropy.io import fits

def features(p, compute_glow=True):
 with fits.open(p,memmap=False) as h:
  a=h['SCI'].data.astype(float); dq=h['DQ'].data
  mode='IFU' if 'IFU' in h[0].header['TEMPLATE'] else 'MSA'
  a[((dq&1)!=0) | ~np.isfinite(a)]=np.nan
  a[:4]=np.nan;a[-4:]=np.nan;a[:,:4]=np.nan;a[:,-4:]=np.nan
  noise=np.nanmedian(np.abs(a[:,1:]-a[:,:-1]))/0.953872
  # Lower quantiles over row bands suppress bright horizontal spectral traces.
  bands=np.nanpercentile(a.reshape(16,128,2048),20,axis=1)
  base=median_filter(np.nan_to_num(bands,nan=np.nanmedian(bands)),size=(1,101))
  residual=np.nan_to_num(bands-base, nan=0.0)
  # Fourth largest same-sign band requires structure in several separated bands.
  pos=np.sort(residual,axis=0)[-4]; neg=-np.sort(residual,axis=0)[3]
  col=np.maximum(pos,neg)/max(noise,1e-10)
  if not compute_glow:
   return mode,float(noise),col,None
  # Median tiles suppress point sources and thin spectra.
  b=np.nanmedian(a.reshape(256,8,256,8),axis=(1,3))
  b=np.nan_to_num(b,nan=np.nanmedian(b))
  # Remove globally coherent columns, then smooth in the vertical direction.
  b-=np.nan_to_num(np.nanmedian(a,axis=0), nan=np.nanmedian(a)).reshape(256,8).mean(axis=1)[None,:]
  broad=median_filter(b,size=(7,1))
  broad-=np.median(broad,axis=1)[:,None]
  g=gaussian_filter(broad,2)/max(noise,1e-10)
  return mode,float(noise),col,g

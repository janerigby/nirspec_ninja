"""Extended positive glow measured only in supplied background geometry."""
from pathlib import Path
import warnings
import numpy as np
from astropy.io import fits
from scipy.ndimage import label

INPUT_SUFFIXES=('_rate.fits','_rate_after_cleaning.fits','_clean_flicker_noise.fits')

def exposure_stem(path):
 """Use the same exposure identity for rate and full-frame cleaned products."""
 name=Path(path).name
 for suffix in INPUT_SUFFIXES:
  if name.endswith(suffix): return name[:-len(suffix)]
 raise ValueError(f'Unsupported input filename: {name}; expected one of {INPUT_SUFFIXES}')

def mask_name(path):
 return exposure_stem(path)+'_mask.fits'

def measure(rate, mask_dir, noise, spectrum_percentile=90, min_rows=30, image=None):
 name=mask_name(rate)
 if isinstance(mask_dir,dict):
  matches=mask_dir.get(name,[])
  if len(matches)>1:
   matches=sorted(matches)
   keys=['TEMPLATE','DETECTOR','GRATING','FILTER','SUBSTRT1','SUBSTRT2','SUBSIZE1','SUBSIZE2','FASTAXIS','SLOWAXIS','OBS_ID','EXPSTART']
   with fits.open(matches[0],memmap=False) as first:
    for other in matches[1:]:
     with fits.open(other,memmap=False) as second:
      if not np.array_equal(first['SCI'].data,second['SCI'].data) or any(first[0].header.get(k)!=second[0].header.get(k) for k in keys):
       raise ValueError(f'Conflicting masks for basename: {name}')
  if not matches: raise FileNotFoundError(f'No matching illumination mask: {name}')
  path=matches[0]
 else:
  path=Path(mask_dir)/name
 if not path.is_file():
  raise FileNotFoundError(f'No matching illumination mask: {path}')
 if image is None:
  with fits.open(rate,memmap=False) as h:
   a=h['SCI'].data.astype(float); dq=h['DQ'].data; rh=h[0].header
 else:
  data,dq,rh=image
  a=np.asarray(data,dtype=float)
 valid=np.isfinite(a)&((dq&1)==0)
 with fits.open(path,memmap=False) as h:
  mask=h['SCI'].data; mh=h[0].header
  if mask.shape!=a.shape or not np.all(np.isin(mask,[0,1])):
   raise ValueError('Mask must match SCI dimensions and contain only 0/1')
  for key in ['TEMPLATE','DETECTOR','GRATING','FILTER','SUBSTRT1','SUBSTRT2','SUBSIZE1','SUBSIZE2','FASTAXIS','SLOWAXIS','OBS_ID','EXPSTART']:
   if rh.get(key)!=mh.get(key): raise ValueError(f'Mask metadata mismatch: {key}')
 valid[:4]=False;valid[-4:]=False;valid[:,:4]=False;valid[:,-4:]=False
 bg=valid&(mask==1); sky=valid&(mask==0)
 if not bg.any() or not sky.any(): raise ValueError('Mask needs valid background and spectral-footprint pixels')
 baseline=float(np.percentile(a[bg],20))
 brightness=max(float(np.percentile(a[sky],spectrum_percentile))-baseline,noise)
 # Only background contributes: no interpolation/smoothing across spectra.
 with warnings.catch_warnings():
  warnings.simplefilter('ignore',RuntimeWarning)
  tiles=np.nanmedian(np.where(bg,a,np.nan).reshape(256,8,256,8),axis=(1,3))
 coverage=bg.reshape(256,8,256,8).sum(axis=(1,3))
 excess=tiles-baseline
 eligible=(coverage>=16)&np.isfinite(excess)&(excess>=brightness)
 components,n=label(eligible) # edge-connected background tiles
 scores=np.zeros((256,256)); max_span=0; count=0
 for k in range(1,n+1):
  yy,xx=np.where(components==k)
  if (yy.max()-yy.min()+1)*8<min_rows: continue
  # Count actual detector rows containing bright background pixels in this patch.
  occupied=np.zeros(2048,dtype=bool)
  for ty,tx in zip(yy,xx):
   sl=np.s_[ty*8:(ty+1)*8,tx*8:(tx+1)*8]
   occupied[ty*8:(ty+1)*8] |= np.any(bg[sl]&(a[sl]-baseline>=brightness),axis=1)
  rows=int(occupied.sum())
  if rows<min_rows: continue
  count+=1;max_span=max(max_span,rows)
  scores[yy,xx]=excess[yy,xx]/noise
 return scores,dict(mask_file=str(path.resolve()),glow_status='ok',glow_background_baseline=baseline,glow_spectrum_brightness=brightness,glow_min_rows=min_rows,glow_component_count=count,glow_component_rows=max_span),mask

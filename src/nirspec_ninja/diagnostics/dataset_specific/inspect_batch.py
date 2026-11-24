"""Dataset-specific diagnostic; original inputs/output paths retained."""
if __name__ == "__main__":
    import os
    from pathlib import Path
    O=Path('/Users/jrrigby1/AI_experiments/Coding/Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/horizontal_stripes');os.environ['MPLCONFIGDIR']=str(O/'mplconfig')
    from astropy.io import fits
    import numpy as np,json
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    P=Path('/Users/jrrigby1/AI_experiments/Temp/batch2');PREV=Path('/Users/jrrigby1/AI_experiments/Temp/leggos1');CONTROL=Path('/Users/jrrigby1/AI_experiments/Coding/JWST_pipeline/diagnostics/leggos_nsclean')
    for f in P.glob('*.fits'):
     with fits.open(f) as h:
      print(f.name,[(x.name,None if x.data is None else x.data.shape) for x in h if x.name!='ASDF'])
      print({k:h[0].header.get(k) for k in ['DATE','CAL_VER','CRDS_CTX','EXP_TYPE','DETECTOR','NGROUPS','S_CLNFNS']})
    for detector in ['nrs1','nrs2']:
     stem='jw04125018001_02101_00001_'+detector
     with fits.open(P/(stem+'_clean_flicker_noise.fits')) as h:
      ramp=h['SCI'].data.astype('float32');t=h[0].header['TGROUP']
      proxy=(ramp[0,-1]-ramp[0,0])/((ramp.shape[1]-1)*t)
      diff=(ramp[0,1]-ramp[0,0])/t
     rate=fits.getdata(P/(stem+'_rate.fits'))
     fig,axs=plt.subplots(1,3,figsize=(16,6),layout='constrained')
     for ax,data,title in zip(axs,[ramp[0,0],diff,rate],['Cleaned ramp, group 0','Cleaned group difference 1−0 / TGROUP','Final RATE']):
      lo,hi=np.nanpercentile(data,[5,95]) if title.startswith('Cleaned ramp') else (-.02,.02)
      ax.imshow(data,origin='lower',cmap='RdBu_r',vmin=lo,vmax=hi);ax.set_title(detector+' '+title,fontsize=10);ax.set_xlabel('Detector x');ax.set_ylabel('Detector y')
     fig.savefig(O/(detector+'_overview.png'),dpi=150)
     if detector=='nrs1':
      previous=fits.getdata(PREV/(stem+'_rate.fits'));candidate=fits.getdata(CONTROL/'det1_new_clip3_both_rate.fits');mask=fits.getdata(P/(stem+'_mask.fits'));expected=fits.getdata(CONTROL/'new_clip3_mask.fits')
      print('NRS1 mask difference vs tested candidate',np.count_nonzero(mask!=expected),'rate difference vs candidate',np.nanstd(rate-candidate),np.nanmax(np.abs(rate-candidate)))
      fig,axs=plt.subplots(1,3,figsize=(16,6),layout='constrained')
      for ax,data,title in zip(axs,[rate-previous,rate-candidate,mask],['New − prior RATE','New − controlled candidate RATE','New mask, white = fit']):
       ax.imshow(data,origin='lower',cmap='RdBu_r' if not title.startswith('New mask') else 'gray',vmin=0 if title.startswith('New mask') else -.01,vmax=1 if title.startswith('New mask') else .01);ax.set_title(title);ax.set_xlabel('Detector x')
      fig.savefig(O/'nrs1_comparisons.png',dpi=150)

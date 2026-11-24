"""Dataset-specific diagnostic; original inputs/output paths retained."""
if __name__ == "__main__":
    import os
    from pathlib import Path
    O=Path('/Users/jrrigby1/AI_experiments/Coding/Custom_masks_for_NIRSpec_CleanFlickerNoise/diagnostics/horizontal_stripes');os.environ['MPLCONFIGDIR']=str(O/'mplconfig')
    import numpy as np,json
    from astropy.io import fits
    from scipy.ndimage import median_filter
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    P=Path('/Users/jrrigby1/AI_experiments/Temp/batch2');C=Path('/Users/jrrigby1/AI_experiments/Coding/JWST_pipeline/diagnostics/leggos_nsclean');OLD=Path('/Users/jrrigby1/AI_experiments/Temp/leggos1');s='jw04125018001_02101_00001_nrs1'
    with fits.open(P/(s+'_clean_flicker_noise.fits')) as h, fits.open(C/'preclean_ramp.fits') as f, fits.open(OLD/(s+'_clean_flicker_noise.fits')) as old:
     t=h[0].header['TGROUP'];metrics={}
     fig,axs=plt.subplots(2,3,figsize=(17,10),layout='constrained')
     for row,g in enumerate([1,14]):
      before=f['SCI'].data[0,g];after=h['SCI'].data[0,g];prior=old['SCI'].data[0,g]
      base=f['SCI'].data[0,0];duration=g*t
      arrays=[(before-base)/duration,(after-base)/duration,(before-after)/duration]
      for ax,d,title in zip(axs[row],arrays,['Uncleaned accumulated slope proxy','Cleaned accumulated slope proxy','Subtracted correction / elapsed time']):
       ax.imshow(d,origin='lower',cmap='RdBu_r',vmin=-.01,vmax=.01);ax.set_title(f'Group {g}: '+title,fontsize=10);ax.set_xlabel('Detector x');ax.set_ylabel('Detector y')
      metrics[str(g)]={'correction_rms':float(np.nanstd(arrays[2][10:-10,10:1400])),'new_vs_prior_ramp_rms_DN':float(np.nanstd((after-prior)[10:-10,10:1400]))}
      fits.writeto(O/f'group{g}_subtracted_correction_dn_s.fits',arrays[2],overwrite=True)
     fig.savefig(O/'nrs1_ramp_corrections.png',dpi=150)
     fig,axs=plt.subplots(1,3,figsize=(17,6),layout='constrained')
     for ax,d,title in zip(axs,[arrays[1],arrays[2],(after-prior)/duration],['New cleaned, group14','Subtracted correction, group14','New − previous cleaned, group14']):
      smoothed=median_filter(np.nan_to_num(d),size=(1,21));ax.imshow(smoothed,origin='lower',cmap='RdBu_r',vmin=-.005,vmax=.005);ax.set_title(title+'\n21-pixel horizontal median for display',fontsize=10);ax.set_xlabel('Detector x')
     fig.savefig(O/'nrs1_horizontal_structure.png',dpi=150)
     (O/'ramp_metrics.json').write_text(json.dumps(metrics,indent=2));print(json.dumps(metrics,indent=2))

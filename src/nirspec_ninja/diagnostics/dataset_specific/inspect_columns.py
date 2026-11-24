"""Dataset-specific diagnostic; original inputs/output paths retained."""
if __name__ == "__main__":
    from pathlib import Path
    import warnings
    import numpy as np
    from astropy.io import fits
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from jwst.clean_flicker_noise.nsclean import NSClean
    warnings.filterwarnings('ignore', message='All-NaN slice encountered')
    root=Path('/Users/jrrigby1/AI_experiments/Coding/JWST_pipeline')
    fig,axs=plt.subplots(3,3,figsize=(16,12),layout='constrained')
    for r,p in enumerate(sorted((root/'L2a').glob('*_rate.fits'))):
     a=fits.getdata(p,'SCI'); dq=fits.getdata(p,'DQ')
     m=fits.getdata(next((root/'L2a/custom_masks').rglob(p.name.replace('_rate','_mask'))),'SCI')!=0
     good=m&np.isfinite(a)&((dq&1)==0)
     v=np.where(good,a,np.nan)
     axs[r,0].imshow(a[:,800:1400],origin='lower',aspect='auto',extent=(801,1400,.5,2048.5),cmap='RdBu_r',vmin=-.012,vmax=.008)
     axs[r,0].set_title(p.stem+' (DN/s)');axs[r,0].axvspan(970,1265,color='yellow',alpha=.1)
     axs[r,1].imshow(m[:,800:1400],origin='lower',aspect='auto',extent=(801,1400,.5,2048.5),cmap='gray')
     axs[r,1].set_title('Mask: white=usable')
     for ylo,yhi in [(4,200),(850,950),(1110,1200),(1800,2044)]:
      axs[r,2].plot(np.arange(1,2049),np.nanmedian(v[ylo:yhi],axis=0),label=f'rows {ylo+1}–{yhi}',lw=.65)
     axs[r,2].set_xlim(800,1400);axs[r,2].set_ylim(-.02,.015);axs[r,2].legend(fontsize=8);axs[r,2].set_title('Background column medians by row band')
     for ax in axs[r]: ax.set_xlabel('DS9 Image column (1-based)')
     # Condition numbers of the exact weighted Fourier matrix used by NSClean,
     # before group-dependent clipping and DQ exclusions.
     cleaner=NSClean('NRS1',m)
     cond=[]
     for x in range(2048):
      mask=cleaner.mask[x]; ys=np.flatnonzero(mask).reshape(-1,1)
      basis=np.exp(2*np.pi*1j*ys*np.arange(cleaner.nvec)/2048)/len(ys)
      A=(basis*cleaner.p_matrix[x,mask,None]).astype(np.complex64)
      cond.append(np.linalg.cond(A.conj().T@A))
     cond=np.array(cond)
     np.savez(root/'diagnostics'/f'{p.stem}_mask_condition.npz',normal_matrix_condition=cond,usable_pixels=m.sum(axis=0))
     print(p.name,'nvec',cleaner.nvec,'condition median/max all',np.median(cond),cond.max(),'condition median/max 970:1265',np.median(cond[969:1265]),cond[969:1265].max(),'worst col',cond.argmax()+1,flush=True)
     for lo,hi in [(1,600),(970,1265),(1400,1800)]:
      print('cols',lo,hi,'usable count median/min',np.median(m[:,lo-1:hi].sum(axis=0)),m[:,lo-1:hi].sum(axis=0).min(),flush=True)
    fig.savefig(root/'diagnostics/affected_columns.png',dpi=140)

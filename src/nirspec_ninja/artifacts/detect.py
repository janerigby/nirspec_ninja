#!/usr/bin/env python3
"""Independent-exposure PRISM/NRS1 artifact candidates; scores are not probabilities."""
import argparse, csv, json, os, sys, warnings, subprocess
from pathlib import Path
os.environ.setdefault('MPLCONFIGDIR', str(Path(__file__).resolve().parent/'mpl_config'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import maximum_filter
from astropy.io import fits
from nirspec_ninja.artifacts.features import features
from nirspec_ninja.artifacts.masked_glow import measure, exposure_stem, mask_name, INPUT_SUFFIXES

VERSION='0.3'

def inputs(paths):
 out=set()
 for name in paths:
  p=Path(name)
  if p.is_dir():
   for suffix in INPUT_SUFFIXES: out.update(p.rglob('*'+suffix))
  elif p.is_file(): out.add(p)
  else: raise ValueError(f'Input does not exist: {p}')
 if not out: raise ValueError('No FITS files found')
 return sorted(out)

def metadata(p):
 h=fits.getheader(p)
 t=h.get('TEMPLATE','')
 if 'IFU' in t: mode='IFU'
 elif 'MultiObject' in t: mode='MSA'
 else: raise ValueError(f'Unsupported TEMPLATE: {t!r}')
 if h.get('GRATING')!='PRISM' or h.get('DETECTOR')!='NRS1':
  raise ValueError('Only PRISM/NRS1 is supported')
 with fits.open(p) as hd:
  if hd['SCI'].data.shape!=(2048,2048): raise ValueError('Only full-frame 2048x2048 files are supported')
 return mode,h.get('S_CLNFNS','MISSING'),h.get('BUNIT',fits.getheader(p,'SCI').get('BUNIT',''))

def extract(p, compute_glow=True):
 metadata(p)
 with warnings.catch_warnings():
  warnings.filterwarnings('ignore',message='All-NaN slice encountered')
  m,n,c,g=features(p,compute_glow=compute_glow)
 if not np.isfinite(n) or n<=0: raise ValueError('Cannot estimate a positive finite noise scale')
 if not np.all(np.isfinite(c)) or (g is not None and not np.all(np.isfinite(g))): raise ValueError('Non-finite feature map')
 return m,n,dict(columns=c,glow=g)

def fit(args):
 raise ValueError('Clean-reference fitting is obsolete. Use score; no reference is required.')


def peaks(g,n=5):
 local=g==maximum_filter(g,size=13)
 iy,ix=np.where(local)
 order=np.argsort(g[iy,ix])[::-1][:n]
 return [(int(ix[i]*8+4),int(iy[i]*8+4),float(g[iy[i],ix[i]])) for i in order]

def diagnostic(p,dest,c,g,noise,row,level,illumination=None,data=None):
 if data is None:
  with fits.open(p,memmap=False) as h: a=h['SCI'].data.astype(float)
 else: a=np.asarray(data,dtype=float)
 finite=a[np.isfinite(a)]
 lo,hi=np.percentile(finite,[5,98])
 fig,axs=plt.subplots(2,2,figsize=(14,12))
 extent=(.5,2048.5,.5,2048.5)
 axs[0,0].imshow(a,origin='lower',extent=extent,cmap='gray',vmin=lo,vmax=hi)
 axs[0,0].set_title('SCI: 5–98 percentile stretch')
 axs[0,1].imshow(np.arcsinh((a-np.median(finite))/noise),origin='lower',extent=extent,cmap='RdBu_r',vmin=-5,vmax=5)
 axs[0,1].set_title('SCI: signed asinh stretch in noise units')
 if illumination is not None:
  axs[0,0].contour(np.arange(2048)+1,np.arange(2048)+1,illumination,levels=[.5],colors=['lime'],linewidths=.4)
  axs[0,0].set_title('SCI; green boundaries mark excluded spectral footprints')
 for x in np.where(c>level)[0]: axs[0,1].axvline(x+1,color='cyan',alpha=.25,lw=.6)
 if np.min(g)<level<np.max(g):
  axs[0,1].contour(np.arange(256)*8+4.5,np.arange(256)*8+4.5,g,levels=[level],colors=['orange'],linewidths=.7)
 im=axs[1,0].imshow(np.log10(1+g),origin='lower',extent=extent,cmap='magma',vmin=0,vmax=max(1,np.log10(1+np.max(g))))
 fig.colorbar(im,ax=axs[1,0],label='log10(1 + broad-structure score)')
 axs[1,0].set_title('Extended bright background glow' if row.get('glow_status')=='ok' else 'Glow unavailable: matching mask missing')
 for x,y,s in (peaks(g,3) if np.max(g)>0 else []):
  axs[1,0].plot(x+.5,y+.5,'+',color='cyan');axs[1,0].text(x+.5,y+.5,f' ({x+.5:g},{y+.5:g})',color='cyan',fontsize=8)
 axs[1,1].plot(np.arange(2048)+1,c,color='teal',lw=.8)
 axs[1,1].set_yscale('symlog',linthresh=1)
 axs[1,1].axhline(level,color='gray',ls='--',label=f'Candidate display level {level:g}')
 axs[1,1].set_title('Column coherence across ≥4 of 16 row bands')
 axs[1,1].set_xlabel('Column x (1-based)');axs[1,1].set_ylabel('Column score');axs[1,1].legend(fontsize=8)
 for ax in [axs[0,0],axs[0,1],axs[1,0]]:
  ax.set_xlabel('x (1-based FITS/DS9)');ax.set_ylabel('y (1-based FITS/DS9)')
 fig.suptitle(f'{p.name}\n{row["mode"]} | glow={row.get("glow_score", "")}, columns={row["column_score"]:.3g} | S_CLNFNS={row["clean_flicker_noise"]}',fontsize=12)
 fig.tight_layout();fig.savefig(dest,dpi=125);plt.close(fig)

def run(args):
 dest=Path(args.output);dest.mkdir(parents=True,exist_ok=True)
 mask_index={}
 if args.mask_dir:
  for mask_path in Path(args.mask_dir).rglob('*_mask.fits'):
   mask_index.setdefault(mask_path.name,[]).append(mask_path)
 rows=[]; errors=0
 for i,p in enumerate(inputs(args.inputs),1):
  row={'file':str(p.resolve()),'status':'ok','error':''}
  try:
   mode,history,unit=metadata(p);m,noise,maps=extract(p,compute_glow=False)
   c=maps['columns']
   g=np.zeros((256,256));illumination=None
   try:
    selected_index=mask_index
    if not args.mask_dir:
     name=mask_name(p)
     exposure=exposure_stem(p)
     mask_roots=[p.parent/'custom_masks']
     if p.parent.name=='L2b_pipeclean':
      mask_roots.append(p.parent.parent/'L2a'/'custom_masks')
     candidates=[]
     for mask_root in mask_roots:
      for candidate in mask_root.rglob('*.fits'):
       if candidate.name==name or (candidate.name=='mask.fits' and (candidate.parent.name==exposure or candidate.parent.name.startswith(exposure+'_'))):
        candidates.append(candidate)
     selected_index={name:candidates}
    g,glow_info,illumination=measure(p,selected_index,noise,args.spectrum_percentile,args.min_glow_rows)
   except FileNotFoundError as exc:
    glow_info=dict(glow_status='missing_mask',mask_file='',glow_error=str(exc))
    print(f'WARNING {p.name}: {exc}. Supply the exposure illumination mask via --mask-dir, adjacent custom_masks, or sibling L2a/custom_masks for L2b_pipeclean inputs; glow is unavailable.',file=sys.stderr,flush=True)
   row.update(glow_info)
   row.update(mode=mode,clean_flicker_noise=history,unit=unit,noise_scale=noise,column_score=float(np.max(c)),glow_score=float(np.max(g)),column_peak_x=int(np.argmax(c)+1),glow_peak_x=float(np.argmax(g)%256*8+4.5),glow_peak_y=float(np.argmax(g)//256*8+4.5))
   if not np.max(g):
    row.update(glow_peak_x='',glow_peak_y='')
   if row['glow_status']!='ok':
    row.update(glow_score='',glow_peak_x='',glow_peak_y='')
   row['decision']='unthresholded' if row['glow_status']=='ok' else 'review'
   if args.glow_threshold is not None or args.column_threshold is not None:
    row['decision']='reject' if (args.glow_threshold is not None and row['glow_status']=='ok' and row['glow_score']>=args.glow_threshold) or (args.column_threshold is not None and row['column_score']>=args.column_threshold) else ('keep' if row['glow_status']=='ok' else 'review')
   # Sequence prefix prevents collisions between files with the same basename.
   row['diagnostic']=''
   if args.PNG:
    png=dest/f'{i:04d}_{p.stem}.png';row['diagnostic']=str(png.resolve())
    diagnostic(p,png,c,g,noise,row,args.candidate_level,illumination)
   if args.save_maps: np.savez_compressed(dest/f'{i:04d}_{p.stem}_maps.npz',column_score=c,glow_score=g,glow_available=np.array(row['glow_status']=='ok'))
   print(f'[{i}] {p.name}: glow={row.get("glow_score", "")}, columns={row["column_score"]:.3g}, S_CLNFNS={history}',flush=True)
  except Exception as e:
   errors+=1;row.update(status='error',error=str(e),decision='review');print(f'ERROR {p}: {e}',file=sys.stderr,flush=True)
  rows.append(row)
 fields=['file','mode','column_score','glow_score','status','glow_status','decision','error']
 if args.detailed_csv:
  with (dest/'scores_detailed.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in rows for k in row)));w.writeheader();w.writerows(rows)
 with (dest/'scores.csv').open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,extrasaction='ignore');w.writeheader();w.writerows(rows)
 (dest/'run.json').write_text(json.dumps({'version':VERSION,'column_method':'noise-normalized coherent column amplitude; no clean reference','inputs':args.inputs,'glow_threshold':args.glow_threshold,'column_threshold':args.column_threshold,'candidate_display_level':args.candidate_level,'PNG':args.PNG,'files':len(rows),'errors':errors,'mask_dir':str(Path(args.mask_dir).resolve()) if args.mask_dir else None,'file_list':str(Path(args.file_list).resolve()) if args.file_list else None,'mask_lookup':'explicit recursive basename' if args.mask_dir else 'adjacent custom_masks plus sibling L2a/custom_masks for L2b_pipeclean inputs','min_glow_rows':args.min_glow_rows,'spectrum_percentile':args.spectrum_percentile},indent=2))
 return 1 if errors else 0

def boolean_value(value):
 if value.lower()=='true': return True
 if value.lower()=='false': return False
 raise argparse.ArgumentTypeError('Expected True or False')

def main():
 parser=argparse.ArgumentParser(description=__doc__); sub=parser.add_subparsers(dest='command',required=True)
 f=sub.add_parser('fit',help='Obsolete; clean-reference fitting is no longer used');f.add_argument('--clean',nargs='+',required=True);f.add_argument('--output',required=True)
 r=sub.add_parser('score',help='Score individual files or recursively scan directories');r.add_argument('inputs',nargs='*');r.add_argument('--reference',help='Obsolete; ignored. No clean reference is used');r.add_argument('--output',required=True);r.add_argument('--glow-threshold',type=float);r.add_argument('--column-threshold',type=float);r.add_argument('--candidate-level',type=float,default=5,help='Diagnostic overlay level only; does not decide keep/reject');r.add_argument('--save-maps',action='store_true')
 r.add_argument('--mask-dir',help='Optional global directory of named *_mask.fits; otherwise search adjacent custom_masks and sibling L2a/custom_masks for L2b_pipeclean inputs')
 r.add_argument('--detailed-csv',action='store_true',help='Also write scores_detailed.csv with all diagnostic metadata')
 r.add_argument('--file-list',help='Text file with one absolute RATE or clean_flicker_noise FITS path per line; blank lines and # comments ignored')
 r.add_argument('--plot-baseline',help='Obsolete; ignored. Old reference-based scores are not compared')
 r.add_argument('--PNG',type=boolean_value,default=True,metavar='{True,False}',help='Generate diagnostic PNGs (default True); False writes CSV without PNGs')
 r.add_argument('--min-glow-rows',type=int,default=30)
 r.add_argument('--spectrum-percentile',type=float,default=90)
 args=parser.parse_args()
 if args.command=='score':
  if args.file_list:
   try:
    listed=[line.strip() for line in Path(args.file_list).read_text().splitlines() if line.strip() and not line.lstrip().startswith('#')]
   except OSError as exc: parser.error(str(exc))
   if any(not Path(line).is_absolute() for line in listed): parser.error('--file-list requires absolute paths, one per line')
   args.inputs.extend(listed)
  if not args.inputs: parser.error('Supply RATE/clean_flicker_noise paths/directories or --file-list')
  if args.reference or args.plot_baseline: print('WARNING: --reference and --plot-baseline are obsolete and ignored.',file=sys.stderr)
 if args.command=='score' and (args.min_glow_rows<1 or not 0<args.spectrum_percentile<100): parser.error('Require min-glow-rows > 0 and 0 < spectrum-percentile < 100')
 try:
  if args.command=='fit': fit(args);return 0
  return run(args)
 except Exception as e: print(f'ERROR: {e}',file=sys.stderr);return 1
if __name__=='__main__': sys.exit(main())

"""Required helpers moved from supplied jrr; existing names retained."""
import numpy as np
from astropy.io import fits
from astropy.time import Time
from scipy.interpolate import interp1d

def get_nirspec_good_fixedslit_names():
    return('S200A1', 'S200A2', 'S400A1', 'S1600A1')

def get_keywords_from_jwst_header2(thisfile):   # common keywords I need all the time
    keys1     = ['targ_ra', 'targ_dec', 'EFFEXPTM',  'expmid', 'PROGRAM', 'OBS_ID', 'FILTER', 'GRATING', 'READPATT', 'TARGTYPE', 'TARGNAME', 'READPATT']
    extens1   = [  0,            0,          0,           0,          0,      0,        0,         0,         0,           0,          0,       0]
    keys2     = ['DATE-OBS', 'velosys', 'SLTNAME']
    extens2   = [  0,           1,           1   ] 
    keys     = keys1 + keys2
    extens  = extens1 + extens2
    result  = {}
    for ii, thiskey in enumerate(keys):
        try:
            result[thiskey] =  gethead(thisfile, thiskey, extension=extens[ii])
        except:  result[thiskey] = '-99'
    result['mjd'] = result['expmid']  # synonym
    # convert MJD of midpoint of observation to more useful formats
    result['datestamp']    = Time(result['mjd'], format='mjd').to_value('isot')
    result['decimalyear']  = Time(result['mjd'], format='mjd').to_value('decimalyear')
    from pathlib import Path
    import json
    from nirspec_ninja.artifacts.pipeline import read_score_metadata
    sidecar = Path(thisfile).with_suffix('.artifacts.json')
    combined = sidecar.exists() and 'input_count' in json.loads(sidecar.read_text())
    with fits.open(thisfile) as hdus:
        result.update(read_score_metadata(hdus[0].header, combined=combined))
    return(result)
def mask_nans(data) :  # return a numpy masked array, with nans masked
    return(np.ma.masked_array(data,np.isnan(data)))

def gethead(imagefile, header_keyword, extension=1) :
    # return a header keyword from a file, similar to wcstools gethead
    with fits.open(imagefile) as hdu:
        return hdu[extension].header[header_keyword]

def find_science_extensions(infile, kind='SCI'):  # much faster
    with fits.open(infile) as hdulist:
        names = [hdu.name for hdu in hdulist]
        indexes = [i for i, x in enumerate(names) if x == kind]
    return(indexes)
def rebin_spec_new(wave, specin, new_wave, fill=np.nan, return_masked=False):
    # Rebin spectra (wave, specin) to new wavelength aray new_wave.  Fill values are nan.  If return_masked, then mask the nans
    f = interp1d(wave, specin, bounds_error=False, fill_value=fill)  # With these settings, writes NaN to extrapolated regions
    new_spec = f(new_wave)
    if return_masked :  return(mask_nans(new_spec))
    else :              return(new_spec)

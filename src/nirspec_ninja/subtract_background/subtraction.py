"""Subtract a supplied expected-background spectrum from an IFU cube.

Moved from Brian's LEGGOS reduction workflow; original function name and
calculation preserved. The calling workflow selects the background file.
"""
import numpy as np
from astropy.io import fits
from scipy.interpolate import interp1d

def exp_bg_sub(cubefile, bgfile, outfile=None):
    # first, load the expected background data:
    exp = np.loadtxt(bgfile)
    expwl = np.array([line[0] for line in exp])
    exptot = np.array([line[1] for line in exp])
    # and make an interp1d object sdo we can get the background at any given point:
    bginterp = interp1d(expwl, exptot)
    # now load the data cube:
    hdu = fits.open(cubefile)
    data = hdu[1].data
    head1 = hdu[1].header
    cubewl = np.arange(head1['CRVAL3'],
                 head1['CRVAL3']+(head1['CDELT3']*len(data)),
                 head1['CDELT3'])
    if len(cubewl) > len(data):
        cubewl = cubewl[:-1]
    # next, evaluate background on data wavelength points:
    bgcalculated = bginterp(cubewl)
    bgcube = np.tile(bgcalculated, (data.shape[2],data.shape[1],1)).T # cubeify!
    # and subtract!
    #print(len(data), len(cubewl))
    bgsub_cube = data - bgcube
    # save or return result:
    if outfile:
        hdu[1].data = bgsub_cube
        hdu[1].header['comment'] = 'Expected Background Subtracted'
        hdu.writeto(outfile, overwrite=True)
        hdu.close()
    else:
        hdu.close()
        return bgsub_cube


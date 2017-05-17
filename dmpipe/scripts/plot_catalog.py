#!/usr/bin/env python
#

"""

"""

__facility__ = "plot_catalog"
__abstract__ = __doc__
__author__   = "E. Charles"
__date__     = "$Date: 2016/07/20 23:11:17 $"
__version__  = "$Revision: 1.1 $"
__release__  = "$Name:  $"


import os
import numpy as np
import yaml

from dmpipe import dmp_roi
from astropy import coordinates
import astropy.wcs as pywcs
import astropy.io.fits as pf
from astropy import table
from fermipy import utils



def create_wcs(skydir, coordsys='CEL', projection='AIT',
               cdelt=1.0, crpix=(1.,1.), naxis=2, energies=None):
    
    w = pywcs.WCS(naxis=naxis)
    
    if coordsys == 'CEL':
        w.wcs.ctype[0] = 'RA---%s' % (projection)
        w.wcs.ctype[1] = 'DEC--%s' % (projection)
        w.wcs.crval[0] = skydir.icrs.ra.deg
        w.wcs.crval[1] = skydir.icrs.dec.deg
    elif coordsys == 'GAL':
        w.wcs.ctype[0] = 'GLON-%s' % (projection)
        w.wcs.ctype[1] = 'GLAT-%s' % (projection)
        w.wcs.crval[0] = skydir.galactic.l.deg
        w.wcs.crval[1] = skydir.galactic.b.deg
    else:
        raise Exception('Unrecognized coordinate system.')

    w.wcs.crpix[0] = crpix[0]
    w.wcs.crpix[1] = crpix[1]
    w.wcs.cdelt[0] = -cdelt
    w.wcs.cdelt[1] = cdelt

    w = pywcs.WCS(w.to_header())
    if naxis == 3 and energies is not None:
        w.wcs.crpix[2] = 1
        w.wcs.crval[2] = 10 ** energies[0]
        w.wcs.cdelt[2] = 10 ** energies[1] - 10 ** energies[0]
        w.wcs.ctype[2] = 'Energy'

    return w





if __name__ == '__main__':
    
    roi_size_deg = 8.0 # in deg                                              

    npix = (1800,900)
    out_array = np.zeros(npix)

    skydir_gc = coordinates.SkyCoord(0.,0.,frame=coordinates.Galactic,unit="deg")
    wcs_gc = create_wcs(skydir_gc,coordsys='GAL',projection="AIT",cdelt=0.2,crpix=(900.5,450.5))
    
    t_all = table.Table.read("new_srcs_filtered.fits")
    glat_all = t_all['GLAT']
    glon_all = t_all['GLON']
    pix_crds_all = wcs_gc.wcs_world2pix(glon_all,glat_all,1)

    clusters = yaml.load(open("clustered_idx_dict.yaml"))

    from matplotlib import pyplot as plt

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xlim(0,npix[0])
    ax.set_ylim(0,npix[1])
    
    for k,v in clusters.items():
        if len(v) < 10:
            continue
        clust = [k] + v
        print(clust)
        ax.plot(pix_crds_all[0][clust],pix_crds_all[1][clust])

    ax.plot(pix_crds_all[0],pix_crds_all[1],'r,')

    
    
        

    

    

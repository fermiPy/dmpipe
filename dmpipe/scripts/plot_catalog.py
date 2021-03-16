#!/usr/bin/env python
#


# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Top level script to make a castro plot in mass / sigmav space
"""


import yaml

from astropy import coordinates
import astropy.wcs as pywcs
from astropy import table

from matplotlib import pyplot as plt


def create_wcs(skydir, coordsys='CEL', projection='AIT',
               cdelt=1.0, crpix=(1., 1.), naxis=2, energies=None):
    """ Build a WCS object """
    wout = pywcs.WCS(naxis=naxis)

    if coordsys == 'CEL':
        wout.wcs.ctype[0] = 'RA---%s' % (projection)
        wout.wcs.ctype[1] = 'DEC--%s' % (projection)
        wout.wcs.crval[0] = skydir.icrs.ra.deg
        wout.wcs.crval[1] = skydir.icrs.dec.deg
    elif coordsys == 'GAL':
        wout.wcs.ctype[0] = 'GLON-%s' % (projection)
        wout.wcs.ctype[1] = 'GLAT-%s' % (projection)
        wout.wcs.crval[0] = skydir.galactic.l.deg
        wout.wcs.crval[1] = skydir.galactic.b.deg
    else:
        raise Exception('Unrecognized coordinate system.')

    wout.wcs.crpix[0] = crpix[0]
    wout.wcs.crpix[1] = crpix[1]
    wout.wcs.cdelt[0] = -cdelt
    wout.wcs.cdelt[1] = cdelt

    wout = pywcs.WCS(wout.to_header())
    if naxis == 3 and energies is not None:
        wout.wcs.crpix[2] = 1
        wout.wcs.crval[2] = 10 ** energies[0]
        wout.wcs.cdelt[2] = 10 ** energies[1] - 10 ** energies[0]
        wout.wcs.ctype[2] = 'Energy'

    return wout


def main():
    """ Hook for command line interface """
    npix = (1800, 900)

    skydir_gc = coordinates.SkyCoord(0., 0., frame=coordinates.Galactic, unit="deg")
    wcs_gc = create_wcs(skydir_gc, coordsys='GAL', projection="AIT",
                        cdelt=0.2, crpix=(900.5, 450.5))

    t_all = table.Table.read("new_srcs_filtered.fits")
    glat_all = t_all['GLAT']
    glon_all = t_all['GLON']
    pix_crds_all = wcs_gc.wcs_world2pix(glon_all, glat_all, 1)

    clusters = yaml.load(open("clustered_idx_dict.yaml"))

    fig = plt.figure()
    axout = fig.add_subplot(111)
    axout.set_xlim(0, npix[0])
    axout.set_ylim(0, npix[1])

    for key, val in list(clusters.items()):
        if len(val) < 10:
            continue
        clust = [key] + val
        print(clust)
        axout.plot(pix_crds_all[0][clust], pix_crds_all[1][clust])

    axout.plot(pix_crds_all[0], pix_crds_all[1], 'r,')


if __name__ == '__main__':
    main()

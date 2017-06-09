#!/usr/bin/env python
#

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Top levels scripts to extract data from an all-sky analysis and produce a map
"""

import argparse

import numpy as np

from astropy import coordinates
import astropy.io.fits as pf

from fermipy import skymap
from fermipy import wcs_utils

from dmpipe import dmp_roi

def main():
    """ Hook for command line interface """
    npix = (1800, 900)

    # Argument defintion
    usage = "usage: %(prog)s [input]"
    description = "Collect all the new source"

    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument('--input', '-i', default='roi_set.yaml', help='ROI set definition file.')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'), help='Output file.')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file.')
    parser.add_argument('--filestr', '-f', default="tscube.fits",
                        help="Name of file within each ROI sub-directory")
    parser.add_argument('--extname', '-e', default=None,
                        help="Name of the HDU with the column to extract")
    parser.add_argument('--colname', '-c', default=None, help="Name of the column to extract")
    parser.add_argument('--nbin', type=int, default=None, help="Number of bins per pixel")
    parser.add_argument('--size', '-s', type=float, default=0.2, help='Pixel size (degrees)')
    parser.add_argument('--proj', '-p', default="AIT", help='Projection')

    # Argument parsing
    args = parser.parse_args()

    # Make a WCS, and get the corresponding mesh of coordinates
    # Make a WCS, and get the corresponding mesh of coordinates
    if args.nbin is None:
        npix = (int(np.ceil(360. / args.size)), int(np.ceil(180. / args.size)))
    else:
        npix = (int(np.ceil(360. / args.size)), int(np.ceil(180. / args.size)), args.nbin)
    crpix = ((npix[0] + 1.) / 2., (npix[1] + 1.) / 2.)

    crpix = ((npix[0] + 1.) / 2., (npix[1] + 1.) / 2.)

    skydir_gc = coordinates.SkyCoord(0., 0., frame=coordinates.Galactic, unit="deg")
    wcs_gc = wcs_utils.create_wcs(skydir_gc, coordsys='GAL',
                                  projection=args.proj, cdelt=args.size, crpix=crpix)

    pix_dirs = np.dstack(np.meshgrid(np.arange(1, npix[0] + 1),
                                     np.arange(1, npix[1] + 1))).swapaxes(0, 1).reshape((npix[0] * npix[1], 2))
    sky_crds = wcs_gc.wcs_pix2world(pix_dirs, 1)

    # Make the roi_set object
    roi_set, basedir = dmp_roi.DMROISet.create_from_yaml(args.input)

    kwargs = dict(tabname=args.extname,
                  colname=args.colname,
                  nbin=args.nbin)

    # extract the data
    out_array = roi_set.extract_column_from_tables(npix, sky_crds, args.filestr,
                                                   basedir=basedir, **kwargs)

    # Write the output
    outmap = skymap.Map(out_array.T, wcs_gc)
    hdulist = pf.HDUList([outmap.create_primary_hdu()])
    hdulist.writeto(args.output)


if __name__ == '__main__':
    main()

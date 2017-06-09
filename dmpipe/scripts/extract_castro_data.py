#!/usr/bin/env python
#

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Top levels scripts to extract castro data from an all-sky analysis
"""

import os
import argparse

import numpy as np
import yaml

from astropy import table

from fermipy import fits_utils

from dmpipe import dmp_roi
from dmpipe.dm_target import DMTargetFactory


def read_targets(filepath):
    """ Read a set of targets from a fits file """
    return read_targets_from_fits(filepath)


def read_targets_from_fits(fitsfile):
    """ Read a set of targets from a fits file """
    tab = table.Table.read(fitsfile)
    mask = np.zeros(len(tab), bool)
    key_col = tab['key']
    for i in range(len(tab)):
        mask[i] = key_col[i].find('point') != -1
    tab_mask = tab[mask]
    coords = np.ndarray((len(tab_mask), 2))
    coords[:, 0] = tab_mask['glon']
    coords[:, 1] = tab_mask['glat']

    out_dict = {'targets': tab_mask['target'],
                'coordinates': coords}
    return out_dict


def read_targets_from_yaml(yamlfile):
    """ Read a set of targets from a yaml file """
    din = yaml.load(yamlfile)
    coords = np.ndarray((len(din), 2))
    for i, (key, val) in enumerate(din.items()):
        coords[i, 0] = val['l']
        coords[i, 1] = val['b']

    out_dict = {'targets': din.keys(),
                'coordinates': coords}
    return out_dict


def add_columns(out_table, in_table, col_names):
    """ Add columnes to a table """
    for col in col_names:
        out_table.add_column(in_table[col])



def main():
    """ Hook for command line access """
    # Argument defintion
    usage = "usage: %(prog)s [input]"
    description = "Collect all the new source"

    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument('--input', '-i', default='roi_set.yaml', help='ROI set definition file.')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'), help='Output file.')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file.')
    parser.add_argument('--targets', '-t', type=str, help='Target file.')
    parser.add_argument('--filestr', '-f', default="tscube.fits",
                        help="Name of file within each ROI sub-directory")

    # Argument parsing
    args = parser.parse_args()

    # Read the target file
    targ_type = os.path.splitext(args.targets)[1]
    print targ_type
    if targ_type in ['.fit', '.fits']:
        targets = DMTargetFactory.read_table(args.targets)
        roster = None
    else:
        targets, roster = DMTargetFactory.make_table([args.targets])

    # Get the sky_crds
    sky_crds = DMTargetFactory.get_sky_crds(targets)

    # Make the roi_set object
    roi_set, basedir = dmp_roi.DMROISet.create_from_yaml(args.input)

    # extract the data
    out_table = roi_set.extract_table_data(sky_crds, args.filestr,
                                           basedir=basedir, tables=["SCANDATA", "FITDATA"])

    # add_names_column(out_table,targets['name'])
    col_names = ['name', 'ra', 'dec', 'distance', 'proftype', 'glat', 'glon', 'j_integ', 'd_integ']
    add_columns(out_table, targets, col_names)

    ebounds_table = roi_set.extract_single_table(args.filestr, basedir=basedir, table="EBOUNDS")

    # Write the output
    fits_utils.write_tables_to_fits(args.output, [out_table, ebounds_table],
                                    clobber=args.clobber, namelist=["SCANDATA", "EBOUNDS"])


if __name__ == '__main__':
    main()

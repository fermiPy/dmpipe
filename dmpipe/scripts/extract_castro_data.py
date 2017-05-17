#!/usr/bin/env python
#

"""

"""

__facility__ = "extract_castro_data"
__abstract__ = __doc__
__author__   = "E. Charles"
__date__     = "$Date: 2016/07/20 23:11:17 $"
__version__  = "$Revision: 1.1 $"
__release__  = "$Name:  $"


import sys
import os
import numpy as np
import argparse
import yaml

from astropy import coordinates
from astropy import table
import astropy.wcs as pywcs
import astropy.io.fits as pf

from fermipy import utils
from fermipy import fits_utils
from fermipy import skymap
from fermipy import wcs_utils

from dmpipe import dmp_roi
from dmfit.dm_target import dm_target_factory


def read_targets(filepath):
    return read_targets_from_fits(filepath)

def read_targets_from_fits(fitsfile):
    """
    """
    tab = table.Table.read(fitsfile)
    mask = np.zeros(len(tab),bool)
    key_col = tab['key']
    for i in range(len(tab)):
        mask[i] = key_col[i].find('point') != -1
    tab_mask = tab[mask]
    coords = np.ndarray((len(tab_mask),2))
    coords[:,0] = tab_mask['glon']
    coords[:,1] = tab_mask['glat']
    
    out_dict = {'targets':tab_mask['target'],
                'coordinates':coords}
    return out_dict


def read_targets_from_yaml(yamlfile):
    """
    """
    d = yaml.load(yamlfile)
    coords = np.ndarray((len(d),2))
    for i,(k,v) in enumerate(d.items()):
        coords[i,0] = v['l']
        coords[i,1] = v['b']
        pass

    out_dict = {'targets':d.keys(),
                'coordinates':coords}
    return out_dict


def add_columns(out_table,in_table,col_names):
    """
    """
    for c in col_names:
        out_table.add_column(in_table[c])


if __name__ == '__main__':
    
    # Argument defintion
    usage = "usage: %(prog)s [input]"    
    description = "Collect all the new source"
    
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('--input', '-i', default='roi_set.yaml',help='ROI set definition file.')
    parser.add_argument('--output', '-o', type=argparse.FileType('w'), help='Output file.')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file.')
    parser.add_argument('--targets','-t', type=str, help='Target file.')
    parser.add_argument('--filestr','-f', default="tscube.fits", help="Name of file within each ROI sub-directory")
    
    # Argument parsing
    args = parser.parse_args()

    # Read the target file
    targ_type = os.path.splitext(args.targets)[1]
    print targ_type
    if targ_type in ['.fit','.fits']:
        targets = dm_target_factory.read_table(args.targets)
        roster = None
    else:
        targets,roster = dm_target_factory.make_table([args.targets])

    # Get the sky_crds
    sky_crds =  dm_target_factory.get_sky_crds(targets)
    

    # Make the roi_set object
    roi_set,basedir = dmp_roi.dmp_roi_set.create_from_yaml(args.input)

    # extract the data
    out_table = roi_set.extract_table_data(sky_crds,args.filestr,   
                                           basedir=basedir,tables=["SCANDATA","FITDATA"])

    #add_names_column(out_table,targets['name'])
    col_names = ['name','ra','dec','distance','proftype','glat','glon','j_integ','d_integ']
    add_columns(out_table,targets,col_names)

    ebounds_table = roi_set.extract_single_table(args.filestr,basedir=basedir,table="EBOUNDS")

    # Write the output
    fits_utils.write_tables_to_fits(args.output,[out_table,ebounds_table],
                                    clobber=args.clobber,namelist=["SCANDATA","EBOUNDS"])



  
    
        

    

    

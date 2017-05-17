#!/usr/bin/env python
#

"""

"""

__facility__ = "extract_castro_dsphs"
__abstract__ = __doc__
__author__   = "E. Charles"
__date__     = "$Date: 2016/07/27 18:57:04 $"
__version__  = "$Revision: 1.2 $"
__release__  = "$Name:  $"


import sys
import os
import numpy as np
import argparse
import yaml

from astropy import table
from fermipy import castro
from dmfit.dm_target import dm_target_factory
from fermipy import fits_utils


def make_castro_from_name(name):
    """
    """
    path = os.path.join('matt_results',name,'pwl','2.0',"%s_sed.yaml"%name)
    return make_castro_from_dsph(path)


def make_castro_from_dsph(yamlfile):
    """
    """
    d = yaml.load(open(yamlfile))
    emins = np.array([b['emin'] for b in d])
    emaxs = np.array([b['emax'] for b in d])
    #convs = np.array([b['eflux2npred'] for b in d])
    normlist = np.array([b['eflux'] for b in d])    
    loglist = -1.*np.array([b['logLike'] for b in d])
    ref = np.ones(emins.shape)

    sd = castro.SpecData(emins,emaxs,ref,ref,ref,ref)
    norm_type = 'EFLUX'
    
    return castro.CastroData(normlist,loglist,sd,norm_type)


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
    parser.add_argument('--output', '-o', type=argparse.FileType('w'), help='Output file.')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file.')
    parser.add_argument('--targets','-t', type=str, help='Target file.')
    
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

    tab_e = None
    tab_s = None
    for t in targets:
        c = make_castro_from_name(t['name'])
        tab_s_targ = c.build_scandata_table()
        print (t['name'],len(tab_s_targ))
        if tab_s is None:
            tab_s = tab_s_targ
        else:
            tab_s.add_row(tab_s_targ[0])
        if tab_e is None:
            tab_e = c.refSpec.build_ebound_table()

    col_names = ['name','ra','dec','distance','proftype','glat','glon','j_integ','d_integ']
    add_columns(tab_s,targets,col_names)

    fits_utils.write_tables_to_fits(args.output,[tab_s,tab_e],namelist=["SCANDATA","EBOUNDS"],clobber=args.clobber)


  
    
        

    

    

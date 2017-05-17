#!/usr/bin/env python
#

"""
Interface to Dark Matter spectra
"""

import os
import sys
import numpy as np

from dmfit.dm_spectral import dm_spec_table, dm_castro_data, make_castro_dm
from dmfit.dm_fit_spec import DMFitFunction


import argparse
import yaml
from astropy import table

# Argument defintion
usage = "usage: %(prog)s [input]"    
description = "Collect all the new source"
    
parser = argparse.ArgumentParser(usage=usage,description=description)
parser.add_argument('--spec', '-s', default="dm_spec_2.fits",help='Spectra table')
parser.add_argument('--castro', '-c', default="dsph_castro.fits",help='Target castro data')
parser.add_argument('--output', '-o', default=None, help='Output file.')
parser.add_argument('--clobber', action='store_true', help='Overwrite output file.')

# Argument parsing
args = parser.parse_args()

#CHANNELS = ['ee','mumu','tautau','bb','tt','gg','ww','zz','cc','uu','dd','ss']
CHANNELS = ['bb']
norm_type = 'EFLUX'

spec_table = dm_spec_table.create_from_fits(args.spec)
castro_shape = (13,100)

tab_s = table.Table.read(args.castro,"SCANDATA")
tab_e = table.Table.read(args.castro,"EBOUNDS")

targ_names = tab_s['name']
nlist = [t for t in targ_names]
nlist.append('stacked')

for chan in CHANNELS:
    print "Channel %s: "%chan
    chan_idx = DMFitFunction.channel2int(chan)

    clist_flux = {}
    clist_no_prior = {}
    clist_prior = {}

    tlist_no_prior = []
    tlist_prior = []

    for tname in targ_names:
        sys.stdout.write('.')
        sys.stdout.flush()
        r = make_castro_dm(tname,spec_table,tab_s,tab_e,channel=4)
    
        clist_flux[tname] = r[1]
        clist_no_prior[tname] = r[2]
        clist_prior[tname] = r[3]
        
        tlist_no_prior.append(r[2].build_scandata_table())
        tlist_prior.append(r[3].build_scandata_table())
        
    sys.stdout.write("!\n")

    stacked_no_prior = dm_castro_data.create_from_stack(castro_shape,clist_no_prior.values())
    stacked_prior = dm_castro_data.create_from_stack(castro_shape,clist_prior.values())

    tlist_no_prior.append(stacked_no_prior)
    tlist_prior.append(stacked_prior)
    
    fits_utils.write_tables_to_fits("%s_%s_no_prior.fits"%(args.output,chan),tlist_no_prior,namelist=nlist,clobber=args.clobber)
    fits_utils.write_tables_to_fits("%s_%s_prior.fits"%(args.output,chan),tlist_prior,namelist=nlist,clobber=args.clobber)

    

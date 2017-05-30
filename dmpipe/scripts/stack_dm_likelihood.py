#!/usr/bin/env python
#

"""
Interface to Dark Matter spectra
"""

import os
import sys
import argparse

import yaml

from astropy.table import Table

from fermipy import fits_utils
from fermipy.castro import CastroData

from fermipy.utils import load_yaml, write_yaml

from dmpipe.dm_spectral import DMCastroData


def stack_roster(roster_name, rost, basedir, channels, jprior_key):
    """
    """
    component_dict = {}
    out_dict = {}
    for chan in channels:
        component_dict[chan] = []    

    for target_key in rost:
        tokens = target_key.split(':')
        target_name = tokens[0]
        target_version = tokens[1]
        target_dir = os.path.join(basedir, target_name)
        dmlike_path = os.path.join(target_dir, "dmlike_%s_%s.fits"%(target_version, jprior_key))
        tab_m = Table.read(dmlike_path, hdu="MASSES")
        for chan in channels:
            try:
                tab_s = Table.read(dmlike_path, hdu=chan)
            except KeyError:
                continue
            dm_castro = DMCastroData.create_from_tables(tab_s, tab_m)
            component_dict[chan].append(dm_castro)
            
    for chan, comps in component_dict.items():
        if len(comps) == 0:
            continue
        stacked = DMCastroData.create_from_stack(comps)
        out_dict[chan] = stacked

    return out_dict


def write_stacked(basedir, roster_name, stacked_dict, jprior_key, clobber):
    """
    """
    outdir = os.path.join(basedir, "stacked")
    try:
        os.makedirs(outdir)
    except OSError:
        pass        
    outpath = os.path.join(outdir, "results_%s_%s.fits"%(roster_name, jprior_key))
    print("Writing stacked results %s"%outpath)
    channels = stacked_dict.keys()
    t_list = []
    n_list = []
    mass_table = None
    for chan in channels:
        stacked = stacked_dict[chan]
        if mass_table is None:
            mass_table = stacked.build_mass_table()
        t_list.append(stacked.build_scandata_table())
        n_list.append(chan)
    t_list.append(mass_table)
    n_list.append("MASSES")
    fits_utils.write_tables_to_fits(outpath, t_list,
                                    clobber=clobber, namelist=n_list)


def stack_rosters(roster_dict, basedir, channels, jprior_key, clobber):
    """
    """
    for roster_name, rost in roster_dict.items():
        stacked_dict = stack_roster(roster_name, rost, basedir, channels, jprior_key)
        write_stacked(basedir, roster_name, stacked_dict, jprior_key, clobber)


        
def main():
    """
    """
    # Argument defintion
    usage = "usage: %(prog)s [input]"    
    description = "Collect all the new source"
    
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('--jprior', '-j', default=None, type=str, help='Type of Prior on J-factor')
    parser.add_argument('--dir', '-d', required=True, help='Name of top-level directory')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file.')
    
    # Argument parsing
    args = parser.parse_args()
    channels = ['ee','mumu','tautau','bb','tt','gg','ww','zz','cc','uu','dd','ss']

    roster_dict = load_yaml(os.path.join(args.dir, 'roster_list.yaml'))
    stack_rosters(roster_dict, args.dir, channels, args.jprior, args.clobber)   
  

if __name__=='__main__':
    main()


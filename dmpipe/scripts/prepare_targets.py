#!/usr/bin/env python
#

"""
Interface to Dark Matter spectra
"""

import os
import sys
import numpy as np


import argparse
import yaml

from dmsky.roster import RosterLibrary
from fermipy.utils import load_yaml, write_yaml


def write_target_dirs(basedir, roster_dict, base_config):
    """
    """
    target_dict = {}
    
    target_info_dict = {}
    roster_info_dict = {}

    try:
        os.makedirs(basedir)
    except OSError:
        pass

    for roster_name, rost in roster_dict.items():
        tlist = []
        for target_name, target in rost.items():
            target_key = "%s:%s"%(target_name, target.version)
            print("Writing %s"%(target_key))
            tlist.append(target_key)
            if target_info_dict.has_key(target_name):
                target_info_dict[target_name].append(target.version)
            else:
                target_info_dict[target_name] = [target.version]
            target_dir = os.path.join(basedir, target_name)
            target_config_path = os.path.join(target_dir, 'config_baseline.yaml')    
            jmap_path = os.path.join(target_dir, 'profile_%s.fits'%target.version)
            profile_path = os.path.join(target_dir, 'profile_%s.yaml'%target.version)

            if target_dict.has_key(target_name):
                # Already made the config for this target
                target_config = target_dict[target_name]
            else:
                # Make the config for this target
                try:
                    os.makedirs(target_dir)
                except OSError:
                    pass
                target_config = base_config.copy()
                target_config['selection']['ra'] = target.ra
                target_config['selection']['dec'] = target.dec
                target_dict[target_name] = target_config
                write_yaml(target_config, target_config_path)
               

            profile_data = target.profile.copy()
            target.write_jmap_wcs(jmap_path, clobber=True)
            profile_data['j_integ'] = target.j_integ
            profile_data['j_sigma'] = target.j_sigma
            profile_data['j_map_file'] = jmap_path
            write_yaml(profile_data, profile_path)

        roster_info_dict[roster_name] = tlist

    write_yaml(roster_info_dict, os.path.join(basedir, 'roster_list.yaml'))
    write_yaml(target_info_dict, os.path.join(basedir, 'target_list.yaml'))


def main():
    
    # Argument defintion
    usage = "usage: %(prog)s [input]"    
    description = ""
    
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('--roster', '-r', required=True, help='Roster to build targets for')
    parser.add_argument('--baseconfig', '-b', required=True, help='Baseline configuration yaml file')
    parser.add_argument('--topdir', '-t', required=True, help='Top level output directory.')

    # Argument parsing
    args = parser.parse_args()
    
    roster_lib = RosterLibrary()
    roster_dict  = {}
    rost = roster_lib.create_roster(args.roster)
    roster_dict[args.roster] = rost

    base_config = load_yaml(args.baseconfig)

    write_target_dirs(args.topdir, roster_dict, base_config)


if __name__ == "__main__":
    main()
    

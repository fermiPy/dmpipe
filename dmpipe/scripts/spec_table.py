#!/usr/bin/env python
#

"""
Build a table with DM spectra
"""

import os
import sys
import numpy as np


import argparse
import yaml

from dmpipe.dm_spectral import DMSpecTable


def main():
    
    # Argument defintion
    usage = "usage: %(prog)s [input]"    
    description = "Build a table with DM spectra"
    
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('--config', '-c', required=True, help='Baseline configuration yaml file')
    parser.add_argument('--outfile', '-o', required=True, help='Output file.')
    parser.add_argument('--clobber', action='store_true', help='Overwrite output file.')

    args = parser.parse_args()

    channels = ['ee','mumu','tautau','bb','tt','gg','ww','zz','cc','uu','dd','ss']
    masses = np.logspace(1,6,21)

    dm_spec_table = DMSpecTable.create_from_config(args.config, channels, masses)
    dm_spec_table.write_fits(args.outfile, args.clobber)


if __name__ == "__main__":
    main()
    

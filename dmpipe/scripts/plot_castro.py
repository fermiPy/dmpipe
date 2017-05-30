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

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro


        
def main():
    """
    """
    # Argument defintion
    usage = "usage: %(prog)s [input]"    
    description = "Collect all the new source"
    
    parser = argparse.ArgumentParser(usage=usage,description=description)
    parser.add_argument('--input', '-i', required=True, help='Input FITS file')
    
    # Argument parsing
    args = parser.parse_args()

    castro_data = CastroData.create_from_sedfile(args.input)
    ylims = [1e-8, 1e-5]

    plot = plotCastro(castro_data, ylims)
    return plot


if __name__=='__main__':
    plot = main()


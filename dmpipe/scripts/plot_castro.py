#!/usr/bin/env python
#

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Top level script to make a castro plot
"""
from __future__ import absolute_import, division, print_function


import argparse

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro


def main():
    """ Hook for command line interface
    """
    # Argument defintion
    usage = "usage: %(prog)s [input]"
    description = "Collect all the new source"

    parser = argparse.ArgumentParser(usage=usage, description=description)
    parser.add_argument('--input', '-i', required=True, help='Input FITS file')
    parser.add_argument('--output', '-o', default=None, type=str, help='Output file')
 
    # Argument parsing
    args = parser.parse_args()

    castro_data = CastroData.create_from_sedfile(args.input)
    ylims = [1e-8, 1e-5]

    plot = plotCastro(castro_data, ylims)
    if args.output:
        plot[0].savefig(args.output)    
        return None
    return plot


if __name__ == '__main__':
    PLOT = main()

#!/usr/bin/env python
#

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Top level script to make a plot of the various DM spectra
"""
from __future__ import absolute_import, division, print_function

import sys
import os
import argparse
import numpy as np

from astropy.table import Table

from fermipy.utils import init_matplotlib_backend, load_yaml
init_matplotlib_backend()

from fermipy.jobs.chain import Link

from dmpipe.dm_spectral import DMSpecTable
from dmpipe.dm_plotting import plot_limits_from_arrays, plot_mc_truth, plot_dm_spectra_by_mass, plot_dm_spectra_by_channel

from dmpipe import defaults


class DMSpectraPlotter(Link):
    """Small class wrap an analysis script.
    
    This is useful for parallelizing analysis using the fermipy.jobs module.
    """

    default_options = dict(infile=defaults.generic['infile'],
                           outfile=defaults.generic['outfile'],
                           chan=defaults.common['chan'],
                           mass=defaults.common['mass'],
                           spec_type=defaults.common['spec_type'])

    def __init__(self, **kwargs):
        """C'tor
        """
        parser = argparse.ArgumentParser(usage="dmpipe-plot-dm-spectra [options]",
                                         description="Make a plot of the DM spectra")
        Link.__init__(self, kwargs.pop('linkname', 'plot-dm-spectra'),
                      parser=parser,
                      appname='dmpipe-plot-dm-spectra',
                      options=DMSpectraPlotter.default_options.copy(),
                      **kwargs)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)
        
        dm_spec_table = DMSpecTable.create_from_fits(args.infile)
        dm_plot_by_mass = plot_dm_spectra_by_mass(dm_spec_table, chan=args.chan, spec_type=args.spec_type)
        dm_plot_by_chan = plot_dm_spectra_by_channel(dm_spec_table, mass=args.mass, spec_type=args.spec_type)
        
        if args.outfile:
            dm_plot_by_mass[0].savefig(args.outfile.replace('.png', '_%s.png'%args.chan))
            dm_plot_by_chan[0].savefig(args.outfile.replace('.png', '_%1.FGeV.png'%args.mass))


def create_link_plot_dm_spectra(**kwargs):
    """Build and return a `Link` object that can invoke CastroPlotter"""
    spectra_plotter = DMSpectraPlotter(**kwargs)
    return spectra_plotter


def main_single():
    """ Hook for command line interface
    """
    spectra_plotter = DMSpectraPlotter()
    spectra_plotter.run_analysis(sys.argv[1:])

if __name__ == '__main__':
    main_single()

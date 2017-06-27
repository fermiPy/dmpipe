#!/usr/bin/env python
#

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Top level script to make a castro plot
"""
from __future__ import absolute_import, division, print_function

import sys
import os
import argparse

from fermipy.utils import init_matplotlib_backend, load_yaml
init_matplotlib_backend('Agg')

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro

from fermipy.jobs.chain import Link
from fermipy.jobs.scatter_gather import ConfigMaker
from fermipy.jobs.lsf_impl import build_sg_from_link


class CastroPlotter(Link):
    """Small class wrap an analysis script.
    
    This is useful for parallelizing analysis using the fermipy.jobs module.
    """

    default_options = dict(input=(None, 'Input file path', str),
                           output=(None, 'Output file path', str),)

    def __init__(self, **kwargs):
        """C'tor
        """
        parser = argparse.ArgumentParser(usage="dmpipe-plot-castro [options]",
                                         description="Make castro plots")
        Link.__init__(self, kwargs.pop('linkname', 'plot-castro'),
                      parser=parser,
                      appname='dmpipe-plot-castro',
                      options=CastroPlotter.default_options.copy(),
                      **kwargs)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)
        castro_data = CastroData.create_from_sedfile(args.input)
        ylims = [1e-8, 1e-5]
        
        plot = plotCastro(castro_data, ylims)
        if args.output:
            plot[0].savefig(args.output)    
            return None
        return plot
  


class ConfigMaker_PlotCastro(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(targetlist=('target_list.yaml', 'Yaml file with list of targets', str),
                           topdir=(None, 'Top level directory', str))

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_PlotCastro.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        input_config = {}
        job_configs = {}
        output_config = {}

        topdir = args['topdir']
        targets_yaml = os.path.join(topdir, args['targetlist'])

        try:
            targets = load_yaml(targets_yaml)
        except IOError:
            targets = {}

        for target_name, target_list in targets.items():
            for targ_prof in target_list:
                targ_key = "%s_%s"%(target_name, targ_prof)
                input_path = os.path.join(topdir, target_name, 'sed_%s.fits'%targ_prof)
                output_path = os.path.join(topdir, target_name, 'sed_%s.png'%targ_prof)
                logfile = os.path.join(topdir, target_name, 'plot_castro_%s.log'%targ_prof)
                job_config = dict(input=input_path,
                                  output=output_path)
                job_configs[targ_key] = job_config
                
        return input_config, job_configs, output_config



def create_link_plot_castro(**kwargs):
    """Build and return a `Link` object that can invoke CastroPlotter"""
    castro_plotter = CastroPlotter(**kwargs)
    return castro_plotter


def create_sg_plot_castro(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""

    link = CastroPlotter(**kwargs)
    link.linkname = kwargs.pop('linkname', link.linkname)
    appname = kwargs.pop('appname', 'dmpipe-plot-castro-sg')

    lsf_args = {'W': 50,
                'R': 'rhel60'}

    usage = "%s [options]" % (appname)
    description = "Make castro plots for set of targets"

    config_maker = ConfigMaker_PlotCastro(link)
    lsf_sg = build_sg_from_link(link, config_maker,
                                lsf_args=lsf_args,
                                usage=usage,
                                description=description,
                                appname=appname,
                                **kwargs)
    return lsf_sg


def main_single():
    """ Hook for command line interface
    """
    castro_plotter = CastroPlotter()
    castro_plotter.run_analysis(sys.argv[1:])


def main_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_plot_castro()
    lsf_sg(sys.argv)



if __name__ == '__main__':
    PLOT = main_single()

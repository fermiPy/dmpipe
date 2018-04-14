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
from fermipy.jobs.scatter_gather import ConfigMaker, build_sg_from_link
from fermipy.jobs.lsf_impl import make_nfs_path, get_lsf_default_args, LSF_Interface

from dmpipe.name_policy import NameFactory
from dmpipe import defaults

NAME_FACTORY = NameFactory(basedir='.')


class CastroPlotter(Link):
    """Small class wrap an analysis script.
    
    This is useful for parallelizing analysis using the fermipy.jobs module.
    """

    default_options = dict(infile=defaults.generic['infile'],
                           outfile=defaults.generic['outfile'])

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
        castro_data = CastroData.create_from_sedfile(args.infile)
        ylims = [1e-8, 1e-5]
        
        plot = plotCastro(castro_data, ylims)
        if args.outfile:
            plot[0].savefig(args.outfile)   
            return None
        return plot
  


class ConfigMaker_PlotCastro(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_PlotCastro.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args)
        if targets_yaml is None:
            return job_configs

        targets = load_yaml(targets_yaml)

        for target_name, target_list in targets.items():
            for targ_prof in target_list:
                name_keys = dict(target_type=ttype,
                                 target_name=target_name,
                                 profile=targ_prof,
                                 fullpath=True)
                targ_key = "%s_%s"%(target_name, targ_prof)
                input_path = NAME_FACTORY.sedfile(**name_keys)
                output_path = input_path.replace('.fits', '.png')
                logfile = make_nfs_path(input_path.replace('.fits', '.log'))
                job_config = dict(infile=input_path,
                                  outfile=output_path,
                                  logfile=logfile)
                job_configs[targ_key] = job_config
                
        return job_configs



def create_link_plot_castro(**kwargs):
    """Build and return a `Link` object that can invoke CastroPlotter"""
    castro_plotter = CastroPlotter(**kwargs)
    return castro_plotter


def create_sg_plot_castro(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    appname = kwargs.pop('appname', 'dmpipe-plot-castro-sg')
    link = create_link_plot_castro(**kwargs)
    link.linkname = kwargs.pop('linkname', link.linkname)

    batch_args = get_lsf_default_args()    
    batch_args['lsf_args']['W'] = 50
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Make castro plots for set of targets"

    config_maker = ConfigMaker_PlotCastro(link)
    lsf_sg = build_sg_from_link(link, config_maker,
                                interface=batch_interface,
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

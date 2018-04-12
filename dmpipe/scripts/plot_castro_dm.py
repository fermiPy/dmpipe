#!/usr/bin/env python
#

# Licensed under a 3-clause BSD style license - see LICENSE.rst

"""
Top level script to make a castro plot in mass / sigmav space
"""
from __future__ import absolute_import, division, print_function

import sys
import os
import argparse

from astropy.table import Table

from fermipy.utils import init_matplotlib_backend, load_yaml
init_matplotlib_backend()

from fermipy.jobs.chain import Link
from fermipy.jobs.scatter_gather import ConfigMaker, build_sg_from_link
from fermipy.jobs.lsf_impl import get_lsf_default_args, LSF_Interface

from dmpipe.dm_spectral import DMCastroData
from dmpipe.dm_plotting import plot_dm_castro

from dmpipe.name_policy import NameFactory

NAME_FACTORY = NameFactory(basedir='.')

class CastroPlotterDM(Link):
    """Small class wrap an analysis script.
    
    This is useful for parallelizing analysis using the fermipy.jobs module.
    """

    default_options = dict(input=(None, 'Input file path', str),
                           output=(None, 'Output file path', str),
                           chan=(None, 'Analysis channel', str),
                           jprior=(None, 'Type of Prior on J-factor', str),)

    def __init__(self, **kwargs):
        """C'tor
        """
        parser = argparse.ArgumentParser(usage="dmpipe-plot-dm [options]",
                                         description="Make castro plots")
        Link.__init__(self, kwargs.pop('linkname', 'plot-dm'),
                      parser=parser,
                      appname='dmpipe-plot-dm',
                      options=CastroPlotterDM.default_options.copy(),
                      **kwargs)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)
        tab_m = Table.read(args.input, hdu="MASSES")
        tab_s = Table.read(args.input, hdu=args.chan)
        dm_castro = DMCastroData.create_from_tables(tab_s, tab_m)
        
        dm_plot = plot_dm_castro(dm_castro)
        if args.output:
            dm_plot[0].savefig(args.output)
            return None
        return dm_plot


 
class ConfigMaker_PlotCastroDM(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(targetlist=('target_list.yaml', 'Yaml file with list of targets', str),
                           chan=(None, 'Analysis channel', str),
                           topdir=(None, 'Top level directory', str),
                           jprior=(None, 'Type of Prior on J-factor', str),)

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_PlotCastroDM.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        topdir = args['topdir']
        targets_yaml = os.path.join(topdir, args['targetlist'])

        try:
            targets = load_yaml(targets_yaml)
        except IOError:
            targets = {}

        if args['jprior'] is None or args['jprior'] == 'None' or args['jprior'] == 'none':
            j_prior_key = 'None'
        else:           
            j_prior_key = args['jprior']
            
        chan = args['chan']

        for target_name, target_list in targets.items():
            for targ_prof in target_list:
                name_keys = dict(target_type=topdir,
                                 target_name=target_name,
                                 profile=targ_prof,
                                 jprior=j_prior_key,
                                 fullpath=True)
                targ_key = "%s_%s"%(target_name, targ_prof)
                input_path = NAME_FACTORY.dmlikefile(**name_keys)
                output_path = input_path.replace('.fits', '.png')
                logfile = make_nfs_path(input_path.replace('.fits', 'log'))
                job_config = dict(input=input_path,
                                  output=output_path,
                                  j_prior=j_prior_key,
                                  logfile=logfile,
                                  chan=chan)
                job_configs[targ_key] = job_config
                
        return job_configs



def create_link_plot_dm(**kwargs):
    """Build and return a `Link` object that can invoke CastroPlotter"""
    castro_plotter = CastroPlotterDM(**kwargs)
    return castro_plotter


def create_sg_plot_dm(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    appname = kwargs.pop('appname', 'dmpipe-plot-dm-sg')
    link = create_link_plot_dm(**kwargs)
    link.linkname = kwargs.pop('linkname', link.linkname)

    batch_args = get_lsf_default_args()    
    batch_args['lsf_args']['W'] = 50
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Make castro plots for set of targets"

    config_maker = ConfigMaker_PlotCastroDM(link)
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
    castro_plotter = CastroPlotterDM()
    castro_plotter.run_analysis(sys.argv[1:])


def main_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_plot_dm()
    lsf_sg(sys.argv)


if __name__ == '__main__':
    DM_PLOT = main_single()

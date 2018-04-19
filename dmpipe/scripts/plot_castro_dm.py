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
from fermipy.jobs.lsf_impl import make_nfs_path, get_lsf_default_args, LSF_Interface

from dmpipe.dm_spectral import DMCastroData
from dmpipe.dm_plotting import plot_dm_castro

from dmpipe.name_policy import NameFactory
from dmpipe import defaults

NAME_FACTORY = NameFactory(basedir='.')

class CastroPlotterDM(Link):
    """Small class wrap an analysis script.
    
    This is useful for parallelizing analysis using the fermipy.jobs module.
    """

    default_options = dict(infile=defaults.generic['infile'],
                           outfile=defaults.generic['outfile'],
                           chan=defaults.common['chan'])

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
        tab_m = Table.read(args.infile, hdu="MASSES")
        tab_s = Table.read(args.infile, hdu=args.chan)
        dm_castro = DMCastroData.create_from_tables(tab_s, tab_m)
        
        dm_plot = plot_dm_castro(dm_castro)
        if args.outfile:
            dm_plot[0].savefig(args.outfile)
            return None
        return dm_plot


 
class ConfigMaker_PlotCastroDM(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           channels=defaults.common['channels'],
                           jpriors=defaults.common['jpriors'],
                           dry_run=defaults.common['dry_run'])

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

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args)
        if targets_yaml is None:
            return job_configs

        targets = load_yaml(targets_yaml)

        jpriors = args['jpriors']
        channels = args['channels']

        for target_name, target_list in targets.items():
            for targ_prof in target_list:
                for jprior in jpriors:
                    name_keys = dict(target_type=ttype,
                                     target_name=target_name,
                                     profile=targ_prof,
                                     jprior=jprior,
                                     fullpath=True)
                    input_path = NAME_FACTORY.dmlikefile(**name_keys)
                    for chan in channels:
                        targ_key = "%s:%s:%s:%s"%(target_name, targ_prof, jprior, chan)
                        output_path = input_path.replace('.fits', '_%s.png'%chan)
                        logfile = make_nfs_path(output_path.replace('.png', '.log'))
                        job_config = dict(infile=input_path,
                                          outfile=output_path,
                                          jprior=jprior,
                                          logfile=logfile,
                                          chan=chan)
                        job_configs[targ_key] = job_config
                
        return job_configs


class ConfigMaker_PlotStackedCastroDM(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           rosterlist=defaults.common['targetlist'],
                           channels=defaults.common['channels'],
                           jpriors=defaults.common['jpriors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],                           
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_PlotStackedCastroDM.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}
        
        ttype = args['ttype']
        (roster_yaml, sim) = NAME_FACTORY.resolve_rosterfile(args)
        if roster_yaml is None:
            return job_configs

        roster_dict = load_yaml(roster_yaml)
 
        jpriors = args['jpriors']
        channels = args['channels']

        for roster_name in roster_dict.keys():
            for jprior in jpriors:
                name_keys = dict(target_type=ttype,
                                 roster_name=roster_name,
                                 jprior=jprior,
                                 sim_name=sim,
                                 fullpath=True)
                
                for chan in channels:
                    targ_key = "%s:%s:%s"%(roster_name, jprior, chan)
                
                    if sim is not None:
                        seedlist = range(args['seed'], args['seed']+args['nsims'])
                    else:
                        seedlist = [None]
                        
                    for seed in seedlist:
                        if seed is not None:
                            name_keys['seed'] = "%06i"%seed                    
                            input_path = NAME_FACTORY.sim_resultsfile(**name_keys)
                            targ_key += "_%06i"%seed   
                        else:
                            input_path = NAME_FACTORY.resultsfile(**name_keys)

                        output_path = input_path.replace('.fits', '_%s.png'%chan)
                        logfile = make_nfs_path(output_path.replace('.png', '.log'))
                        job_config = dict(infile=input_path,
                                          outfile=output_path,
                                          jprior=jprior,
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

def create_sg_plot_stacked_dm(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    appname = kwargs.pop('appname', 'dmpipe-plot-stacked-dm-sg')
    link = create_link_plot_dm(**kwargs)
    link.linkname = kwargs.pop('linkname', link.linkname)

    batch_args = get_lsf_default_args()    
    batch_args['lsf_args']['W'] = 50
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Make castro plots for set of targets"

    config_maker = ConfigMaker_PlotStackedCastroDM(link)
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

def main_batch_stacked():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_plot_stacked_dm()
    lsf_sg(sys.argv)


if __name__ == '__main__':
    DM_PLOT = main_single()

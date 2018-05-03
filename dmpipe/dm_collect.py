#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Collect information for simulated realizations of an analysis
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import argparse
import numpy as np

from shutil import copyfile

#from dmsky.roster import RosterLibrary
from astropy.table import Table, Column, vstack

from fermipy.utils import load_yaml, write_yaml, init_matplotlib_backend

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro

from fermipy.jobs.utils import is_null, is_not_null
from fermipy.jobs.link import Link
from fermipy.jobs.scatter_gather import ConfigMaker, build_sg_from_link
from fermipy.jobs.slac_impl import make_nfs_path
from fermipy import fits_utils

from fermipy.jobs.target_collect import fill_output_table, vstack_tables, collect_summary_stats, add_summary_stats_to_table

from dmpipe.name_policy import NameFactory
from dmpipe import defaults

init_matplotlib_backend('Agg')

NAME_FACTORY = NameFactory(basedir=('.'))


def summarize_limits_results(limit_table):
    """Build a stats summary table for a table that has all the SED results """
    del_cols = ['UL_0.68', 'UL_0.95', 'MLES']
    stats_cols = ['UL_0.95', 'MLES']

    table_out = Table(limit_table[0])
    table_out.remove_columns(del_cols)
    add_summary_stats_to_table(limit_table, table_out, stats_cols)
    return table_out


class CollectLimits(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    appname = 'dmpipe-collect-limits'
    linkname_default = 'collect-limits'
    usage = '%s [options]'%(appname)
    description = "Collect Limits from simulations"

    default_options = dict(limitfile=defaults.generic['limitfile'],
                           specconfig=defaults.common['specconfig'],
                           outfile=defaults.generic['outfile'],
                           summaryfile=defaults.generic['summaryfile'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(CollectLimits, self).__init__(linkname, **init_dict)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        limitfile = args.limitfile
        first = args.seed
        last = first + args.nsims
        flist = [ limitfile.replace("_SEED.fits", "_%06i.fits"%seed) for seed in range(first, last) ]

        spec_config = load_yaml(args.specconfig)
        channels = spec_config['channels']

        outfile = args.outfile
        summaryfile = args.summaryfile

        hdus = channels + ['MASSES']

        out_tables, out_names = vstack_tables(flist, hdus) 

        if is_not_null(outfile):
            fits_utils.write_tables_to_fits(outfile, out_tables, namelist=out_names)

        if is_not_null(summaryfile):
            summary_tables = [summarize_limits_results(ot) for ot in out_tables[0:-1]]
            summary_tables.append( Table(out_tables[-1][0])  )
            fits_utils.write_tables_to_fits(summaryfile, summary_tables, namelist=out_names)



class CollectLimits_SG(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-collect-limits-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = CollectLimits

    job_time = 120

    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           specconifg=defaults.common['specconfig'],
                           jpriors=defaults.common['jpriors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           write_full=defaults.collect['write_full'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        super(CollectLimits_SG, self).__init__(link,
                                               options=kwargs.get('options',
                                                                  self.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args, require_sim_name=True)
        if targets_yaml is None:
            return job_configs

        specconfig = NAME_FACTORY.resolve_specconfig(args)

        jpriors = args['jpriors']

        write_full = args.get('write_full', False)

        targets = load_yaml(targets_yaml)
        for target_name, profile_list in targets.items():
            for profile in profile_list:
                for jprior in jpriors:
                    if is_null(jprior):
                        jprior = 'none'
                    full_key = "%s:%s:%s:%s" % (target_name, profile, sim, jprior)
                    name_keys = dict(target_type=ttype, 
                                     target_name=target_name,
                                     sim_name=sim,
                                     profile=profile,
                                     jprior=jprior, 
                                     fullpath=True)
                    limitfile = NAME_FACTORY.sim_dmlimitsfile(**name_keys)
                    first = args['seed']
                    last = first + args['nsims'] - 1
                    outfile = limitfile.replace('_SEED.fits','_collected_%06i_%06i.fits'%(first, last))
                    logfile = make_nfs_path(outfile.replace('.fits', '.log'))
                    if not write_full:
                        outfile = None
                    summaryfile = limitfile.replace('_SEED.fits','_summary_%06i_%06i.fits'%(first, last))
                    job_config = dict(limitfile=limitfile,
                                      specconfig=specconfig,
                                      jprior=jprior, 
                                      outfile=outfile,
                                      summaryfile=summaryfile,
                                      logfile=logfile,
                                      nsims=args['nsims'],
                                      seed=args['seed'],
                                      dry_run=args['dry_run'])
                    job_configs[full_key] = job_config

        return job_configs


class CollectStackedLimits_SG(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-collect-stacked-limits-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = CollectLimits

    job_time = 120

    default_options = dict(ttype=defaults.common['ttype'],
                           rosterlist=defaults.common['targetlist'],
                           jpriors=defaults.common['jpriors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           write_full=defaults.collect['write_full'],
                           write_summary=defaults.collect['write_summary'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        super(CollectStackedLimits_SG, self).__init__(link,
                                                      options=kwargs.get('options',
                                                                         self.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (roster_yaml, sim) = NAME_FACTORY.resolve_rosterfile(args, require_sim_name=True)
        if roster_yaml is None:
            return job_configs

        specconfig = NAME_FACTORY.resolve_specconfig(args)

        jpriors = args['jpriors']

        write_full = args['write_full']

        roster_dict = load_yaml(roster_yaml)
        for roster_name in roster_dict.keys():
            for jprior in jpriors:
                if is_null(jprior):
                    jprior = 'none'
                full_key = "%s:%s:%s" % (roster_name, sim, jprior)
                name_keys = dict(target_type=ttype, 
                                 roster_name=roster_name,
                                 sim_name=sim,
                                 jprior=jprior, 
                                 fullpath=True)

                limitfile = NAME_FACTORY.sim_stackedlimitsfile(**name_keys)
                first = args['seed']
                last = first + args['nsims'] - 1
                outfile = limitfile.replace('_SEED.fits','_collected_%06i_%06i.fits'%(first, last))
                logfile = make_nfs_path(outfile.replace('.fits', '.log'))
                if not write_full:
                    outfile = None
                summaryfile = limitfile.replace('_SEED.fits','_summary_%06i_%06i.fits'%(first, last))
                job_config = dict(limitfile=limitfile,
                                  specconfig=specconfig,
                                  jprior=jprior, 
                                  outfile=outfile,
                                  summaryfile=summaryfile,
                                  logfile=logfile,
                                  nsims=args['nsims'],
                                  seed=args['seed'],
                                  dry_run=args['dry_run'])
                job_configs[full_key] = job_config

        return job_configs

def register_classes():
    CollectLimits.register_class()
    CollectLimits_SG.register_class()
    CollectStackedLimits_SG.register_class()

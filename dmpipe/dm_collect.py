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

from fermipy.jobs.chain import Link
from fermipy.jobs.scatter_gather import ConfigMaker, build_sg_from_link
from fermipy.jobs.lsf_impl import make_nfs_path, get_lsf_default_args, LSF_Interface
from fermipy import fits_utils

from dmpipe.name_policy import NameFactory
from dmpipe import defaults

init_matplotlib_backend('Agg')

try:
    from fermipy.gtanalysis import GTAnalysis
    HAVE_ST = True
except ImportError:
    HAVE_ST = False

NAME_FACTORY = NameFactory(basedir=('.'))


def fill_output_table(filelist, hdu, collist, nbins):
    """Fill the arrays from the files in filelsit"""
    nfiles = len(filelist)
    shape = (nbins, nfiles)
    outdict = {}
    for c in collist:
        outdict[c['name']] = np.ndarray(shape)

    sys.stdout.write('Working on %i files: '%nfiles)
    sys.stdout.flush()
    for i, f in enumerate(filelist):
        sys.stdout.write('.')
        sys.stdout.flush()
        tab = Table.read(f, hdu)
        for c in collist:
            cname = c['name']
            outdict[cname][:,i] = tab[cname].data
    sys.stdout.write('!\n')
    outcols = []
    for c in collist:
        cname = c['name']
        if c.has_key('unit'):
            col = Column(data=outdict[cname], name=cname, dtype=np.float, shape=nfiles, unit=c['unit'])
        else:
            col = Column(data=outdict[cname], name=cname, dtype=np.float, shape=nfiles)  
        outcols.append(col)
    tab = Table(data=outcols)
    return tab


def vstack_tables(filelist, hdus):
    """Fill the arrays from the files in filelsit"""
    nfiles = len(filelist)
    out_tables = []
    out_names = []
    for hdu in hdus:
        sys.stdout.write('Working on %i files for %s: '%(nfiles, hdu))
        sys.stdout.flush()
        tlist = []
        for f in filelist:
            try:
                tab = Table.read(f, hdu)
                tlist.append(tab)
                sys.stdout.write('.')
            except KeyError:
                sys.stdout.write('x')
            sys.stdout.flush()            
        sys.stdout.write('!\n')
        if len(tlist) > 0:
            out_table = vstack(tlist)
            out_tables.append(out_table) 
            out_names.append(hdu)             
    return (out_tables, out_names)


def collect_summary_stats(data):
    """Collect summary statisitics from an array"""
    mean = np.mean(data, axis=0)
    std = np.std(data, axis=0)
    median = np.median(data, axis=0)
    q02,q16,q84,q97 = np.percentile(data, [2.5,16,84,97.5], axis=0)

    o = dict(mean=mean,
             std=std,
             median=median,
             q02=q02,
             q16=q16,
             q84=q84,
             q97=q97)
    
    return o


def add_summary_stats_to_table(table_in, table_out, colnames):
    """Collect summary statisitics from an input table and add them to an output table """
    for col in colnames:
        col_in = table_in[col]
        stats = collect_summary_stats(col_in.data)
        for k, v in stats.items():
            out_name = "%s_%s"%(col, k)
            col_out = Column(data=np.vstack([v]), name=out_name, dtype=col_in.dtype, shape=v.shape, unit=col_in.unit)
            table_out.add_column(col_out)

def summarize_sed_results(sed_table):
    """Build a stats summary table for a table that has all the SED results """
    del_cols = ['dnde','dnde_err','dnde_errp','dnde_errn','dnde_ul',
                'e2dnde','e2dnde_err','e2dnde_errp','e2dnde_errn','e2dnde_ul',
                'norm','norm_err','norm_errp','norm_errn','norm_ul',
                'ts']
    stats_cols = ['dnde','dnde_ul',
                  'e2dnde','e2dnde_ul',
                  'norm','norm_ul']

    table_out = Table(sed_table[0])
    table_out.remove_columns(del_cols)
    add_summary_stats_to_table(sed_table, table_out, stats_cols)
    return table_out

def summarize_limits_results(limit_table):
    """Build a stats summary table for a table that has all the SED results """
    del_cols = ['UL_0.68', 'UL_0.95', 'MLES']
    stats_cols = ['UL_0.95', 'MLES']

    table_out = Table(limit_table[0])
    table_out.remove_columns(del_cols)
    add_summary_stats_to_table(limit_table, table_out, stats_cols)
    return table_out


class CollectSEDResults(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    default_options = dict(sed_file=defaults.common['sed_file'],
                           outfile=defaults.generic['outfile'],
                           summaryfile=defaults.generic['summaryfile'],
                           nsims=defaults.sims['nsims'], 
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'])

    collist = [dict(name='e_min', unit = 'MeV'),
               dict(name='e_ref', unit = 'MeV'),
               dict(name='e_max', unit = 'MeV'),
               dict(name='ref_dnde_e_min', unit = 'cm-2 MeV-1 ph s-1'),
               dict(name='ref_dnde_e_max', unit = 'cm-2 MeV-1 ph s-1'),
               dict(name='ref_dnde', unit = 'cm-2 MeV-1 ph s-1'),
               dict(name='ref_flux', unit = 'cm-2 ph s-1'),
               dict(name='ref_eflux', unit = 'cm-2 MeV s-1'),
               dict(name='ref_npred'),
               dict(name='dnde', unit = 'cm-2 MeV-1 ph s-1'),
               dict(name='dnde_err', unit = 'cm-2 MeV-1 ph s-1'),
               dict(name='dnde_errp', unit = 'cm-2 MeV-1 ph s-1'),
               dict(name='dnde_errn', unit = 'cm-2 MeV-1 ph s-1'),
               dict(name='dnde_ul', unit = 'cm-2 MeV-1 ph s-1'),
               dict(name='e2dnde', unit = 'cm-2 MeV s-1'),
               dict(name='e2dnde_err', unit = 'cm-2 MeV s-1'),
               dict(name='e2dnde_errp', unit = 'cm-2 MeV s-1'),
               dict(name='e2dnde_errn', unit = 'cm-2 MeV s-1'),
               dict(name='e2dnde_ul', unit = 'cm-2 MeV s-1'),
               dict(name='norm'),
               dict(name='norm_err'),
               dict(name='norm_errp'),
               dict(name='norm_errn'),
               dict(name='norm_ul'),
               dict(name='ts')]

    def __init__(self, **kwargs):
        """C'tor
        """
        parser = argparse.ArgumentParser(usage="dmpipe-collect-sed [options]",
                                         description="Collect SED results from simulations")
        Link.__init__(self, kwargs.pop('linkname', 'collect-sed'),
                      parser=parser,
                      appname='dmpipe-collect-sed',
                      options=CollectSEDResults.default_options.copy(),
                      **kwargs)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        sedfile = args.sed_file
        first = args.seed
        last = first + args.nsims
        flist = [ sedfile.replace("_SEED.fits", "_%06i.fits"%seed) for seed in range(first, last) ]
        outfile = args.outfile
        summaryfile = args.summaryfile

        outtable = fill_output_table(flist, "SED", CollectSEDResults.collist, nbins=12)

        if outfile not in [None, 'none', 'None']:
            outtable.write(outfile)

        if summaryfile not in [None, 'none', 'None']:
            summary = summarize_sed_results(outtable)
            summary.write(summaryfile)


class CollectLimits(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
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
        parser = argparse.ArgumentParser(usage="dmpipe-collect-limits [options]",
                                         description="Collect Limits from simulations")
        Link.__init__(self, kwargs.pop('linkname', 'collect-limits'),
                      parser=parser,
                      appname='dmpipe-collect-limits',
                      options=CollectLimits.default_options.copy(),
                      **kwargs)

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

        if outfile not in [None, 'none', 'None']:
            fits_utils.write_tables_to_fits(outfile, out_tables, namelist=out_names)

        if summaryfile not in [None, 'none', 'None']:
            summary_tables = [summarize_limits_results(ot) for ot in out_tables[0:-1]]
            summary_tables.append( Table(out_tables[-1][0])  )
            fits_utils.write_tables_to_fits(summaryfile, summary_tables, namelist=out_names)



class ConfigMaker_CollectSED(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           write_full=defaults.collect['write_full'],
                           write_summary=defaults.collect['write_summary'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_CollectSED.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args, require_sim_name=True)
        if targets_yaml is None:
            return job_configs

        write_full = args['write_full']

        targets = load_yaml(targets_yaml)
        for target_name, profile_list in targets.items():
            for profile in profile_list:
                full_key = "%s:%s:%s" % (target_name, profile, sim)
                name_keys = dict(target_type=ttype, 
                                 target_name=target_name,
                                 sim_name=sim,
                                 profile=profile,
                                 fullpath=True)
                sed_file = NAME_FACTORY.sim_sedfile(**name_keys)
                first = args['seed']
                last = first + args['nsims'] - 1
                outfile = sed_file.replace('_SEED.fits','_collected_%06i_%06i.fits'%(first, last))
                logfile = make_nfs_path(outfile.replace('.fits', '.log'))
                if not write_full:
                    outfile = None
                summaryfile = sed_file.replace('_SEED.fits','_summary_%06i_%06i.fits'%(first, last))
                job_config = dict(sed_file=sed_file,
                                  outfile=outfile,
                                  summaryfile=summaryfile,
                                  logfile=logfile,
                                  nsims=args['nsims'],
                                  seed=args['seed'],
                                  dry_run=args['dry_run'])
                job_configs[full_key] = job_config

        return job_configs



class ConfigMaker_CollectLimits(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           specconifg=defaults.common['specconfig'],
                           jprior=defaults.common['jprior'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           write_full=defaults.collect['write_full'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_CollectLimits.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args, require_sim_name=True)
        if targets_yaml is None:
            return job_configs

        specconfig = NAME_FACTORY.resolve_specconfig(args)

        j_prior = args.get('jprior', 'none')
        if j_prior in [None, 'None', 'none']:
            j_prior = 'None'
        write_full = args.get('write_full', False)

        targets = load_yaml(targets_yaml)
        for target_name, profile_list in targets.items():
            for profile in profile_list:
                full_key = "%s:%s:%s" % (target_name, profile, sim)
                name_keys = dict(target_type=ttype, 
                                 target_name=target_name,
                                 sim_name=sim,
                                 profile=profile,
                                 jprior=j_prior, 
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
                                  jprior=j_prior, 
                                  outfile=outfile,
                                  summaryfile=summaryfile,
                                  logfile=logfile,
                                  nsims=args['nsims'],
                                  seed=args['seed'],
                                  dry_run=args['dry_run'])
                job_configs[full_key] = job_config

        return job_configs


class ConfigMaker_CollectStackedLimits(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           rosterlist=defaults.common['targetlist'],
                           jprior=defaults.common['jprior'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           write_full=defaults.collect['write_full'],
                           write_summary=defaults.collect['write_summary'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_CollectLimits.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (roster_yaml, sim) = NAME_FACTORY.resolve_rosterfile(args, require_sim_name=True)
        if roster_yaml is None:
            return job_configs

        specconfig = NAME_FACTORY.resolve_specconfig(args)

        j_prior = args.get('jprior', 'None')
        if j_prior in [None, 'None', 'none']:
            j_prior = 'None'

        write_full = args['write_full']

        roster_dict = load_yaml(roster_yaml)
        for roster_name in roster_dict.keys():
            full_key = "%s:%s:%s" % (roster_name, sim, j_prior)
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


def create_link_collect_sed(**kwargs):
    """Build and return a `Link` object that can invoke TargetAnalysis"""
    collect_sed = CollectSEDResults(**kwargs)
    return collect_sed


def create_link_collect_limits(**kwargs):
    """Build and return a `Link` object that can invoke SEDAnalysis"""
    collect_limits = CollectLimits(**kwargs)
    return collect_limits


def create_sg_collect_sed(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    collect_sed = CollectSEDResults(**kwargs)
    link = collect_sed

    appname = kwargs.pop('appname', 'dmpipe-collect-sed-sg')

    batch_args = get_lsf_default_args()    
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Collect SED results for lots of simulations"

    config_maker = ConfigMaker_CollectSED(link)
    sg = build_sg_from_link(link, config_maker,
                            interface=batch_interface,
                            usage=usage,
                            description=description,
                            appname=appname,
                            **kwargs)
    return sg


def create_sg_collect_limits(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    collect_limits = CollectLimits(**kwargs)
    link = collect_limits

    appname = kwargs.pop('appname', 'dmpipe-collect-limits-sg')

    batch_args = get_lsf_default_args()    
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Collect limtis for lots of simulations"

    config_maker = ConfigMaker_CollectLimits(link)
    sg = build_sg_from_link(link, config_maker,
                            interface=batch_interface,
                            usage=usage,
                            description=description,
                            appname=appname,
                            **kwargs)
    return sg


def create_sg_collect_stacked_limits(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    collect_limits = CollectLimits(**kwargs)
    link = collect_limits

    appname = kwargs.pop('appname', 'dmpipe-collect-stacked-limits-sg')

    batch_args = get_lsf_default_args()    
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Collect limtis for lots of simulations"

    config_maker = ConfigMaker_CollectStackedLimits(link)
    sg = build_sg_from_link(link, config_maker,
                            interface=batch_interface,
                            usage=usage,
                            description=description,
                            appname=appname,
                            **kwargs)
    return sg



def main_collect_sed_single():
    """ Entry point for analysis of a single ROI """
    collect_sed = CollectSEDResults()
    return collect_sed.run_analysis(sys.argv[1:])

def main_collect_limits_single():
    """ Entry point for analysis of a single ROI """
    collect_limits = CollectLimits()
    return collect_limits.run_analysis(sys.argv[1:])

def main_collect_sed_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_collect_sed()
    lsf_sg(sys.argv)

def main_collect_limits_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_collect_limits()
    lsf_sg(sys.argv)

def main_collect_stacked_limits_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_collect_stacked_limits()
    lsf_sg(sys.argv)


if __name__ == "__main__":
    #tables = main_collect_sed_single()
    tables = main_collect_limits_single()

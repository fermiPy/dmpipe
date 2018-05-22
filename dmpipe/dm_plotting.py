#!/usr/bin/env python
#

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Top level scripts to make castro plot and limits plots in mass / sigmav space
"""
from __future__ import absolute_import, division, print_function

import os
from os.path import splitext
import numpy as np

from astropy.table import Table

from fermipy.utils import init_matplotlib_backend, load_yaml
from fermipy.jobs.utils import is_not_null
from fermipy.jobs.link import Link
from fermipy.jobs.scatter_gather import ScatterGather
from fermipy.jobs.slac_impl import make_nfs_path

from dmpipe.dm_spectral_utils import DMCastroData, DMSpecTable
from dmpipe.dm_plotting_utils import plot_dm_castro
from dmpipe.dm_plotting_utils import plot_dm_spectra_by_mass, plot_dm_spectra_by_channel
from dmpipe.dm_plotting_utils import plot_limits_from_arrays, plot_mc_truth

from dmpipe.name_policy import NameFactory
from dmpipe import defaults

init_matplotlib_backend()
NAME_FACTORY = NameFactory(basedir='.')


def get_ul_bands(table, prefix):
    """ Get the upper limit bands a table

    Parameters
    ----------

    table : `astropy.table.Table`
        Table to get the limits from.

    prefix : str
        Prefix to append to the column names for the limits


    Returns
    -------

    output : dict
        A dictionary with the limits bands

    """
    o = dict(q02=np.squeeze(table["%s_q02" % prefix]),
             q16=np.squeeze(table["%s_q16" % prefix]),
             q84=np.squeeze(table["%s_q84" % prefix]),
             q97=np.squeeze(table["%s_q97" % prefix]),
             median=np.squeeze(table["%s_median" % prefix]))
    return o


class PlotDMSpectra(Link):
    """Small class to plot the DM spectra from pre-computed tables.

    """
    appname = 'dmpipe-plot-dm-spectra'
    linkname_default = 'plot-dm-spectra'
    usage = '%s [options]' % (appname)
    description = "Plot the DM spectra stored in pre-computed tables"

    default_options = dict(infile=defaults.generic['infile'],
                           outfile=defaults.generic['outfile'],
                           chan=defaults.common['chan'],
                           mass=defaults.common['mass'],
                           spec_type=defaults.common['spec_type'])

    __doc__ += Link.construct_docstring(default_options)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        dm_spec_table = DMSpecTable.create_from_fits(args.infile)
        dm_plot_by_mass = plot_dm_spectra_by_mass(
            dm_spec_table, chan=args.chan, spec_type=args.spec_type)
        dm_plot_by_chan = plot_dm_spectra_by_channel(
            dm_spec_table, mass=args.mass, spec_type=args.spec_type)

        if args.outfile:
            dm_plot_by_mass[0].savefig(
                args.outfile.replace(
                    '.png', '_%s.png' %
                    args.chan))
            dm_plot_by_chan[0].savefig(
                args.outfile.replace(
                    '.png', '_%1.FGeV.png' %
                    args.mass))


class PlotLimits(Link):
    """Small class to Plot DM limits on <sigma v> versus mass.

    """
    appname = 'dmpipe-plot-limits'
    linkname_default = 'plot-limits'
    usage = '%s [options]' % (appname)
    description = "Plot DM limits on <sigma v> versus mass"

    default_options = dict(infile=defaults.generic['infile'],
                           outfile=defaults.generic['outfile'],
                           chan=defaults.common['chan'],
                           bands=defaults.collect['bands'],
                           sim=defaults.sims['sim'])

    __doc__ += Link.construct_docstring(default_options)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if is_not_null(args.infile):
            tab_m = Table.read(args.infile, hdu="masses")
            tab_s = Table.read(args.infile, hdu=args.chan)

            xvals = tab_m['masses'][0]
            yvals = tab_s['ul_0.95'][0]
            ldict = dict(limits=(xvals, yvals))
        else:
            ldict = {}

        if is_not_null(args.bands):
            tab_b = Table.read(args.bands, hdu=args.chan)
            tab_bm = Table.read(args.bands, hdu="masses")
            bands = get_ul_bands(tab_b, 'ul_0.95')
            bands['masses'] = tab_bm['masses'][0]
        else:
            bands = None

        if is_not_null(args.sim):
            sim_srcs = load_yaml(args.sim)
            injected_src = sim_srcs.get('injected_source', None)
        else:
            injected_src = None

        xlims = (1e1, 1e4)
        ylims = (1e-28, 1e-22)

        dm_plot = plot_limits_from_arrays(ldict, xlims, ylims, bands)

        if injected_src is not None:
            mc_model = injected_src['source_model']
            plot_mc_truth(dm_plot[1], mc_model)

        if args.outfile:
            dm_plot[0].savefig(args.outfile)
            return None
        return dm_plot


class PlotDM(Link):
    """Small class to plot the likelihood vs <sigma v> and DM particle mass

    """
    appname = 'dmpipe-plot-dm'
    linkname_default = 'plot-dm'
    usage = "%s [options]" % (appname)
    description = "Plot the likelihood vs <sigma v> and DM particle mass"

    default_options = dict(infile=defaults.generic['infile'],
                           outfile=defaults.generic['outfile'],
                           chan=defaults.common['chan'])

    __doc__ += Link.construct_docstring(default_options)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)
        exttype = splitext(args.infile)[-1]
        if exttype in ['.fits']:
            dm_castro = DMCastroData.create_from_fitsfile(args.infile, args.chan)
        elif exttype in ['.yaml']:
            dm_castro = DMCastroData.create_from_yamlfile(args.infile, args.chan)
        else:
            raise ValueError("Can not read file type %s for SED" % extype)

        dm_plot = plot_dm_castro(dm_castro)
        if args.outfile:
            dm_plot[0].savefig(args.outfile)
            return None
        return dm_plot


class PlotLimits_SG(ScatterGather):
    """Small class to generate configurations for `PlotLimits`

    This does a triple nested loop over targets, profiles and j-factor priors
    """
    appname = 'dmpipe-plot-limits-sg'
    usage = "%s [options]" % (appname)
    description = "Make castro plots for set of targets"
    clientclass = PlotLimits

    job_time = 60

    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           channels=defaults.common['channels'],
                           jpriors=defaults.common['jpriors'],
                           dry_run=defaults.common['dry_run'])

    __doc__ += Link.construct_docstring(default_options)

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args)
        if targets_yaml is None:
            return job_configs

        jpriors = args['jpriors']
        channels = args['channels']

        base_config = dict(bands=None,
                           sim=sim)

        targets = load_yaml(targets_yaml)
        for target_name, target_list in targets.items():
            for targ_prof in target_list:
                for jprior in jpriors:
                    name_keys = dict(target_type=ttype,
                                     target_name=target_name,
                                     profile=targ_prof,
                                     jprior=jprior,
                                     fullpath=True)
                    input_path = NAME_FACTORY.dmlimitsfile(**name_keys)
                    for chan in channels:
                        targ_key = "%s:%s:%s:%s" % (
                            target_name, targ_prof, jprior, chan)

                        output_path = input_path.replace(
                            '.fits', '_%s.png' % chan)
                        logfile = make_nfs_path(
                            output_path.replace('.png', '.log'))
                        job_config = base_config.copy()
                        job_config.update(dict(infile=input_path,
                                               outfile=output_path,
                                               jprior=jprior,
                                               logfile=logfile,
                                               chan=chan))
                        job_configs[targ_key] = job_config

        return job_configs


class PlotStackedLimits_SG(ScatterGather):
    """Small class to generate configurations for `PlotStackedLimits`

    This does a double nested loop over rosters and j-factor priors
    """
    appname = 'dmpipe-plot-stacked-limits-sg'
    usage = "%s [options]" % (appname)
    description = "Make castro plots for set of targets"
    clientclass = PlotLimits

    job_time = 60

    default_options = dict(ttype=defaults.common['ttype'],
                           rosterlist=defaults.common['targetlist'],
                           bands=defaults.collect['bands'],
                           channels=defaults.common['channels'],
                           jpriors=defaults.common['jpriors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'])

    __doc__ += Link.construct_docstring(default_options)

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
                    targ_key = "%s:%s:%s" % (roster_name, jprior, chan)
                    if sim is not None:
                        seedlist = range(
                            args['seed'], args['seed'] + args['nsims'])
                        sim_path = os.path.join('config', 'sim_%s.yaml' % sim)
                    else:
                        seedlist = [None]
                        sim_path = None

                    for seed in seedlist:
                        if seed is not None:
                            name_keys['seed'] = "%06i" % seed
                            input_path = NAME_FACTORY.sim_stackedlimitsfile(
                                **name_keys)
                            full_targ_key = "%s_%06i" % (targ_key, seed)
                        else:
                            input_path = NAME_FACTORY.stackedlimitsfile(
                                **name_keys)
                            full_targ_key = targ_key

                        output_path = input_path.replace(
                            '.fits', '_%s.png' % chan)
                        logfile = make_nfs_path(
                            output_path.replace('.png', '.log'))
                        job_config = dict(infile=input_path,
                                          outfile=output_path,
                                          jprior=jprior,
                                          logfile=logfile,
                                          sim=sim_path,
                                          chan=chan)
                        job_configs[full_targ_key] = job_config

        return job_configs


class PlotDM_SG(ScatterGather):
    """Small class to generate configurations for `PlotDM`

    This does a quadruple nested loop over targets, profiles,
    j-factor priors and channels
    """
    appname = 'dmpipe-plot-dm-sg'
    usage = "%s [options]" % (appname)
    description = "Make castro plots for set of targets"
    clientclass = PlotDM

    job_time = 60

    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           channels=defaults.common['channels'],
                           jpriors=defaults.common['jpriors'],
                           dry_run=defaults.common['dry_run'])

    __doc__ += Link.construct_docstring(default_options)

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
                        targ_key = "%s:%s:%s:%s" % (
                            target_name, targ_prof, jprior, chan)
                        output_path = input_path.replace(
                            '.fits', '_%s.png' % chan)
                        logfile = make_nfs_path(
                            output_path.replace('.png', '.log'))
                        job_config = dict(infile=input_path,
                                          outfile=output_path,
                                          jprior=jprior,
                                          logfile=logfile,
                                          chan=chan)
                        job_configs[targ_key] = job_config

        return job_configs


class PlotStackedDM_SG(ScatterGather):
    """Small class to generate configurations for `PlotDM`

    This does a triple loop over rosters, j-factor priors and channels
    """
    appname = 'dmpipe-plot-stacked-dm-sg'
    usage = "%s [options]" % (appname)
    description = "Make castro plots for set of targets"
    clientclass = PlotDM

    job_time = 60

    default_options = dict(ttype=defaults.common['ttype'],
                           rosterlist=defaults.common['targetlist'],
                           channels=defaults.common['channels'],
                           jpriors=defaults.common['jpriors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'])

    __doc__ += Link.construct_docstring(default_options)

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
                    targ_key = "%s:%s:%s" % (roster_name, jprior, chan)

                    if sim is not None:
                        seedlist = range(
                            args['seed'], args['seed'] + args['nsims'])
                    else:
                        seedlist = [None]

                    for seed in seedlist:
                        if seed is not None:
                            name_keys['seed'] = "%06i" % seed
                            input_path = NAME_FACTORY.sim_resultsfile(
                                **name_keys)
                            full_targ_key = "%s_%06i" % (targ_key, seed)
                        else:
                            input_path = NAME_FACTORY.resultsfile(**name_keys)
                            full_targ_key = targ_key

                        output_path = input_path.replace(
                            '.fits', '_%s.png' % chan)
                        logfile = make_nfs_path(
                            output_path.replace('.png', '.log'))
                        job_config = dict(infile=input_path,
                                          outfile=output_path,
                                          jprior=jprior,
                                          logfile=logfile,
                                          chan=chan)
                        job_configs[full_targ_key] = job_config

        return job_configs


def register_classes():
    """Register these classes with the `LinkFactory` """
    PlotDMSpectra.register_class()
    PlotLimits.register_class()
    PlotLimits_SG.register_class()
    PlotDM.register_class()
    PlotDM_SG.register_class()
    PlotStackedDM_SG.register_class()
    PlotStackedLimits_SG.register_class()

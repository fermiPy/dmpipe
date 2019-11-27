#!/usr/bin/env python
#

"""
Interface to Dark Matter spectra
"""
from __future__ import absolute_import, division, print_function

import sys
import os

import numpy as np

from astropy.table import Table

from fermipy import fits_utils
from fermipy.castro import CastroData

from fermipy.utils import load_yaml
from fermipy.jobs.utils import is_null, is_not_null
from fermipy.jobs.link import Link
from fermipy.jobs.scatter_gather import ScatterGather
from fermipy.jobs.slac_impl import make_nfs_path

from fermipy.spectrum import DMFitFunction

from dmpipe.dm_spectral_utils import DMSpecTable, DMCastroData
from dmpipe.name_policy import NameFactory
from dmpipe import defaults

NAME_FACTORY = NameFactory(basedir='.')

REF_J = 1.0e19
REF_SIGV = 1.0e-26


class ConvertCastro(Link):
    """Small class to convert SED to DM space.

    """

    appname = 'dmpipe-convert-castro'
    linkname_default = 'convert-castro'
    usage = '%s [options]' % (appname)
    description = "Convert SED to DMCastroData"

    default_options = dict(specfile=defaults.common['specfile'],
                           sed_file=defaults.common['sed_file'],
                           astro_value_file=defaults.common['astro_value_file'],
                           astro_prior=defaults.common['astro_prior'],
                           outfile=defaults.generic['outfile'],
                           limitfile=defaults.generic['limitfile'],
                           # Note that this defaults to -1
                           nsims=defaults.common['nsims'],
                           seed=defaults.sims['seed'],
                           clobber=defaults.common['clobber'])

    __doc__ += Link.construct_docstring(default_options)

    @staticmethod
    def is_decay_sed(sedfile):
        tokens = os.path.splitext(os.path.basename(sedfile))[0].split('_')
        return tokens[2] in ['point', 'dmap', 'dradial', 'hpxdmap']
        
    @staticmethod
    def is_ann_sed(sedfile):
        tokens = os.path.splitext(os.path.basename(sedfile))[0].split('_')
        return tokens[2] in ['point', 'map', 'radial', 'hpxmap']

    @staticmethod
    def select_channels(channels, sedfile):
        sed_ok_decay = ConvertCastro.is_decay_sed(sedfile)
        sed_ok_ann = ConvertCastro.is_ann_sed(sedfile)
        ochans = []
        for chan in channels:
            chan_is_decay = chan.find('_decay') >= 0
            if chan_is_decay:
                if sed_ok_decay:
                    ochans.append(chan)
            else:
                if sed_ok_ann:
                    ochans.append(chan)
        return ochans
    

    @staticmethod
    def convert_sed_to_dm(spec_table, sed, channels, **kwargs):
        """ Convert an SED file to a DMCastroData object

        Parameters
        ----------

        spec_table : `DMSpecTable`
            Object with all the DM spectra

        sed : `CastroData`
            Object with the SED data

        channels : list
            List of the channels to convert

        
        Keyword arguments
        -----------------

        norm_type : str
            Normalization type to use

        j_val : dict
            Dictionary with information about the J-factor

        d_val : dict
            Dictionary with information about the D-factor


        Returns
        -------

        castro_list : list
            List of the `DMCastroData` objects with the Likelihood data

        table_list : list
            List of `astropy.table.Table` objects with the Likelihood data

        name_list : list
            List of names

        """
        c_list = []
        t_list = []
        n_list = []

        mass_table = None
        norm_type = kwargs.get('norm_type', 'eflux')

        j_val = kwargs.get('j_val', None)
        d_val = kwargs.get('d_val', None)

        for chan in channels:
            chan_idx = DMFitFunction.channel_rev_map[chan]
            if chan_idx >= 100:
                astro_val = d_val
            else:
                astro_val = j_val

            dm_castro = spec_table.convert_castro_data(sed, chan_idx, norm_type, astro_val)
            tab_castro = dm_castro.build_scandata_table()

            if mass_table is None:
                mass_table = dm_castro.build_mass_table()
            # except IndexError, msg:
            #    raise IndexError("Skipping channel %s" % msg)
            #    continue
            c_list.append(dm_castro)
            t_list.append(tab_castro)
            n_list.append(chan)

        t_list.append(mass_table)
        n_list.append("MASSES")
        return c_list, t_list, n_list

    @staticmethod
    def extract_dm_limits(dm_castro_list, channels, alphas, mass_table):
        """Extract limits from a series of `DMCastroData` objects
        for a set of channels and masses

        Parameters
        ----------

        dm_castro_lsit : list
            `DMCastroData` objects with all the DM spectra

        channels : list
            List of the channels to convert

        alphas : list
            List of the confidence level threshold to extract limits

        mass_table : `astropy.table.Table`
            Table with the masses.  This just gets appended to the lists of output tables.

        Returns
        -------

        castro_list : list
            List of the `DMCastroData` objects with the Likelihood data

        table_list : list
            List of `astropy.table.Table` objects with the Likelihood data

        name_list : list
            List of names

        """
        l_list = []
        t_list = []
        n_list = []

        for castro_data, chan in zip(dm_castro_list, channels):
            mles = castro_data.mles()
            limit_dict = dict(mles=mles)
            for alpha in alphas:
                if castro_data.decay:
                    limits = castro_data.getLimits(1.-alpha)
                    limit_dict['ll_%.02f' % alpha] = limits
                else:
                    limits = castro_data.getLimits(alpha)
                    limit_dict['ul_%.02f' % alpha] = limits

            tab_limits = castro_data.build_limits_table(limit_dict)
            l_list.append(limit_dict)
            t_list.append(tab_limits)
            n_list.append(chan)

        t_list.append(mass_table)
        n_list.append("MASSES")
        return l_list, t_list, n_list


    @staticmethod
    def convert_sed(spec_table, channels, sed_file, outfile, limitfile, **kwargs):
        """Convert a single SED to DM space.

        Parameters
        ----------

        spec_table : `DMSpecTable`
            Object with all the DM spectra

        channels : list
            List of the channels to convert

        sed_file : str
            Path to the SED file
 
        outfile : str
            Path to write the output `DMCastroData` object to

        limitfile : str
            Path to write the output limits to.


        Keyword arguments
        -----------------
 
        norm_type : str
            Normalization type to use

        j_factor : dict
            Dictionary with information about the J-factor

        d_factor : dict
            Dictionary with information about the J-factor

        clobber : bool
            Flag to overwrite existing files.



        """
        norm_type = kwargs.get('norm_type', 'eflux')
        j_factor = kwargs.get('j_factor', None)
        d_factor = kwargs.get('d_factor', None)
        clobber = kwargs.get('clobber', False)        

        use_chans = ConvertCastro.select_channels(channels, sed_file)

        exttype = os.path.splitext(sed_file)[-1]
        if exttype in ['.fits', '.npy']:
            sed = CastroData.create_from_sedfile(sed_file, norm_type)
        elif exttype in ['.yaml']:
            sed = CastroData.create_from_yamlfile(sed_file)
        else:
            raise ValueError("Can not read file type %s for SED" % extype)
        
        c_list, t_list, n_list = ConvertCastro.convert_sed_to_dm(spec_table, sed, use_chans, 
                                                                 norm_type=norm_type, 
                                                                 j_val=j_factor, d_val=d_factor)

        if is_not_null(outfile):
            fits_utils.write_tables_to_fits(outfile, t_list, clobber=clobber, namelist=n_list)

        if is_not_null(limitfile):
            mass_table = t_list[-1]
            c_list_lim, t_list_lim, n_list_lim = ConvertCastro.extract_dm_limits(
                c_list, use_chans, [0.68, 0.95], mass_table)
            fits_utils.write_tables_to_fits(limitfile, t_list_lim,
                                            clobber=clobber, namelist=n_list_lim)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        norm_type = 'eflux'

        spec_table = DMSpecTable.create_from_fits(args.specfile)
        profile = load_yaml(args.astro_value_file)

        channels = spec_table.channel_names

        j_factor = None
        d_factor = None

        j_value = profile.get('j_integ', None)
        j_sigma = profile.get('j_sigma', None)
 
        if is_null(j_value):
            j_factor = None
        elif is_null(args.astro_prior) or is_null(j_sigma) or j_sigma == 0.0:
            j_factor = j_value
        else:
            j_factor = dict(functype=args.astro_prior,
                            j_value=j_value,
                            mu=j_value, sigma=j_sigma)

        d_value = profile.get('d_integ', None)
        d_sigma = profile.get('d_sigma', None)        

        if is_null(d_value):
            d_factor = None
        elif is_null(args.astro_prior) or is_null(d_sigma) or d_sigma == 0.0:
            d_factor = d_value
        else:
            d_factor = dict(functype=args.astro_prior,
                            d_value=d_value,
                            mu=d_value, sigma=d_sigma)

        if args.nsims < 0:
            seedlist = [None]
        else:
            seedlist = range(args.seed, args.seed + args.nsims)

        for seed in seedlist:
            sedfile = args.sed_file
            outfile = args.outfile
            limitfile = args.limitfile
            if seed is not None:
                sedfile = sedfile.replace('_SEED.fits', '_%06i.fits' % seed)
                if is_not_null(outfile):
                    outfile = outfile.replace(
                        '_SEED.fits',
                        '_%06i.fits' %
                        seed)
                if is_not_null(limitfile):
                    limitfile = limitfile.replace(
                        '_SEED.fits', '_%06i.fits' % seed)

            self.convert_sed(spec_table, channels, sedfile, outfile, limitfile,
                             norm_type=norm_type, j_factor=j_factor,
                             d_factor=d_factor, clober=args.clobber)


class SpecTable(Link):
    """Small class to build a table with all the DM spectra for this analysis

    """
    appname = 'dmpipe-spec-table'
    linkname_default = 'spec-table'
    usage = '%s [options]' % (appname)
    description = "Build a table with the spectra for DM signals"

    default_options = dict(ttype=defaults.common['ttype'],
                           config=defaults.common['config'],
                           specconfig=defaults.common['specconfig'],
                           specfile=defaults.common['specfile'],
                           clobber=defaults.common['clobber'])

    __doc__ += Link.construct_docstring(default_options)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if args.ttype is not None:
            name_keys = dict(target_type=args.ttype,
                             fullpath=True)
            config_file = NAME_FACTORY.ttypeconfig(**name_keys)
            spec_config = NAME_FACTORY.specconfig(**name_keys)
            spec_file = NAME_FACTORY.specfile(**name_keys)
        else:
            config_file = None
            spec_config = None
            spec_file = None

        if is_not_null(args.config):
            config_file = args.config
        if is_not_null(args.specconfig):
            spec_config = args.specconfig
        if is_not_null(args.specfile):
            spec_file = args.specfile

        if config_file is None:
            sys.stderr.write('No input configuration file is specified')
            return -1

        if spec_config is None:
            sys.stderr.write('No input spectra configurate file is specified')
            return -1

        if spec_file is None:
            sys.stderr.write('No output spectra file is specified')
            return -1

        spec_config = load_yaml(spec_config)
        channels = spec_config['channels']
        
        if isinstance(spec_config['masses'], dict):
            masses = np.logspace(np.log10(spec_config['masses']['mass_min']),
                                 np.log10(spec_config['masses']['mass_max']),
                                 spec_config['masses']['mass_nstep'])
        elif isinstance(spec_config['masses'], list):
            masses = spec_config['masses']

        dm_spec_table = DMSpecTable.create_from_config(config_file, channels, masses)
        dm_spec_table.write_fits(spec_file, args.clobber)
        return 0


class StackLikelihood(Link):
    """Small class to stack likelihoods that were written to `DMCastroData` objects.

    """
    appname = 'dmpipe-stack-likelihood'
    linkname_default = 'stack-likelihood'
    usage = '%s [options]' % (appname)
    description = "Stack the likelihood from a set of targets"

    default_options = dict(ttype=defaults.common['ttype'],
                           specconfig=defaults.common['specconfig'],
                           rosterlist=defaults.common['rosterlist'],
                           astro_prior=defaults.common['astro_prior'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           clobber=defaults.common['clobber'])

    __doc__ += Link.construct_docstring(default_options)

    @staticmethod
    def is_decay_roster(roster):
        tokens = roster.split('_')
        return tokens[-1] in ['point', 'dmap', 'dradial']
        
    @staticmethod
    def is_ann_roster(roster):
        tokens = roster.split('_')
        return tokens[-1] in ['point', 'map', 'radial']

    @staticmethod
    def select_channels(channels, roster):
        roster_ok_decay = StackLikelihood.is_decay_roster(roster)
        roster_ok_ann = StackLikelihood.is_ann_roster(roster)
        ochans = []
        for chan in channels:
            chan_is_decay = chan.find('_decay') >= 0
            if chan_is_decay:
                if roster_ok_decay:
                    ochans.append(chan)
            else:
                if roster_ok_ann:
                    ochans.append(chan)
        return ochans
    

    @staticmethod
    def stack_roster(rost, ttype,
                     channels, astro_prior_key, sim, seed):
        """ Stack all of the DMCastroData in a roster

        Parameters
        ----------

        rost : list
            List of the targets

        ttype : str
            Type of target, used for bookkeeping and file names

        channels : list
            List of the channels to convert

        j_prior_key : str
            String that identifies the type of prior on the J-factor

        sim : str
            String that specifies the simulation scenario

        seed : int or None
            Key for the simulation instance, used for bookkeeping and file names

        Returns
        -------

        output : dict
            Dictionary of `DMCastroData` objects, keyed by channel

        """
        component_dict = {}
        out_dict = {}
        for chan in channels:
            component_dict[chan] = []

        for target_key in rost:
            tokens = target_key.split(':')
            name_keys = dict(target_type=ttype,
                             target_name=tokens[0],
                             profile=tokens[1],
                             fullpath=True,
                             sim_name=sim,
                             seed="%06i" % seed,
                             astro_prior=astro_prior_key)

            if is_not_null(sim):
                dmlike_path = NAME_FACTORY.sim_dmlikefile(**name_keys)
            else:
                dmlike_path = NAME_FACTORY.dmlikefile(**name_keys)

            tab_m = Table.read(dmlike_path, hdu="MASSES")

            for chan in channels:
                if chan.find('_decay') >= 0:
                    decay = True
                    norm_type = 'tau'
                else:
                    decay = False
                    norm_type = 'sigmav'
                try:
                    tab_s = Table.read(dmlike_path, hdu=chan)
                except KeyError:
                    continue
                dm_castro = DMCastroData.create_from_tables(tab_s, tab_m, norm_type, decay)
                component_dict[chan].append(dm_castro)

        for chan, comps in component_dict.items():
            if not comps:
                continue
            decay = chan.find('_decay') >= 0
            stacked = DMCastroData.create_from_stack(comps, decay=decay)
            out_dict[chan] = stacked

        return out_dict

    @staticmethod
    def write_fits_files(stacked_dict, resultsfile, limitfile, clobber=False):
        """ Write the stacked DMCastroData object and limits a FITS files

        Parameters
        ----------

        stacked_dict : dict
            Dictionary of `DMCastroData` objects, keyed by channel

        resultsfile : str
            Path to the output file to write the `DMCastroData` objects to

        limitfile : str
            Path to write the upper limits to

        clobber : bool
            Overwrite existing files

        """
        channels = stacked_dict.keys()
        t_list = []
        n_list = []
        lim_list = []
        lim_table_list = []
        mass_table = None
        alphas = [0.68, 0.95]
        for chan in channels:
            is_decay = chan.find('_decay') >= 0
            stacked = stacked_dict[chan]
            mles = stacked.mles()
            limit_dict = dict(mles=mles)
            for alpha in alphas:
                if is_decay:
                    limits = stacked.getLimits(1. - alpha)
                    limit_dict['ll_%.02f' % alpha] = limits
                else:
                    limits = stacked.getLimits(alpha)
                    limit_dict['ul_%.02f' % alpha] = limits
            tab_limits = stacked.build_limits_table(limit_dict)
            if mass_table is None:
                mass_table = stacked.build_mass_table()
            t_list.append(stacked.build_scandata_table())
            n_list.append(chan)
            lim_list.append(limit_dict)
            lim_table_list.append(tab_limits)

        t_list.append(mass_table)
        lim_table_list.append(mass_table)
        n_list.append("MASSES")
        fits_utils.write_tables_to_fits(resultsfile, t_list,
                                        clobber=clobber, namelist=n_list)
        fits_utils.write_tables_to_fits(limitfile, lim_table_list,
                                        clobber=clobber, namelist=n_list)
       


    @staticmethod
    def write_stacked(ttype, roster_name, stacked_dict,
                      astro_prior_key, sim, seed, clobber):
        """ Write the stacked DMCastroData object to a FITS file

        Parameters
        ----------

        ttype : str
            Type of target, used for bookkeeping and file names

        roster_name : str
            Name of the roster, used for bookkeeping and file names

        stacked_dict : dict
            Dictionary of `DMCastroData` objects, keyed by channel

        astro_prior_key : str
            String that identifies the type of prior on the J-factor

        sim : str
            String that specifies the simulation scenario

        seed : int or None
            Key for the simulation instance, used for bookkeeping and file names

        clobber : bool
            Flag to overwrite existing files.

        """
        name_keys = dict(target_type=ttype,
                         target_name="stacked",
                         fullpath=True,
                         roster_name=roster_name,
                         sim_name=sim,
                         seed="%06i" % seed,
                         astro_prior=astro_prior_key)

        if is_not_null(sim):
            outdir = NAME_FACTORY.sim_targetdir(**name_keys)
            outpath = NAME_FACTORY.sim_resultsfile(**name_keys)
        else:
            outdir = NAME_FACTORY.targetdir(**name_keys)
            outpath = NAME_FACTORY.resultsfile(**name_keys)

        try:
            os.makedirs(outdir)
        except OSError:
            pass

        limitfile = outpath.replace('results', 'limits')
        print("Writing stacked results %s" % outpath)
        StackLikelihood.write_fits_files(stacked_dict, outpath, limitfile, clobber)
        

    @staticmethod
    def stack_rosters(roster_dict, ttype, channels,
                      astro_prior_key, sim, seed, clobber):
        """ Stack all of the DMCastroData in a dictionary of rosters

        Parameters
        ----------

        roster_dict : dict
            Dictionary of all the roster being used.

        ttype : str
            Type of target, used for bookkeeping and file names

        channels : list
            List of the channels to convert

        astro_prior_key : str
            String that identifies the type of prior on the J-factor

        sim : str
            String that specifies the simulation scenario

        seed : int or None
            Key for the simulation instance, used for bookkeeping and file names

        clobber : bool
            Flag to overwrite existing files.

        """
        for roster_name, rost in roster_dict.items():
            rost_chans = StackLikelihood.select_channels(channels, roster_name)
            stacked_dict = StackLikelihood.stack_roster(rost, ttype,
                                                        rost_chans, astro_prior_key, sim, seed)
            StackLikelihood.write_stacked(ttype, roster_name, stacked_dict,
                                          astro_prior_key, sim, seed, clobber)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if args.ttype is None:
            raise RuntimeError('Target type must be specified')

        name_keys = dict(target_type=args.ttype,
                         rosterlist='roster_list.yaml',
                         sim_name=args.sim,
                         fullpath=True)

        spec_config = NAME_FACTORY.specconfig(**name_keys)
        if is_not_null(args.specconfig):
            spec_config = args.specconfig

        spec_config = load_yaml(spec_config)
        channels = spec_config['channels']

        if is_not_null(args.sim):
            roster_file = NAME_FACTORY.sim_rosterfile(**name_keys)
            sim_name = args.sim
            is_sim = True
        else:
            roster_file = NAME_FACTORY.rosterfile(**name_keys)
            is_sim = False
            sim_name = None

        if is_not_null(args.rosterlist):
            roster_file = args.rosterlist

        roster_dict = load_yaml(roster_file)

        if is_sim:
            seedlist = range(args.seed, args.seed + args.nsims)
        else:
            seedlist = [0]

        astro_prior = args.astro_prior
        if is_null(astro_prior):
            astro_prior = 'none'

        for seed in seedlist:
            StackLikelihood.stack_rosters(roster_dict, args.ttype, channels,
                                          astro_prior, sim_name, seed, args.clobber)


class ConvertCastro_SG(ScatterGather):
    """Small class to generate configurations for the `ConvertCastro` script

    This does a triple loop over targets, spatial profiles and J-factor priors
    """
    appname = 'dmpipe-convert-castro-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = ConvertCastro

    job_time = 600

    default_options = dict(ttype=defaults.common['ttype'],
                           specfile=defaults.common['specfile'],
                           targetlist=defaults.common['targetlist'],
                           astro_priors=defaults.common['astro_priors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           clobber=defaults.common['clobber'])

    __doc__ += Link.construct_docstring(default_options)

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args)
        if targets_yaml is None:
            return job_configs

        specfile = NAME_FACTORY.resolve_specfile(args)

        targets = load_yaml(targets_yaml)

        astro_priors = args['astro_priors']
        clobber = args['clobber']

        if is_not_null(sim):
            is_sim = True
            nsims = args['nsims']
            seed = args['seed']
        else:
            is_sim = False
            nsims = -1
            seed = -1

        base_config = dict(specfile=specfile,
                           nsims=nsims,
                           seed=seed,
                           clobber=clobber)

        for target_name, profile_list in targets.items():
            for profile in profile_list:
                for astro_prior in astro_priors:
                    full_key = "%s:%s:%s" % (target_name, profile, astro_prior)
                    target_version = profile.split('_')[0]
                    name_keys = dict(target_type=ttype,
                                     target_name=target_name,
                                     target_version=target_version,
                                     profile=profile,
                                     astro_prior=astro_prior,
                                     fullpath=True)
                    if is_sim:
                        name_keys['sim_name'] = sim
                        sed_file = NAME_FACTORY.sim_sedfile(**name_keys)
                        astro_value_yaml = NAME_FACTORY.sim_astro_valuefile(**name_keys)
                        outfile = NAME_FACTORY.sim_dmlikefile(**name_keys)
                        limitfile = NAME_FACTORY.sim_dmlimitsfile(**name_keys)
                        full_key += ":%s" % sim
                    else:
                        sed_file = NAME_FACTORY.sedfile(**name_keys)
                        astro_value_yaml = NAME_FACTORY.astro_valuefile(**name_keys)
                        outfile = NAME_FACTORY.dmlikefile(**name_keys)
                        limitfile = NAME_FACTORY.dmlimitsfile(**name_keys)

                    logfile = make_nfs_path(outfile.replace('.fits', '.log'))
                    job_config = base_config.copy()
                    job_config.update(dict(sed_file=sed_file,
                                           astro_value_file=astro_value_yaml,
                                           astro_prior=astro_prior,
                                           outfile=outfile,
                                           limitfile=limitfile,
                                           logfile=logfile))

                    job_configs[full_key] = job_config

        return job_configs


class StackLikelihood_SG(ScatterGather):
    """Small class to generate configurations for `StackLikelihood`

    This loops over the types of priors on the J-factor
    """
    appname = 'dmpipe-stack-likelihood-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = StackLikelihood

    job_time = 120

    default_options = dict(ttype=defaults.common['ttype'],
                           specconfig=defaults.common['specfile'],
                           rosterlist=defaults.common['rosterlist'],
                           astro_priors=defaults.common['astro_priors'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           clobber=defaults.common['clobber'])

    __doc__ += Link.construct_docstring(default_options)

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        astro_priors = args['astro_priors']
        clobber = args['clobber']
        sim = args['sim']

        if is_not_null(sim):
            is_sim = True
            nsims = args['nsims']
            seed = args['seed']
        else:
            is_sim = False
            nsims = -1
            seed = -1

        base_config = dict(ttype=args['ttype'],
                           specconfig=args['specconfig'],
                           rosterlist=args['rosterlist'],
                           sim=sim,
                           nsims=nsims,
                           seed=seed,
                           clobber=clobber)

        for astro_prior in astro_priors:

            name_keys = dict(target_type=args['ttype'],
                             target_name='stacked',
                             astro_prior=astro_prior,
                             sim_name=sim,
                             fullpath=True)
            if is_sim:
                target_dir = NAME_FACTORY.sim_targetdir(**name_keys)
                full_key = "%s:%s" % (astro_prior, sim)
            else:
                target_dir = NAME_FACTORY.targetdir(**name_keys)
                full_key = astro_prior

            logfile = os.path.join(target_dir, 'stack_%s.log' % astro_prior)

            job_config = base_config.copy()
            job_config.update(dict(astro_prior=astro_prior,
                                   logfile=logfile))

            job_configs[full_key] = job_config

        return job_configs


def register_classes():
    """Register these classes with the `LinkFactory` """
    ConvertCastro.register_class()
    ConvertCastro_SG.register_class()
    StackLikelihood.register_class()
    StackLikelihood_SG.register_class()
    SpecTable.register_class()

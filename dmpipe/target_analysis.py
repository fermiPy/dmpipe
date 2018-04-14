#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Run gtsrcmaps for a single energy plane for a single source

This is useful to parallize the production of the source maps
"""
from __future__ import absolute_import, division, print_function

import os
import sys
import argparse
import numpy as np

from shutil import copyfile

from dmsky.roster import RosterLibrary

from fermipy.utils import load_yaml, write_yaml, init_matplotlib_backend

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro

from fermipy.jobs.chain import Link
from fermipy.jobs.scatter_gather import ConfigMaker, build_sg_from_link
from fermipy.jobs.lsf_impl import make_nfs_path, get_lsf_default_args, LSF_Interface

from dmpipe.name_policy import NameFactory
from dmpipe import defaults

init_matplotlib_backend('Agg')

try:
    from fermipy.gtanalysis import GTAnalysis
    HAVE_ST = True
except ImportError:
    HAVE_ST = False

NAME_FACTORY = NameFactory(basedir=('.'))

class TargetPreparer(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """

    default_options = dict(ttype=defaults.common['ttype'],
                           roster=defaults.common['roster'],
                           config=defaults.common['config'],
                           sim=defaults.sims['sim'],
                           dry_run=defaults.common['dry_run'])

    copyfiles = ['srcmap_00.fits', 'fit_baseline.fits', 'fit_baseline.npy', 'fit_baseline_00.xml']

    def __init__(self, **kwargs):
        """C'tor
        """
        parser = argparse.ArgumentParser(usage="dmpipe-prepare-targets [options]",
                                         description="Prepare directories for target analyses")
        Link.__init__(self, kwargs.pop('linkname', 'prepare-targets'),
                      parser=parser,
                      appname='dmpipe-prepare-targets',
                      options=TargetPreparer.default_options.copy(),
                      **kwargs)

    @staticmethod
    def copy_analysis_files(orig_dir, dest_dir, files):
        """ Copy a list of files from orig_dir to dest_dir"""
        for f in files:
            orig_path = os.path.join(orig_dir, f)
            dest_path = os.path.join(dest_dir, f)
            try:
                copyfile(orig_path, dest_path)
            except IOError:
                sys.stderr.write("WARNING: failed to copy %s\n"%orig_path)


    @staticmethod
    def write_target_dirs(ttype, roster_dict, base_config, sim):
        """ Create and populate directoris for target analysis
        """
        target_dict = {}

        target_info_dict = {}
        roster_info_dict = {}

        if sim is None or sim == 'none':
            is_sim = False
        else:
            is_sim = True

        try:
            if is_sim:
                os.makedirs("%s_sim"%ttype)
            else:
                os.makedirs(ttype)
        except OSError:
            pass

        for roster_name, rost in roster_dict.items():
            tlist = []
            for target_name, target in rost.items():
                if is_sim:
                    target_key = "%s:%s:%s" % (target_name, target.version, sim)
                else:
                    target_key = "%s:%s" % (target_name, target.version)
                print("Writing %s" % (target_key))
                tlist.append(target_key)
                if target_info_dict.has_key(target_name):
                    target_info_dict[target_name].append(target.version)
                else:
                    target_info_dict[target_name] = [target.version]
                name_keys = dict(target_type=ttype,
                                 target_name=target_name,
                                 profile=target.version,
                                 sim_name=sim,
                                 fullpath=True)
                
                target_dir = NAME_FACTORY.targetdir(**name_keys)
                if is_sim:                    
                    profile_path = NAME_FACTORY.sim_profilefile(**name_keys)
                    sim_target_dir = NAME_FACTORY.sim_targetdir(**name_keys)
                    target_config_path = os.path.join(sim_target_dir, 'config.yaml')
                else:
                    profile_path = NAME_FACTORY.profilefile(**name_keys)
                    sim_target_dir = None
                    target_config_path = os.path.join(target_dir, 'config.yaml')

                jmap_path = profile_path.replace('.yaml', 'fits')

                if target_dict.has_key(target_name):
                    # Already made the config for this target
                    target_config = target_dict[target_name]
                else:
                    # Make the config for this target
                    try:
                        if is_sim:
                            os.makedirs(sim_target_dir)
                        else:
                            os.makedirs(target_dir)
                    except OSError:
                        pass
                    target_config = base_config.copy()
                    target_config['selection']['ra'] = target.ra
                    target_config['selection']['dec'] = target.dec
                    if is_sim:
                        TargetPreparer.copy_analysis_files(target_dir, sim_target_dir, TargetPreparer.copyfiles)
                        target_config['gtlike']['bexpmap'] = os.path.abspath(os.path.join(target_dir,'bexpmap_00.fits'))
                        target_config['gtlike']['srcmap'] = os.path.abspath(os.path.join(sim_target_dir,'srcmap_00.fits'))
                        target_config['gtlike']['use_external_srcmap'] = True
                        sim_orig_file = os.path.join('config', 'sim_%s.yaml'%sim)
                        sim_dest_file = os.path.join(sim_target_dir, 'sim_%s.yaml'%sim)
                        copyfile(sim_orig_file, sim_dest_file)

                    target_dict[target_name] = target_config
                    write_yaml(target_config, target_config_path)

                profile_data = target.profile.copy()
                #target.write_jmap_wcs(jmap_path, clobber=True)
                profile_data['j_integ'] = target.j_integ
                profile_data['j_sigma'] = target.j_sigma
                profile_data['j_map_file'] = jmap_path

                write_yaml(profile_data, profile_path)

            roster_info_dict[roster_name] = tlist

        if is_sim:
            roster_file = os.path.join("%s_sim"%ttype, "sim_%s"%sim, 'roster_list.yaml')
            target_file = os.path.join("%s_sim"%ttype, "sim_%s"%sim, 'target_list.yaml')
        else:
            roster_file = os.path.join(ttype, 'roster_list.yaml')
            target_file = os.path.join(ttype, 'target_list.yaml')

        write_yaml(roster_info_dict, roster_file)
        write_yaml(target_info_dict, target_file)


    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)
        roster_lib = RosterLibrary()
        roster_dict = {}

        if args.roster is None or args.roster == 'none':
            sys.stderr.write("You must specify a target roster")
            return -1
        
        if args.ttype is None or args.ttype == 'none':
            sys.stderr.write("You must specify a target type")
            return -1

        name_keys = dict(target_type=args.ttype,
                         fullpath=True)
        config_file = NAME_FACTORY.ttypeconfig(**name_keys)
        if args.config is not None and args.config != 'none':
            config_file = args.config

        rost = roster_lib.create_roster(args.roster)
        roster_dict[args.roster] = rost

        base_config = load_yaml(config_file)
        TargetPreparer.write_target_dirs(args.ttype, roster_dict, base_config, args.sim)


class TargetAnalysis(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    default_options = dict(config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs):
        """C'tor
        """
        parser = argparse.ArgumentParser(usage="dmpipe-analyze-roi [options]",
                                         description="Run analysis of a single ROI")
        Link.__init__(self, kwargs.pop('linkname', 'analyze-roi'),
                      parser=parser,
                      appname='dmpipe-analyze-roi',
                      options=TargetAnalysis.default_options.copy(),
                      **kwargs)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if not HAVE_ST:
            raise RuntimeError("Trying to run fermipy analysis, but don't have ST")
        
        gta = GTAnalysis(args.config, logging={'verbosity': 3},
                         fileio={'workdir_regex': '\.xml$|\.npy$'})

        gta.setup(overwrite=False)
        gta.free_sources(False)
        gta.print_roi()
        gta.optimize()
        gta.print_roi()

        exclude = ['3FGL J1707.8+5626']

        # Localize all point sources
        for src in sorted(gta.roi.sources, key=lambda t: t['ts'], reverse=True):
            #    for s in gta.roi.sources:

            if not src['SpatialModel'] == 'PointSource':
                continue
            if src['offset_roi_edge'] > -0.1:
                continue

            if src.name in exclude:
                continue
            if not '3FGL' in src.name:
                continue

            gta.localize(src.name, nstep=5, dtheta_max=0.5, update=True,
                         prefix='base', make_plots=True)

        gta.optimize()
        gta.print_roi()

        gta.write_roi('base_roi', make_plots=True)

        gta.find_sources(sqrt_ts_threshold=5.0, search_skydir=gta.roi.skydir,
                         search_minmax_radius=[1.0, np.nan])
        gta.optimize()
        gta.print_roi()
        gta.print_params()

        gta.free_sources(skydir=gta.roi.skydir, distance=1.0, pars='norm')
        gta.fit()
        gta.print_roi()
        gta.print_params()

        gta.write_roi('fit_baseline', make_plots=True)


class SEDAnalysis(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    default_options = dict(config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'],
                           profiles=([], 'Profiles to build SED for', list))

    def __init__(self, **kwargs):
        """C'tor
        """
        parser = argparse.ArgumentParser(usage='dmpipe-analyze-sed',
                                         description="Extract the SED for a single target")
        Link.__init__(self, kwargs.pop('linkname', 'analyze-sed'),
                      parser=parser,
                      appname='dmpipe-analyze-sed',
                      options=SEDAnalysis.default_options.copy(),
                      **kwargs)

    @staticmethod
    def _build_profile_dict(basedir, profile_name):
        """
        """
        profile_path = os.path.join(basedir, "profile_%s.yaml" % profile_name)
        #FIXME, use the profile_config to build the profile_dict
        #profile_config = load_yaml(profile_path)

        profile_dict = {}
        profile_dict['SpatialModel'] = 'PointSource'
        profile_dict['SpectrumType'] = 'PowerLaw'
        return profile_name, profile_dict

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if not HAVE_ST:
            raise RuntimeError("Trying to run fermipy analysis, but don't have ST")

        gta = GTAnalysis(args.config,
                         logging={'verbosity': 3},
                         fileio={'workdir_regex': '\.xml$|\.npy$'})
        gta.setup(overwrite=False)
        gta.load_roi('fit_baseline')
        gta.print_roi()

        basedir = os.path.dirname(args.config)
        # This should be a no-op, b/c it was done in the baseline analysis

        gta.free_sources(skydir=gta.roi.skydir, distance=1.0, pars='norm')
        ylims = [1e-8, 1e-5]

        for profile in args.profiles:
            pkey, pdict = SEDAnalysis._build_profile_dict(basedir, profile)
            outfile="sed_%s.fits" % pkey
            outplot="sed_%s.png" % pkey
            # test_case need to be a dict with spectrum and morphology
            gta.add_source(pkey, pdict)
            # refit the ROI
            gta.fit()
            # build the SED
            gta.sed(pkey, outfile=outfile)       
            # plot the SED
            # FIXME, make this optional
            if False:
                castro_data = CastroData.create_from_sedfile(os.path.join(basedir, outfile))
                plot = plotCastro(castro_data, ylims)
                plot[0].savefig(outplot)    
            # remove the source
            gta.delete_source(pkey)
            # put the ROI back to how it was
            gta.load_xml('fit_baseline')

            

        return gta


class ConfigMaker_TargetAnalysis(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_TargetAnalysis.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args)
        if targets_yaml is None:
            return job_configs

        config_yaml = 'config.yaml'
        config_override = args.get('config')
        if config_override is not None and config_override != 'none':
            config_yaml = config_override

        targets = load_yaml(targets_yaml)

        for target_name in targets.keys():
            name_keys = dict(target_type=ttype,
                             target_name=target_name,
                             fullpath=True)
            target_dir = NAME_FACTORY.targetdir(**name_keys)
            config_path = os.path.join(target_dir, config_yaml)
            logfile = make_nfs_path(os.path.join(target_dir, "%s_%s.log"%(self.link.linkname, target_name)))
            job_config = dict(config=config_path, 
                              logfile=logfile)
            job_configs[target_name] = job_config

        return job_configs


class ConfigMaker_SEDAnalysis(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_SEDAnalysis.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}
        
        ttype = args['ttype']
        (targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args)
        if targets_yaml is None:
            return job_configs

        targets = load_yaml(targets_yaml)
        config_yaml = 'config.yaml'
 
        for target_name, target_list in targets.items():
            name_keys = dict(target_type=ttype,
                             target_name=target_name,
                             fullpath=True)
            target_dir = NAME_FACTORY.targetdir(**name_keys)
            config_path = os.path.join(target_dir, config_yaml)
            logfile =  make_nfs_path(os.path.join(target_dir, "%s_%s.log"%(self.link.linkname, target_name)))
            job_config = dict(config=config_path,
                              profiles=target_list,
                              logfile=logfile)
            job_configs[target_name] = job_config

        return job_configs


def create_link_prepare_targets(**kwargs):
    """Build and return a `Link` object that can invoke TargetPreparer"""
    target_prep = TargetPreparer(**kwargs)
    return target_prep


def create_link_roi_analysis(**kwargs):
    """Build and return a `Link` object that can invoke TargetAnalysis"""
    target_analysis = TargetAnalysis(**kwargs)
    return target_analysis


def create_link_sed_analysis(**kwargs):
    """Build and return a `Link` object that can invoke SEDAnalysis"""
    sed_analysis = SEDAnalysis(**kwargs)
    return sed_analysis


def create_sg_roi_analysis(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    roi_analysis = TargetAnalysis(**kwargs)
    link = roi_analysis

    appname = kwargs.pop('appname', 'dmpipe-analyze-roi-sg')

    batch_args = get_lsf_default_args()    
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"

    config_maker = ConfigMaker_TargetAnalysis(link)
    sg = build_sg_from_link(link, config_maker,
                            interface=batch_interface,
                            usage=usage,
                            description=description,
                            appname=appname,
                            **kwargs)
    return sg


def create_sg_sed_analysis(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    sed_analysis = SEDAnalysis(**kwargs)
    link = sed_analysis

    appname = kwargs.pop('appname', 'dmpipe-analyze-sed-sg')

    batch_args = get_lsf_default_args()    
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"

    config_maker = ConfigMaker_SEDAnalysis(link)
    sg = build_sg_from_link(link, config_maker,
                            interface=batch_interface,
                            usage=usage,
                            description=description,
                            appname=appname,
                            **kwargs)
    return sg


def main_prepare_targets():
    """ Entry point for analysis of a single ROI """
    target_prep = TargetPreparer()
    target_prep.run_analysis(sys.argv[1:])


def main_roi_single():
    """ Entry point for analysis of a single ROI """
    target_analysis = TargetAnalysis()
    target_analysis.run_analysis(sys.argv[1:])


def main_sed_single():
    """ Entry point for analysis of a single ROI """
    sed_analysis = SEDAnalysis()
    sed_analysis.run_analysis(sys.argv[1:])


def main_roi_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_roi_analysis()
    lsf_sg(sys.argv)


def main_sed_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_sed_analysis()
    lsf_sg(sys.argv)


if __name__ == "__main__":
    main_sed_batch()

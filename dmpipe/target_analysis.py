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

from dmsky.roster import RosterLibrary

from fermipy.utils import load_yaml, write_yaml, init_matplotlib_backend

from fermipy.jobs.chain import Link
from fermipy.jobs.scatter_gather import ConfigMaker
from fermipy.jobs.lsf_impl import build_sg_from_link

init_matplotlib_backend('Agg')

try:
    from fermipy.gtanalysis import GTAnalysis
    HAVE_ST = True
except ImportError:
    HAVE_ST = False


class TargetPreparer(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """

    default_options = dict(roster=(None, 'Roster to build targets for', str),
                           topdir=(None, 'Top level output directory', str),
                           baseconfig=('config_baseline.yaml', 'Name of config script', str),
                           dry_run=(False, 'Print but do not run commands', bool))

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
    def write_target_dirs(basedir, roster_dict, base_config):
        """ Create and populate directoris for target analysis
        """
        target_dict = {}

        target_info_dict = {}
        roster_info_dict = {}

        try:
            os.makedirs(basedir)
        except OSError:
            pass

        for roster_name, rost in roster_dict.items():
            tlist = []
            for target_name, target in rost.items():
                target_key = "%s:%s" % (target_name, target.version)
                print("Writing %s" % (target_key))
                tlist.append(target_key)
                if target_info_dict.has_key(target_name):
                    target_info_dict[target_name].append(target.version)
                else:
                    target_info_dict[target_name] = [target.version]
                target_dir = os.path.join(basedir, target_name)
                target_config_path = os.path.join(target_dir, 'config_baseline.yaml')
                jmap_path = os.path.join(target_dir, 'profile_%s.fits' % target.version)
                profile_path = os.path.join(target_dir, 'profile_%s.yaml' % target.version)

                if target_dict.has_key(target_name):
                    # Already made the config for this target
                    target_config = target_dict[target_name]
                else:
                    # Make the config for this target
                    try:
                        os.makedirs(target_dir)
                    except OSError:
                        pass
                    target_config = base_config.copy()
                    target_config['selection']['ra'] = target.ra
                    target_config['selection']['dec'] = target.dec
                    target_dict[target_name] = target_config
                    write_yaml(target_config, target_config_path)

                profile_data = target.profile.copy()
                #target.write_jmap_wcs(jmap_path, clobber=True)
                profile_data['j_integ'] = target.j_integ
                profile_data['j_sigma'] = target.j_sigma
                profile_data['j_map_file'] = jmap_path
                write_yaml(profile_data, profile_path)

            roster_info_dict[roster_name] = tlist

        write_yaml(roster_info_dict, os.path.join(basedir, 'roster_list.yaml'))
        write_yaml(target_info_dict, os.path.join(basedir, 'target_list.yaml'))

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)
        roster_lib = RosterLibrary()
        roster_dict = {}
        rost = roster_lib.create_roster(args.roster)
        roster_dict[args.roster] = rost

        base_config = load_yaml(args.baseconfig)

        TargetPreparer.write_target_dirs(args.topdir, roster_dict, base_config)


class TargetAnalysis(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    default_options = dict(config=('config_baseline.yaml', 'Name of config script', str),
                           dry_run=(False, 'Print but do not run commands', bool))

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

        gta.find_sources(sqrt_ts_threshold=5.0)
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
    default_options = dict(config=('config_baseline.yaml', 'Name of config script', str),
                           dry_run=(False, 'Print but do not run commands', bool),
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

        for profile in args.profiles:
            pkey, pdict = SEDAnalysis._build_profile_dict(basedir, profile)
            # test_case need to be a dict with spectrum and morphology
            gta.add_source(pkey, pdict)
            # refit the ROI
            gta.fit()
            # build the SED
            gta.sed(pkey, outfile="sed_%s.fits" % pkey)
            # remove the source
            gta.delete_source(pkey)
            # put the ROI back to how it was
            gta.load_xml('fit_baseline')

        return gta


class ConfigMaker_TargetAnalysis(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(targetlist=('target_list.yaml', 'Yaml file with list of targets', str),
                           config=('config_baseline.yaml', 'Name of configuration file', str),
                           topdir=(None, 'Top level directory', str))

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_TargetAnalysis.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        input_config = {}
        job_configs = {}
        output_config = {}

        topdir = args['topdir']
        targets_yaml = os.path.join(topdir, args['targetlist'])
        config_yaml = args['config']

        logdir = os.path.abspath(topdir).replace('gpfs', 'nfs')

        try:
            targets = load_yaml(targets_yaml)
        except IOError:
            targets = {}

        for target_name in targets.keys():
            config_path = os.path.join(topdir, target_name, config_yaml)
            logfile = os.path.join(logdir, target_name, "%s_%s.log"%(self.link.linkname, target_name))
            job_config = dict(config=config_path, 
                              logfile=logfile)
            job_configs[target_name] = job_config

        return input_config, job_configs, output_config


class ConfigMaker_SEDAnalysis(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(targetlist=('target_list.yaml', 'Yaml file with list of targets', str),
                           config=('config_baseline.yaml', 'Name of configuration file', str),
                           topdir=(None, 'Top level directory', str))

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_SEDAnalysis.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        input_config = {}
        job_configs = {}
        output_config = {}

        topdir = args['topdir']
        targets_yaml = os.path.join(topdir, args['targetlist'])
        config_yaml = args['config']

        try:
            targets = load_yaml(targets_yaml)
        except IOError:
            targets = {}

        logdir = os.path.abspath(topdir).replace('gpfs', 'nfs')

        for target_name, target_list in targets.items():
            config_path = os.path.join(topdir, target_name, config_yaml)
            logfile = os.path.join(logdir, target_name, "%s_%s.log"%(self.link.linkname, target_name))

            job_config = dict(config=config_path,
                              profiles=target_list,
                              logfile=logfile)
            job_configs[target_name] = job_config

        return input_config, job_configs, output_config


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
    link.linkname = kwargs.pop('linkname', link.linkname)
    appname = kwargs.pop('appname', 'dmpipe-analyze-roi-sg')

    lsf_args = {'W': 1500,
                'R': 'rhel60'}

    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"

    config_maker = ConfigMaker_TargetAnalysis(link)
    lsf_sg = build_sg_from_link(link, config_maker,
                                lsf_args=lsf_args,
                                usage=usage,
                                description=description,
                                appname=appname,
                                **kwargs)
    return lsf_sg


def create_sg_sed_analysis(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    sed_analysis = SEDAnalysis(**kwargs)
    link = sed_analysis
    link.linkname = kwargs.pop('linkname', link.linkname)
    appname = kwargs.pop('appname', 'dmpipe-analyze-sed-sg')

    lsf_args = {'W': 1500,
                'R': 'rhel60'}

    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"

    config_maker = ConfigMaker_SEDAnalysis(link)
    lsf_sg = build_sg_from_link(link, config_maker,
                                lsf_args=lsf_args,
                                usage=usage,
                                description=description,
                                appname=appname,
                                **kwargs)
    return lsf_sg


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

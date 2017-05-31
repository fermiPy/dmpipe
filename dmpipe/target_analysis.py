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


from fermipy.utils import load_yaml

from fermipy.jobs.file_archive import FileFlags
from fermipy.jobs.chain import add_argument, Link
from fermipy.jobs.scatter_gather import ConfigMaker
from fermipy.jobs.lsf_impl import build_sg_from_link

from fermipy.gtanalysis import GTAnalysis
from fermipy.catalog import Catalog3FGL


class TargetAnalysis(object):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    NULL_MODEL = 'srcmdls/null.xml'

    default_options = dict(config=('config_baseline.yaml', 'Name of config script', str),
                           dry_run=(False, 'Print but do not run commands', bool))

    def __init__(self, **kwargs):
        """C'tor
        """
        self.parser = TargetAnalysis._make_parser()
        self.link = TargetAnalysis._make_link(self.parser, **kwargs)

    @staticmethod
    def _make_parser():
        """Make an argument parser for this class """
        usage = "dmpipe-target-analysis [options]"
        description = "Run analysis on a target"

        parser = argparse.ArgumentParser(usage=usage, description=description)
        for key, val in TargetAnalysis.default_options.items():
            add_argument(parser, key, val)
        return parser

    @staticmethod
    def _make_link(parser, **kwargs):
        link = Link(kwargs.pop('linkname', 'roi-analysis'),
                    appname='dmpipe-analyze-roi',
                    options=TargetAnalysis.default_options.copy(),
                    file_args=dict())
        return link

    def run(self, argv):
        """Run this analysis"""
        args = self.parser.parse_args(argv)

        gta = GTAnalysis(args.config,logging={'verbosity' : 3},
                         fileio={'workdir_regex' : '\.xml$|\.npy$'})
        
        gta.setup(overwrite=False)
        gta.free_sources(False)
        gta.print_roi()
        gta.optimize()
        gta.print_roi()
        
        exclude = ['3FGL J1707.8+5626']
    
        # Localize all point sources
        for s in sorted(gta.roi.sources, key=lambda t: t['ts'],reverse=True):
            #    for s in gta.roi.sources:

            if not s['SpatialModel'] == 'PointSource':
                continue
            if s['offset_roi_edge'] > -0.1:
                continue

            if s.name in exclude:
                continue
            if not '3FGL' in s.name:
                continue
        
            gta.localize(s.name,nstep=5,dtheta_max=0.5,update=True,
                         prefix='base', make_plots=True)

        gta.optimize()
        gta.print_roi()
        
        gta.write_roi('base_roi',make_plots=True)
        
        gta.find_sources(sqrt_ts_threshold=5.0)
        gta.optimize()
        gta.print_roi()
        gta.print_params()
        
        gta.free_sources(skydir=gta.roi.skydir,distance=1.0, pars='norm')
        gta.fit()
        gta.print_roi()
        gta.print_params()
        
        gta.write_roi('fit_baseline',make_plots=True)



class SEDAnalysis(object):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    NULL_MODEL = 'srcmdls/null.xml'

    default_options = dict(config=('config_baseline.yaml', 'Name of config script', str),
                           dry_run=(False, 'Print but do not run commands', bool),
                           profiles=([], 'Profiles to build SED for', list))

    def __init__(self, **kwargs):
        """C'tor
        """
        self.parser = SEDAnalysis._make_parser()
        self.link = SEDAnalysis._make_link(self.parser, **kwargs)

    @staticmethod
    def _make_parser():
        """Make an argument parser for this class """
        usage = "dmpipe-target-analysis [options]"
        description = "Run analysis on a target"

        parser = argparse.ArgumentParser(usage=usage, description=description)
        for key, val in SEDAnalysis.default_options.items():
            add_argument(parser, key, val)
        return parser

    @staticmethod
    def _make_link(parser, **kwargs):
        link = Link(kwargs.pop('linkname', 'sed-analysis'),
                    appname='dmpipe-analyze-sed',
                    options=SEDAnalysis.default_options.copy(),
                    file_args=dict())
        return link

    @staticmethod
    def _build_profile_dict(basedir, profile_name):
        """
        """
        profile_path = os.path.join(basedir, "profile_%s.yaml"%profile_name)
        profile_config = load_yaml(profile_path)

        profile_dict = {}
        profile_dict['SpatialModel'] = 'PointSource'    
        profile_dict['SpectrumType'] = 'PowerLaw'
        return profile_name, profile_dict


    def run(self, argv):
        """Run this analysis"""
        args = self.parser.parse_args(argv)
        gta = GTAnalysis(args.config,
                         logging={'verbosity' : 3},
                         fileio={'workdir_regex' : '\.xml$|\.npy$'})
        gta.setup(overwrite=False)
        gta.load_roi('fit_baseline')
        gta.print_roi()

        basedir = os.path.dirname(args.config)
        # This should be a no-op, b/c it was done in the baseline analysis
    
        gta.free_sources(skydir=gta.roi.skydir,distance=1.0, pars='norm')
        
        for profile in args.profiles:
            pkey, pdict = SEDAnalysis._build_profile_dict(basedir, profile)
            # test_case need to be a dict with spectrum and morphology
            gta.add_source(pkey, pdict)
            # refit the ROI
            gta.fit()
            # build the SED
            gta.sed(pkey, outfile="sed_%s.fits"%pkey)
            # remove the source
            gta.delete_source(pkey)
            # put the ROI back to how it was
            gta.load_xml('fit_draco_baseline')

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
        
        targets = load_yaml(targets_yaml)
        for target_name in targets.keys():
            config_path = os.path.join(topdir, target_name, config_yaml)
            job_config = dict(config=config_path)
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
        
        targets = load_yaml(targets_yaml)
        for target_name, target_list in targets.items():
            config_path = os.path.join(topdir, target_name, config_yaml)
            job_config = dict(config=config_path,
                              profiles=target_list)
            job_configs[target_name] = job_config

        return input_config, job_configs, output_config



def create_sg_roi_analysis(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    roi_analysis = TargetAnalysis()
    link = roi_analysis.link
    link.linkname = kwargs.pop('linkname', link.linkname)
    appname = kwargs.pop('appname', 'dmpipe-analyze-roi-sg')

    lsf_args = {'W': 1500,
                'R': 'rhel60'}

    usage = "%s [options]"%(appname)
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
    sed_analysis = SEDAnalysis()
    link = sed_analysis.link
    link.linkname = kwargs.pop('linkname', link.linkname)
    appname = kwargs.pop('appname', 'dmpipe-analyze-sed-sg')

    lsf_args = {'W': 1500,
                'R': 'rhel60'}

    usage = "%s [options]"%(appname)
    description = "Run analyses on a series of ROIs"

    config_maker = ConfigMaker_SEDAnalysis(link)
    lsf_sg = build_sg_from_link(link, config_maker,
                                lsf_args=lsf_args,
                                usage=usage,
                                description=description,
                                appname=appname,
                                **kwargs)
    return lsf_sg


def main_roi_single():
    """ Entry point for analysis of a single ROI """
    target_analysis = TargetAnalysis()
    target_analysis.run(sys.argv[1:])

def main_sed_single():
    """ Entry point for analysis of a single ROI """
    sed_analysis = SEDAnalysis()
    sed_analysis.run(sys.argv[1:])

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

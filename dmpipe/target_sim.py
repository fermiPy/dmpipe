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

from fermipy.utils import load_yaml, write_yaml, init_matplotlib_backend

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro

from fermipy import utils
from fermipy.jobs.chain import Link
from fermipy.jobs.scatter_gather import ConfigMaker, build_sg_from_link
from fermipy.jobs.lsf_impl import make_nfs_path, get_lsf_default_args, LSF_Interface

from dmpipe.name_policy import NameFactory

init_matplotlib_backend('Agg')

try:
    from fermipy.gtanalysis import GTAnalysis
    HAVE_ST = True
except ImportError:
    HAVE_ST = False

NAME_FACTORY = NameFactory(basedir=('.'))


class TargetSim(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    default_options = dict(config=('config.yaml', 'Name of config script', str),
                           sim=(None, 'Name of simulation configuration file', str),
                           nsims=(20, 'Number of realizations to simulate', int),
                           seed=(0, 'Seed to use for first realization', int),
                           dry_run=(False, 'Print but do not run commands', bool))

    def __init__(self, **kwargs):
        """C'tor
        """
        parser = argparse.ArgumentParser(usage="dmpipe-simulate-roi [options]",
                                         description="Run simulated analysis of a single ROI")
        Link.__init__(self, kwargs.pop('linkname', 'analyze-roi'),
                      parser=parser,
                      appname='dmpipe-simulate-roi',
                      options=TargetSim.default_options.copy(),
                      **kwargs)

    def run_simulation(self, gta, injected_source, test_source, seed, outfile):
        """Simulate a realization of this analysis"""
        gta.load_roi('fit_baseline')
        gta.set_random_seed(seed)
        if injected_source:            
            gta.add_source(injected_source['name'], injected_source['source_model'])
            print (gta.like)
        gta.simulate_roi()
        if injected_source:
            gta.delete_source(injected_source['name'])
        gta.optimize()
        gta.find_sources(sqrt_ts_threshold=5.0, search_skydir=gta.roi.skydir,
                         search_minmax_radius=[1.0, np.nan])
        gta.optimize()
        gta.free_sources(skydir=gta.roi.skydir, distance=1.0, pars='norm')
        gta.fit()
        gta.add_source(test_source['name'], test_source['source_model'])
        gta.sed(test_source['name'], outfile=outfile)     
        # Set things back to how they were
        gta.delete_source(test_source['name'])    
        return outfile
        

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if not HAVE_ST:
            raise RuntimeError("Trying to run fermipy analysis, but don't have ST")
        
        gta = GTAnalysis(args.config, logging={'verbosity': 3},
                         fileio={'workdir_regex': '\.xml$|\.npy$'})

        profilefile = args.config.replace('config.yaml','profile_default.yaml')
        profile = utils.load_yaml(profilefile)

        sim_config = utils.load_yaml(args.sim)
        
        injected_source = sim_config.get('injected_source', None)
        if injected_source is not None:
            injected_source['source_model']['norm'] = dict(value=profile['j_integ'])
            injected_source['source_model']['ra'] = gta.config['selection']['ra']
            injected_source['source_model']['dec'] = gta.config['selection']['dec']
        test_source = sim_config['test_source']
        
        first = args.seed
        last = first + args.nsims
        for i in range(first, last):
            sedfile = "sed_%s_%06i.fits"%(test_source['name'], i)
            self.run_simulation(gta, injected_source, test_source, i, sedfile)


class ConfigMaker_TargetSim(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(targetlist=('target_list.yaml', 'Yaml file with list of targets', str),
                           sim=(None, 'Name of simulation configuration file', str),
                           nsims=(20, 'Number of realizations to simulate', int),
                           seed=(0, 'Seed to use for first realization', int),
                           topdir=(None, 'Top level directory', str))

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_TargetSim.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        topdir = args['topdir']
        sim = args['sim']
        targets_yaml = os.path.join("%s_sim"%topdir, "sim_%s"%sim, args['targetlist'])
        config_yaml = 'config.yaml'

        try:
            targets = load_yaml(targets_yaml)
        except IOError:
            targets = {}

        for target_name in targets.keys():
            name_keys = dict(target_type=topdir,
                             target_name=target_name,
                             sim_name=sim,
                             fullpath=True)
            simdir = NAME_FACTORY.sim_targetdir(**name_keys)
            simyamlfile = os.path.join(simdir, 'sim_%s.yaml'%sim)
            config_path = os.path.join(simdir, 'config.yaml')
            logfile = make_nfs_path(os.path.join(simdir, "%s_%s.log"%(self.link.linkname, target_name)))
            job_config = dict(config=config_path, 
                              logfile=logfile,
                              sim=simyamlfile,
                              nsims=args['nsims'],
                              seed=args['seed'])
            job_configs[target_name] = job_config

        return job_configs




def create_link_roi_sim(**kwargs):
    """Build and return a `Link` object that can invoke TargetAnalysis"""
    target_analysis = TargetSim(**kwargs)
    return target_analysis



def create_sg_roi_sim(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    roi_analysis = TargetSim(**kwargs)
    link = roi_analysis

    appname = kwargs.pop('appname', 'dmpipe-simulated-roi-sg')

    batch_args = get_lsf_default_args()    
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"

    config_maker = ConfigMaker_TargetSim(link)
    sg = build_sg_from_link(link, config_maker,
                            interface=batch_interface,
                            usage=usage,
                            description=description,
                            appname=appname,
                            **kwargs)
    return sg


def main_roi_single():
    """ Entry point for analysis of a single ROI """
    target_analysis = TargetSim()
    target_analysis.run_analysis(sys.argv[1:])


def main_roi_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_roi_sim()
    lsf_sg(sys.argv)


if __name__ == "__main__":
    main_roi_single()

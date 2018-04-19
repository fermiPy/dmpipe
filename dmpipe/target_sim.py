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

from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic
from gammapy.maps import WcsGeom

from fermipy.utils import load_yaml, write_yaml, init_matplotlib_backend

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro

from fermipy import utils
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

class RandomSkyDirGenerator(Link):
    """Small class to generate random sky directions inside an ROI

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    default_options = dict(config=defaults.common['config'],
                           rand_config=defaults.sims['rand_config'],
                           outfile=defaults.generic['outfile'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs): 
        """C'tor
        """
        parser = argparse.ArgumentParser(usage="dmpipe-random-dir-gen [options]",
                                         description="Generate random sky directions in an ROI")
        Link.__init__(self, kwargs.pop('linkname', 'random-skydirs'),
                      parser=parser,
                      appname='dmpipe-random-dir-gen',
                      options=RandomSkyDirGenerator.default_options.copy(),
                      **kwargs)

    @staticmethod
    def make_wcsgeom_from_config(config):    
        binning = config['binning']
        binsz = binning['binsz']
        coordsys = binning.get('coordsys', 'GAL')
        roiwidth = binning['roiwidth']
        proj = binning.get('proj', 'AIT')
        ra = config['selection']['ra']
        dec = config['selection']['dec']
        npix = int(np.round(roiwidth / binsz))
        skydir = SkyCoord(ra * u.deg, dec * u.deg)
        
        wcsgeom = WcsGeom.create(npix=npix, binsz=binsz,
                                 proj=proj, coordsys=coordsys,
                                 skydir=skydir)
        return wcsgeom
      
    @staticmethod
    def build_skydir_dict(wcsgeom, rand_config):    
        """Build a dictionary of random directions"""
        step_x = rand_config['step_x']
        step_y = rand_config['step_y']
        max_x = rand_config['max_x']
        max_y = rand_config['max_y']
        seed = rand_config['seed']
        nsims = rand_config['nsims']
        
        cdelt = wcsgeom.wcs.wcs.cdelt
        pixstep_x = step_x / cdelt[0]
        pixstep_y = -1. * step_y / cdelt[1]
        pixmax_x = max_x / cdelt[0]
        pixmax_y = max_y / cdelt[0]
        
        nstep_x = int(np.ceil( 2.*pixmax_x / pixstep_x)) + 1
        nstep_y = int(np.ceil( 2.*pixmax_y / pixstep_y)) + 1

        print (nstep_x, nstep_y)

        center = np.array(wcsgeom._center_pix)

        grid = np.meshgrid(np.linspace(-1*pixmax_x, pixmax_x, nstep_x),
                           np.linspace(-1*pixmax_y, pixmax_y, nstep_y))
        grid[0] += center[0]
        grid[1] += center[1]

        test_grid = wcsgeom.pix_to_coord(grid)     
        glat_vals = test_grid[0].flat
        glon_vals = test_grid[1].flat
        conv_vals = SkyCoord(glat_vals * u.deg, glon_vals * u.deg, frame=Galactic).transform_to(ICRS)

        ra_vals = conv_vals.ra.deg[seed:nsims]
        dec_vals = conv_vals.dec.deg[seed:nsims]

        o_dict = {}
        for i, (ra, dec) in enumerate(zip(ra_vals, dec_vals)):
            key = i+seed
            o_dict[key] = dict(ra=ra, dec=dec)
        return o_dict


    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if args.config in [None, 'none', 'None']:
            raise ValueError("Config yaml file must be specified")
        if args.rand_config in [None, 'none', 'None']:
            raise ValueError("Random direction config yaml file must be specified")
        config = load_yaml(args.config)
        rand_config = load_yaml(args.rand_config)

        wcsgeom = RandomSkyDirGenerator.make_wcsgeom_from_config(config)
        dir_dict = RandomSkyDirGenerator.build_skydir_dict(wcsgeom, rand_config)

        if args.outfile not in [None, 'none', 'None']:
            write_yaml(dir_dict, args.outfile)


class TargetSim(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    default_options = dict(config=defaults.common['config'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'])

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

    def run_simulation(self, gta, injected_source, test_source, seed, outfile, mcube_file=None):
        """Simulate a realization of this analysis"""
        gta.load_roi('fit_baseline')
        gta.set_random_seed(seed)
        if injected_source:            
            gta.add_source(injected_source['name'], injected_source['source_model'])
            if mcube_file is not None:
                gta.write_model_map(mcube_file)
                mc_spec_dict = dict(true_counts=gta.model_counts_spectrum(injected_source['name']),
                                    energies=gta.energies,
                                    model=injected_source['source_model'])
                mcspec_file = os.path.join(gta.workdir, "mcspec_%s.yaml"%mcube_file)
                utils.write_yaml(mc_spec_dict, mcspec_file)

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

        workdir = os.path.dirname(args.config)
        profilefile = os.path.join(workdir, 'profile_default.yaml')
        simfile = os.path.join(workdir, 'sim_%s.yaml'%args.sim)
        profile = utils.load_yaml(profilefile)
        sim_config = utils.load_yaml(simfile)
        
        injected_source = sim_config.get('injected_source', None)
        if injected_source is not None:
            injected_source['source_model']['norm'] = dict(value=profile['j_integ'])
            injected_source['source_model']['ra'] = gta.config['selection']['ra']
            injected_source['source_model']['dec'] = gta.config['selection']['dec']
            mcube_file = args.sim
        else:
            mcube_file = None

        test_source = sim_config['test_source']
        
        first = args.seed
        last = first + args.nsims

        for i in range(first, last):
            sedfile = "sed_%s_%06i.fits"%(test_source['name'], i)
            if i == first:
                mcube_out = mcube_file
            else:
                mcube_out = None
            self.run_simulation(gta, injected_source, test_source, i, sedfile, mcube_out)


class ConfigMaker_RandomDirGen(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           config=defaults.common['config'],
                           rand_config=defaults.sims['rand_config'], 
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        ConfigMaker.__init__(self, link,
                             options=kwargs.get('options',
                                                ConfigMaker_RandomDirGen.default_options.copy()))

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
                             sim_name='random',
                             fullpath=True)
            simdir = NAME_FACTORY.sim_targetdir(**name_keys)
            config_path = os.path.join(simdir, config_yaml)
            outfile = os.path.join(simdir, 'skydirs.yaml')
            logfile = make_nfs_path(outfile.replace('yaml', 'log'))
            job_config = dict(config=config_path, 
                              rand_config=args['rand_config'],
                              outfile=outfile,
                              logfile=logfile,
                              dry_run=args['dry_run'])
            job_configs[target_name] = job_config

        return job_configs



class ConfigMaker_TargetSim(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           config=defaults.common['config'],
                           sim=defaults.sims['sim'], 
                           nsims=defaults.sims['nsims'], 
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'])

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
                             sim_name=sim,
                             fullpath=True)
            simdir = NAME_FACTORY.sim_targetdir(**name_keys)
            config_path = os.path.join(simdir, config_yaml)
            logfile = make_nfs_path(os.path.join(simdir, "%s_%s.log"%(self.link.linkname, target_name)))
            job_config = dict(config=config_path, 
                              logfile=logfile,
                              sim=sim,
                              nsims=args['nsims'],
                              seed=args['seed'])
            job_configs[target_name] = job_config

        return job_configs


def create_link_random_dir_gen(**kwargs):
    """Build and return a `Link` object that can invoke RandomSkyDirGenerator"""
    random_dir_gen = RandomSkyDirGenerator(**kwargs)
    return random_dir_gen

def create_link_roi_sim(**kwargs):
    """Build and return a `Link` object that can invoke TargetAnalysis"""
    target_analysis = TargetSim(**kwargs)
    return target_analysis


def create_sg_random_dir_gen(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    random_dir_gen = RandomSkyDirGenerator(**kwargs)
    link = random_dir_gen

    appname = kwargs.pop('appname', 'dmpipe-random-dir-gen-sg')

    batch_args = get_lsf_default_args()    
    batch_interface = LSF_Interface(**batch_args)

    usage = "%s [options]" % (appname)
    description = "Create random directions for a set of ROIs"

    config_maker = ConfigMaker_RandomDirGen(link)
    sg = build_sg_from_link(link, config_maker,
                            interface=batch_interface,
                            usage=usage,
                            description=description,
                            appname=appname,
                            **kwargs)
    return sg



def create_sg_roi_sim(**kwargs):
    """Build and return a ScatterGather object that can invoke this script"""
    roi_analysis = TargetSim(**kwargs)
    link = roi_analysis

    appname = kwargs.pop('appname', 'dmpipe-simulate-roi-sg')

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


def main_random_dir_gen():
    """ Entry point for analysis of a single ROI """
    random_dir_gen = RandomSkyDirGenerator()
    random_dir_gen.run_analysis(sys.argv[1:])

def main_roi_single():
    """ Entry point for analysis of a single ROI """
    target_analysis = TargetSim()
    target_analysis.run_analysis(sys.argv[1:])

def main_random_dir_gen_batch():
    """ Entry point for analysis of a single ROI """
    lsf_sg = create_sg_random_dir_gen()
    lsf_sg(sys.argv)
 
def main_roi_batch():
    """ Entry point for command line use for dispatching batch jobs """
    lsf_sg = create_sg_roi_sim()
    lsf_sg(sys.argv)


if __name__ == "__main__":
    main_roi_single()

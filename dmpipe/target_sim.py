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

from astropy import units as u
from astropy.coordinates import SkyCoord, ICRS, Galactic
from gammapy.maps import WcsGeom

from fermipy.utils import load_yaml, write_yaml, init_matplotlib_backend

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro

from fermipy import utils
from fermipy.jobs.utils import is_null, is_not_null
from fermipy.jobs.link import Link
from fermipy.jobs.scatter_gather import ConfigMaker, build_sg_from_link
from fermipy.jobs.slac_impl import make_nfs_path, get_slac_default_args, Slac_Interface

from dmpipe.name_policy import NameFactory
from dmpipe import defaults

init_matplotlib_backend('Agg')

try:
    from fermipy.gtanalysis import GTAnalysis
    HAVE_ST = True
except ImportError:
    HAVE_ST = False

NAME_FACTORY = NameFactory(basedir=('.'))


class CopyBaseROI(Link):
    """Small class to copy a baseline ROI to a simulation area

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    appname = 'dmpipe-copy-base-roi'
    linkname_default = 'copy-base-roi'
    usage = '%s [options]' %(appname)
    description = "Copy a baseline ROI to a simulation area"
    
    default_options = dict(ttype=defaults.common['ttype'],                           
                           target=defaults.common['target'],
                           profiles=defaults.common['profiles'],
                           sim=defaults.sims['sim'],
                           dry_run=defaults.common['dry_run'])

    copyfiles = ['srcmap_00.fits', 'fit_baseline.fits', 'fit_baseline.npy', 'fit_baseline_00.xml']

    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(CopyBaseRoi, self).__init__(linkname, **init_dict)

    @classmethod
    def copy_analysis_files(cls, orig_dir, dest_dir, files):
        """ Copy a list of files from orig_dir to dest_dir"""
        for f in files:
            orig_path = os.path.join(orig_dir, f)
            dest_path = os.path.join(dest_dir, f)
            try:
                copyfile(orig_path, dest_path)
            except IOError:
                sys.stderr.write("WARNING: failed to copy %s\n"%orig_path)

    @classmethod
    def copy_target_dir(cls, orig_dir, dest_dir, sim, profiles):
        """ Create and populate directoris for target analysis
        """
        target_dict = {}
        
        try:
            os.makedirs(dest_dir)
        except OSError:
            pass

        orig_config_path = os.path.join(orig_dir, 'config.yaml')
        dest_config_path = os.path.join(dest_dir, 'config.yaml')
        self.copy_analysis_files(orig_dir, dest_dir, self.copyfiles)

        dest_config = load_yaml(orig_config_path).copy()
        dest_config['gtlike']['bexpmap'] = os.path.abspath(os.path.join(orig_dir,'bexpmap_00.fits'))
        dest_config['gtlike']['srcmap'] = os.path.abspath(os.path.join(dest_dir,'srcmap_00.fits'))
        dest_config['gtlike']['use_external_srcmap'] = True
        write_yaml(target_config, dest_config_path)

        sim_orig_file = os.path.join('config', 'sim_%s.yaml'%sim)
        sim_dest_file = os.path.join(dest_dir, 'sim_%s.yaml'%sim)
        copyfile(sim_orig_file, sim_dest_file)
               
        for profile in profiles:
            orig_profile_path = os.path.join(orig_dir, 'profile_%s.yaml'%profile)
            dest_profile_path = os.path.join(dest_dir, 'profile_%s.yaml'%profile)
            jmap_path = dest_profile_path.replace('.yaml', 'fits')
            dest_profile = oad_yaml(orig_profile_path).copy()
            dest_profile['j_map_file'] = jmap_path
            write_yaml(dest_profile, dest_profile_path)


    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        name_keys = dict(target_type=args.ttype,
                         target_name=args.target,
                         profiles=args.profiles,
                         sim_name=args.sim,
                         fullpath=True)

        orig_dir = NAME_FACTORY.targetdir(**name_keys)
        dest_dir = NAME_FACTORY.sim_targetdir(**name_keys)
        self.copy_target_dir(orig_dir, dest_dir, args.sim, args.profiles)


class CopyBaseROI_SG(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-copy-base-roi-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = CopyBaseROI

    batch_args = get_slac_default_args()    
    batch_interface = Slac_Interface(**batch_args)

    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           rosterlist=defaults.common['rosterlist'],
                           sim=defaults.sim['sim'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        super(CopyBaseROI_SG, self).__init__(link,
                                            options=kwargs.get('options',
                                                               self.default_options.copy()))

    def build_job_configs(self, args):
        """Hook to build job configurations
        """
        job_configs = {}

        orig_args = args.copy()
        orig_args.pop('sim')
        (data_targets_yaml, dummy) = NAME_FACTORY.resolve_targetfile(orig_args)
        (data_roster_yaml, dummy) = NAME_FACTORY.resolve_rosterfile(orig_args)

        ttype = args['ttype']
        dry_run = args['dry_run']
        (sim_targets_yaml, sim) = NAME_FACTORY.resolve_targetfile(args)
        (sim_roster_yaml, sim) = NAME_FACTORY.resolve_rosterfile(args)
             
        copyfile(data_targets_yaml, sim_targets_yaml)
        copyfile(data_roster_yaml, sim_roster_yaml)

        targets = load_yaml(targets_yaml)

        for target_name, target_list in targets.items():
            job_config = dict(ttype=ttype,
                              target=target_name,
                              profiles=target_list,
                              sim=sim, 
                              dry_run=dry_run)
            job_configs[target_name] = job_config

        return job_configs

class RandomDirGen(Link):
    """Small class to generate random sky directions inside an ROI

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    appname = 'dmpipe-random-dir-gen'
    linkname_default = 'random-dir-gen'
    usage = '%s [options]' %(appname)
    description = "Generate random sky directions in an ROI"

    default_options = dict(config=defaults.common['config'],
                           rand_config=defaults.sims['rand_config'],
                           outfile=defaults.generic['outfile'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs): 
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(RandomDirGen, self).__init__(linkname, **init_dict)

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

        if is_null(args.config):
            raise ValueError("Config yaml file must be specified")
        if is_null(args.rand_config):
            raise ValueError("Random direction config yaml file must be specified")
        config = load_yaml(args.config)
        rand_config = load_yaml(args.rand_config)

        wcsgeom = self.make_wcsgeom_from_config(config)
        dir_dict = self.build_skydir_dict(wcsgeom, rand_config)

        if is_not_null(args.outfile):
            write_yaml(dir_dict, args.outfile)


class SimulateROI(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    appname = 'dmpipe-simulate-roi'
    linkname_default = 'simulate-roi'
    usage = '%s [options]' %(appname)
    description = "Run simulated analysis of a single ROI"

    default_options = dict(config=defaults.common['config'],
                           sim=defaults.sims['sim'],
                           nsims=defaults.sims['nsims'],
                           seed=defaults.sims['seed'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(SimulateROI, self).__init__(linkname, **init_dict)

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


class RandomDirGen_SG(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-random-dir-gen-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = RandomDirGen

    batch_args = get_slac_default_args()    
    batch_interface = Slac_Interface(**batch_args)

    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           config=defaults.common['config'],
                           rand_config=defaults.sims['rand_config'], 
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        super(RandomDirGen_SG, self).__init__(link,
                                              options=kwargs.get('options',
                                                                 self.default_options.copy()))

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
        if is_not_null(config_override):
            config_yaml = config_override

        rand_yaml = NAME_FACTORY.resolve_randconfig(args)        

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
                              rand_config=rand_yaml,
                              outfile=outfile,
                              logfile=logfile,
                              dry_run=args['dry_run'])
            job_configs[target_name] = job_config

        return job_configs



class SimulateROI_SG(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-simulate-roi-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = SimulateROI

    batch_args = get_slac_default_args()    
    batch_interface = Slac_Interface(**batch_args)

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
        super(SimulateROI_SG, self).__init__(link,
                                             options=kwargs.get('options',
                                                                self.default_options.copy()))

 
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
        if is_not_null(config_override):
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

def register_classes():
    CopyBaseROI.register_class()
    CopyBaseROI_SG.register_class()
    SimulateROI.register_class()
    SimulateROI_SG.register_class()
    RandomDirGen.register_class()
    RandomDirGen_SG.register_class()

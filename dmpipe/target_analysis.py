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

class PrepareTargets(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    appname = 'dmpipe-prepare-targets'
    linkname_default = 'prepare-targets'
    usage = '%s [options]' %(appname)
    description = "Prepare directories for target analyses"

    default_options = dict(ttype=defaults.common['ttype'],
                           rosters=defaults.common['rosters'],
                           config=defaults.common['config'],
                           sims=defaults.sims['sims'],
                           dry_run=defaults.common['dry_run'])


    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(PrepareTargets, self).__init__(linkname, **init_dict)


    @classmethod
    def write_target_dirs(cls, ttype, roster_dict, base_config, sims):
        """ Create and populate directoris for target analysis
        """
        target_dict = {}

        target_info_dict = {}
        roster_info_dict = {}

        try:
            os.makedirs(ttype)
        except OSError:
            pass

        for roster_name, rost in roster_dict.items():
            tlist = []
            for target_name, target in rost.items():
                
                if target.ver_key is None:
                    target_verkey = target.version
                else:
                    target_verkey = target.ver_key
                target_key = "%s:%s" % (target_name, target_verkey)
                print("Writing %s" % (target_key))
                tlist.append(target_key)
                name_keys = dict(target_type=ttype,
                                 target_name=target_name,
                                 profile=target_verkey,
                                 fullpath=True)
                
                target_dir = NAME_FACTORY.targetdir(**name_keys)
                profile_path = NAME_FACTORY.profilefile(**name_keys)
                j_val_path = profile_path.replace('profile_','j_val_')
                target_config_path = os.path.join(target_dir, 'config.yaml')
  
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

                if target_info_dict.has_key(target_name):
                    target_info_dict[target_name].append(target_verkey)
                else:
                    target_info_dict[target_name] = [target_verkey]
         
                source_model = dict(SpectrumType='PowerLaw',
                                    RA=target.ra,
                                    DEC=target.dec)
                                                    
                if target.proftype in [None, 'point', 'Point']:                                        
                    source_model.update(dict(SpatialModel='PointSource'))
                else:
                    target.j_rad_file = profile_path.replace('.yaml', '.dat')
                    target.write_j_rad_file()
                    source_model.update(dict(SpatialModel='DiffuseSource',
                                             SpatialType='RadialProfile',
                                             radialprofile=target.j_rad_file))
                    
                profile_dict = dict(name = target_verkey,
                                    source_model=source_model)
                write_yaml(profile_dict, profile_path)

                j_profile_data = target.profile.copy()
                j_profile_data['j_integ'] = target.j_integ
                j_profile_data['j_sigma'] = target.j_sigma

                write_yaml(j_profile_data, j_val_path)

            roster_info_dict[roster_name] = tlist

        roster_file = os.path.join(ttype, 'roster_list.yaml')
        target_file = os.path.join(ttype, 'target_list.yaml')

        write_yaml(roster_info_dict, roster_file)
        write_yaml(target_info_dict, target_file)

        print ('sims',sims)
        for sim in sims:
            sim_dir = os.path.join("%s_sim"%ttype, "sim_%s"%sim)
            sim_roster_file = os.path.join(sim_dir, 'roster_list.yaml')
            sim_target_file = os.path.join(sim_dir, 'target_list.yaml')
            try:
                os.makedirs(sim_dir)
            except OSError:
                pass
            copyfile(roster_file, sim_roster_file)
            copyfile(target_file, sim_target_file)

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)
        roster_lib = RosterLibrary()
        roster_dict = {}

        if len(args.rosters) == 0:
            sys.stderr.write("You must specify at least one target roster")
            return -1
        
        if is_null(args.ttype):
            sys.stderr.write("You must specify a target type")
            return -1

        if is_null(args.sims):
            sims = []
        else:
            sims = args.sims

        name_keys = dict(target_type=args.ttype,
                         fullpath=True)
        config_file = NAME_FACTORY.ttypeconfig(**name_keys)
        if is_not_null(args.config):
            config_file = args.config

        for roster in args.rosters:
            rost = roster_lib.create_roster(roster)
            roster_dict[roster] = rost

        base_config = load_yaml(config_file)
        self.write_target_dirs(args.ttype, roster_dict, base_config, sims)
        

class AnalyzeROI(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    appname = 'dmpipe-analyze-roi'
    linkname_default = 'analyze-roi'
    usage = '%s [options]' %(appname)
    description = "Run analysis of a single ROI"

    default_options = dict(config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(AnalyzeROI, self).__init__(linkname, **init_dict)

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


class AnalyzeSED(Link):
    """Small class wrap an analysis script.

    This is useful for parallelizing analysis using the fermipy.jobs module.
    """
    appname = 'dmpipe-analyze-sed'
    linkname_default = 'analyze-sed'
    usage = '%s [options]' %(appname)
    description = "Extract the SED for a single target"

    default_options = dict(config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'],
                           skydirs=defaults.sims['skydirs'],
                           profiles=([], 'Profiles to build SED for', list))
                           
    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(AnalyzeSED, self).__init__(linkname, **init_dict)
                          
    @staticmethod
    def _build_profile_dict(basedir, profile_name):
        """
        """
        profile_path = os.path.join(basedir, "profile_%s.yaml" % profile_name)
        profile_config = load_yaml(profile_path)
        if profile_name != profile_config['name']:
            sys.stderr.write('Warning, profile name (%s) != name in %s (%s)\n'%(profile_name, profile_config['name'], profile_path))
            
        profile_dict = profile_config['source_model']
        return profile_name, profile_dict

    def run_analysis(self, argv):
        """Run this analysis"""
        args = self._parser.parse_args(argv)

        if not HAVE_ST:
            raise RuntimeError("Trying to run fermipy analysis, but don't have ST")

        if is_null(args.skydirs):
            skydir_dict = None
        else:
            skydir_dict = load_yaml(args.skydirs)

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
            if skydir_dict is None:
                skydir_keys = [None]
            else:
                skydir_keys = sorted(skydir_dict.keys())

            for skydir_key in skydir_keys:                
                if skydir_key is None:
                    pkey, pdict = AnalyzeSED._build_profile_dict(basedir, profile)
                else:                    
                    skydir_val = skydir_dict[skydir_key]
                    pkey, pdict = AnalyzeSED._build_profile_dict(basedir, profile)
                    pdict['ra'] = skydir_val['ra']
                    pdict['dec'] = skydir_val['dec']
                    pkey += "_%06i"%skydir_key

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


class AnalyzeROI_SG(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-analyze-roi-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = AnalyzeROI

    batch_args = get_slac_default_args()    
    batch_interface = Slac_Interface(**batch_args)

    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        super(AnalyzeROI_SG, self).__init__(link,
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
                             fullpath=True)
            target_dir = NAME_FACTORY.targetdir(**name_keys)
            config_path = os.path.join(target_dir, config_yaml)
            logfile = make_nfs_path(os.path.join(target_dir, "%s_%s.log"%(self.link.linkname, target_name)))
            job_config = dict(config=config_path, 
                              logfile=logfile)
            job_configs[target_name] = job_config

        return job_configs


class AnalyzeSED_SG(ConfigMaker):
    """Small class to generate configurations for this script

    This adds the following arguments:
    """
    appname = 'dmpipe-analyze-sed-sg'
    usage = "%s [options]" % (appname)
    description = "Run analyses on a series of ROIs"
    clientclass = AnalyzeSED

    batch_args = get_slac_default_args()    
    batch_interface = Slac_Interface(**batch_args)

    default_options = dict(ttype=defaults.common['ttype'],
                           targetlist=defaults.common['targetlist'],
                           config=defaults.common['config'],
                           skydirs=defaults.sims['skydirs'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, link, **kwargs):
        """C'tor
        """
        super(AnalyzeSED_SG, self).__init__(link,
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

        targets = load_yaml(targets_yaml)
        config_yaml = 'config.yaml'

        if is_not_null(args['skydirs']):
            skydirs = args['skydirs']
        else:
            skydirs = None
 
        for target_name, target_list in targets.items():
            name_keys = dict(target_type=ttype,
                             target_name=target_name,
                             sim_name='random',
                             fullpath=True)
            if skydirs is None:
                target_dir = NAME_FACTORY.targetdir(**name_keys)
                skydir_path = None
            else:
                target_dir = NAME_FACTORY.sim_targetdir(**name_keys)
                skydir_path = os.path.join(target_dir, skydirs)                                          
            config_path = os.path.join(target_dir, config_yaml)
            logfile =  make_nfs_path(os.path.join(target_dir, "%s_%s.log"%(self.link.linkname, target_name)))
            job_config = dict(config=config_path,
                              profiles=target_list,
                              skydirs=skydir_path,
                              logfile=logfile)
            job_configs[target_name] = job_config

        return job_configs

def register_classes():
    AnalyzeROI.register_class()
    AnalyzeROI_SG.register_class()
    AnalyzeSED.register_class()
    AnalyzeSED_SG.register_class()
    PrepareTargets.register_class()


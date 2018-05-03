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
import copy

from shutil import copyfile

from dmsky.roster import RosterLibrary

from fermipy.utils import load_yaml, write_yaml, init_matplotlib_backend

from fermipy.castro import CastroData
from fermipy.sed_plotting import plotCastro

from fermipy.jobs.utils import is_null, is_not_null
from fermipy.jobs.link import Link
from fermipy.jobs.scatter_gather import ConfigMaker, build_sg_from_link
from fermipy.jobs.slac_impl import make_nfs_path, get_slac_default_args, SlacInterface

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
    usage = '%s [options]' % (appname)
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
    def _write_data_target_config(cls, base_config, target, target_dir):
        target_config_path = os.path.join(target_dir, 'config.yaml')
        target_config = base_config.copy()
        target_config['selection']['ra'] = target.ra
        target_config['selection']['dec'] = target.dec
        try:
            os.makedirs(target_dir)
        except OSError:
            pass
        write_yaml(target_config, target_config_path)
        return target_config

    @classmethod
    def _write_sim_target_config(cls, target_config, target_dir, sim_target_dir):
        sim_target_config_path = os.path.join(sim_target_dir, 'config.yaml')
        sim_target_config = copy.deepcopy(target_config)
        try:
            os.makedirs(sim_target_dir)
        except OSError:
            pass

        sim_target_config['gtlike']['bexpmap'] = os.path.abspath(
            os.path.join(target_dir, 'bexpmap_00.fits'))
        sim_target_config['gtlike']['srcmap'] = os.path.abspath(
            os.path.join(target_dir, 'srcmap_00.fits'))
        sim_target_config['gtlike']['use_external_srcmap'] = True

        write_yaml(sim_target_config, sim_target_config_path)
        return sim_target_config

    @classmethod
    def _write_profile_yaml(cls, target, profile_path, target_verkey):

        source_model = dict(SpectrumType='PowerLaw',
                            RA=target.ra,
                            DEC=target.dec)

        if target.proftype in [None, 'point', 'Point']:
            source_model.update(dict(SpatialModel='PointSource'))
        elif target.proftype in ['Map']:
            source_model.update(dict(SpatialModel='DiffuseSource',
                                     SpatialType='SpatialMap',
                                     Spatial_Filename=target.j_map_file))
        else:
            target.j_rad_file = profile_path.replace('.yaml', '.dat')
            target.write_j_rad_file()
            source_model.update(dict(SpatialModel='DiffuseSource',
                                     SpatialType='RadialProfile',
                                     radialprofile=target.j_rad_file))

        profile_dict = dict(name=target_verkey,
                            source_model=source_model)
        write_yaml(profile_dict, profile_path)
        return profile_dict

    @classmethod
    def _write_j_value_yaml(cls, target, j_val_path):

        j_profile_data = target.profile.copy()
        j_profile_data['j_integ'] = target.j_integ
        j_profile_data['j_sigma'] = target.j_sigma

        write_yaml(j_profile_data, j_val_path)
        return j_profile_data

    @classmethod
    def _write_sim_yaml(cls, target, sim, sim_target_dir, target_verkey):

        sim_profile_yaml = os.path.join('config', 'sim_%s.yaml' % sim)
        sim_profile = load_yaml(sim_profile_yaml)
        injected_source = sim_profile.get('injected_source', None)
        if injected_source is not None:
            sim_profile['injected_source']['source_model']['norm']['value'] = target.j_integ
        sim_out_path = os.path.join(sim_target_dir, 'sim_%s_%s.yaml' % (sim, target_verkey))
        write_yaml(sim_profile, sim_out_path)
        return sim_profile

    @classmethod
    def _write_target_dirs(cls, ttype, roster_dict, base_config, sims):
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

                if target_name in target_info_dict:
                    target_info_dict[target_name].append(target_verkey)
                else:
                    target_info_dict[target_name] = [target_verkey]
                tlist.append(target_key)

                name_keys = dict(target_type=ttype,
                                 target_name=target_name,
                                 profile=target_verkey,
                                 fullpath=True)
                target_dir = NAME_FACTORY.targetdir(**name_keys)
                profile_path = NAME_FACTORY.profilefile(**name_keys)
                j_val_path = NAME_FACTORY.j_valuefile(**name_keys)

                write_config = False
                if target_name in target_dict:
                    # Already made the config for this target
                    target_config = target_dict[target_name].copy()
                else:
                    # Make the config for this target
                    target_config = cls._write_data_target_config(base_config, target, target_dir)
                    target_dict[target_name] = target_config
                    write_config = True

                cls._write_profile_yaml(target, profile_path, target_verkey)
                cls._write_j_value_yaml(target, j_val_path)

                for sim in sims:
                    name_keys['sim_name'] = sim
                    sim_target_dir = NAME_FACTORY.sim_targetdir(**name_keys)
                    sim_profile_path = NAME_FACTORY.sim_profilefile(**name_keys)
                    sim_j_val_path = NAME_FACTORY.sim_j_valuefile(**name_keys)
                    if write_config:
                        cls._write_sim_target_config(target_config, target_dir, sim_target_dir)
                    cls._write_profile_yaml(target, sim_profile_path, target_verkey)
                    j_profile = cls._write_j_value_yaml(target, sim_j_val_path)
                    cls._write_sim_yaml(target, sim, sim_target_dir, target_verkey)

            roster_info_dict[roster_name] = tlist

        roster_file = os.path.join(ttype, 'roster_list.yaml')
        target_file = os.path.join(ttype, 'target_list.yaml')

        write_yaml(roster_info_dict, roster_file)
        write_yaml(target_info_dict, target_file)

        for sim in sims:
            sim_dir = os.path.join("%s_sim" % ttype, "sim_%s" % sim)
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
        self._write_target_dirs(args.ttype, roster_dict, base_config, sims)


def register_classes():
    PrepareTargets.register_class()

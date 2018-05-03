#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Scripts to run the all-sky diffuse analysis
"""
from __future__ import absolute_import, division, print_function

import sys
import os
import argparse

from collections import OrderedDict

from fermipy.jobs.job_archive import JobArchive, JobStatus, JobDetails
from fermipy.jobs.link import Link
from fermipy.jobs.chain import Chain, insert_app_config, purge_dict
from fermipy.utils import load_yaml

from dmpipe import defaults


def get_plot_config(plotting_dict, key, plot_channels_default=None, jpriors_default=None):
    sub_dict = plotting_dict.get(key, None)
    if sub_dict is None:
        return None
    channels = sub_dict.get('plot_channels', plot_channels_default)
    jpriors = sub_dict.get('jpriors', jpriors_default)
    return purge_dict(dict(channels=channels,
                           jpriors=jpriors))


def print_jobs(config_dict, indent="", stream=sys.stdout):

    for k, v in config_dict.items():
        appname = v.pop('appname', None)
        if appname is not None:
            stream.write("%s%s" % (indent, appname))
            for kk, vv in v.items():
                if isinstance(vv, list):
                    for vvv in vv:
                        stream.write(" --%s %s" % (kk, vvv))
                else:
                    stream.write(" --%s %s" % (kk, vv))
            stream.write("\n")
        else:
            stream.write("%s%s:\n" % (indent, k))
            print_jobs(v, indent + "  ", stream)


class PipelineData(Chain):
    """Small class to chain together the steps of the dSphs pipeline
    """
    appname = 'dmpipe-pipeline-data'
    linkname_default = 'pipeline-data'
    usage = '%s [options]' % (appname)
    description = 'Data analysis pipeline'

    default_options = dict(config=defaults.common['config'],
                           sim=defaults.sims['sim'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(PipelineData, self).__init__(linkname, **init_dict)

    def _register_link_classes(self):
        from fermipy.jobs.target_analysis import register_classes as register_analysis
        from fermipy.jobs.target_plotting import register_classes as register_plotting
        from dmpipe.dm_plotting import register_classes as register_plotting_dm
        from dmpipe.dm_spectral import register_classes as register_spectral
        register_plotting()
        register_plotting_dm()
        register_spectral()
        register_analysis()

    def _map_arguments(self, input_dict):
        """Map from the top-level arguments to the arguments provided to
        the indiviudal links """
        config_yaml = input_dict['config']
        o_dict = OrderedDict()
        config_dict = load_yaml(config_yaml)
        link_prefix = config_dict.get('link_prefix', '')
        ttype = config_dict.get('ttype')
        config_template = config_dict.get('config_template', None)
        config_localpath = config_dict.get('config_localpath', None)
        specfile = config_dict.get('specfile')
        rosterlist = config_dict.get('rosterlist')
        targetlist = config_dict.get('targetlist')
        jpriors = config_dict.get('jpriors')

        data_plotting = config_dict.get('data_plotting')
        plot_channels_default = config_dict.get('plot_channels', [])

        insert_app_config(o_dict, link_prefix + 'analyze-roi',
                          'fermipy-analyze-roi-sg',
                          ttype=ttype,
                          targetlist=targetlist,
                          config=config_localpath)
        insert_app_config(o_dict, link_prefix + 'analyze-sed',
                          'fermipy-analyze-sed-sg',
                          ttype=ttype,
                          targetlist=targetlist,
                          config=config_localpath)
        insert_app_config(o_dict, link_prefix + 'convert-castro',
                          'dmpipe-convert-castro-sg',
                          ttype=ttype,
                          jpriors=jpriors,
                          targetlist=targetlist,
                          config=config_localpath,
                          specfile=specfile)
        insert_app_config(o_dict, link_prefix + 'stack-likelihood',
                          'dmpipe-stack-likelihood-sg',
                          ttype=ttype,
                          jpriors=jpriors,
                          rosterlist=rosterlist)

        config_plot_castro = get_plot_config(data_plotting, 'plot-castro')
        if config_plot_castro is not None:
            insert_app_config(o_dict, link_prefix + 'plot-castro-sg',
                              'fermipy-plot-castro-sg',
                              ttype=ttype,
                              targetlist=targetlist,
                              **config_plot_castro)

        config_plot_dm = get_plot_config(data_plotting, 'plot-dm',
                                         plot_channels_default, jpriors)
        if config_plot_castro is not None:
            insert_app_config(o_dict, link_prefix + 'plot-dm-sg',
                              'dmpipe-plot-dm-sg',
                              ttype=ttype,
                              targetlist=targetlist,
                              **config_plot_dm)

        config_plot_limits = get_plot_config(data_plotting, 'plot-limits',
                                             plot_channels_default, jpriors)
        if config_plot_limits is not None:
            insert_app_config(o_dict, link_prefix + 'plot-limits-sg',
                              'dmpipe-plot-limits-sg',
                              ttype=ttype,
                              targetlist=targetlist,
                              **config_plot_limits)

        config_plot_stacked_dm = get_plot_config(data_plotting, 'plot-stacked-dm',
                                                 plot_channels_default, jpriors)
        if config_plot_stacked_dm is not None:
            insert_app_config(o_dict, link_prefix + 'plot-stacked-dm-sg',
                              'dmpipe-plot-stacked-dm-sg',
                              linkname='plot-stacked-dm',
                              ttype=ttype,
                              rosterlist=rosterlist,
                              **config_plot_stacked_dm)

        config_plot_stacked_limits = get_plot_config(data_plotting, 'plot-stacked-limits',
                                                     plot_channels_default, jpriors)
        if config_plot_stacked_limits is not None:
            insert_app_config(o_dict, link_prefix + 'plot-stacked-limits-sg',
                              'dmpipe-plot-stacked-limits-sg',
                              linkname='plot-stacked-limits',
                              ttype=ttype,
                              rosterlist=rosterlist,
                              **config_plot_stacked_limits)
        return o_dict


class PipelineSim(Chain):
    """Small class to chain together the steps of the dSphs pipeline
    """
    appname = 'dmpipe-pipeline-sim'
    linkname_default = 'pipeline-sim'
    usage = '%s [options]' % (appname)
    description = 'Run gtselect and gtbin together'

    default_options = dict(config=defaults.common['config'],
                           sim=defaults.sims['sim'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(PipelineSim, self).__init__(linkname, **init_dict)

    def _register_link_classes(self):
        from fermipy.jobs.target_analysis import register_classes as register_analysis
        from fermipy.jobs.target_collect import register_classes as register_collect
        from fermipy.jobs.target_sim import register_classes as register_sim
        from fermipy.jobs.target_plotting import register_classes as register_plotting
        from dmpipe.dm_plotting import register_classes as register_plotting_dm
        from dmpipe.dm_spectral import register_classes as register_spectral
        from dmpipe.dm_collect import register_classes as register_collect_dm
        register_analysis()
        register_collect()
        register_sim()
        register_plotting()
        register_plotting_dm()
        register_spectral()
        register_collect_dm()

    def _map_arguments(self, input_dict):
        """Map from the top-level arguments to the arguments provided to
        the indiviudal links """
        config_yaml = input_dict['config']
        o_dict = OrderedDict()
        config_dict = load_yaml(config_yaml)
        link_prefix = config_dict.get('link_prefix', '')

        sim_name = input_dict['sim']
        sim_dict = config_dict['sims'][sim_name]

        ttype = config_dict.get('ttype')
        config_template = config_dict.get('config_template', None)
        config_localpath = config_dict.get('config_localpath', None)
        specfile = config_dict.get('specfile')
        rosterlist = config_dict.get('rosterlist')
        targetlist = config_dict.get('targetlist')
        jpriors = config_dict.get('jpriors')

        enumbins = config_dict.get('enumbins', 12)
        sim_values = config_dict['sim_defaults']
        sim_values.update(sim_dict)
        seed = sim_values.get('seed', 0)
        nsims = sim_values.get('nsims', 20)

        sim_plotting = config_dict.get('sim_plotting')
        plot_channels_default = config_dict.get('plot_channels', [])

        insert_app_config(o_dict, link_prefix + 'copy-base-roi',
                          'fermipy-copy-base-roi-sg',
                          ttype=ttype,
                          targetlist=targetlist,
                          rosterlist=rosterlist,
                          sim=sim_name,
                          config=config_template)
        insert_app_config(o_dict, link_prefix + 'simulate-roi',
                          'fermipy-simulate-roi-sg',
                          ttype=ttype,
                          sim=sim_name,
                          targetlist=targetlist,
                          config=config_localpath,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'convert-castro',
                          'dmpipe-convert-castro-sg',
                          ttype=ttype,
                          sim=sim_name,
                          jpriors=jpriors,
                          targetlist=targetlist,
                          config=config_localpath,
                          specfile=specfile,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'stack-likelihood',
                          'dmpipe-stack-likelihood-sg',
                          ttype=ttype,
                          sim=sim_name,
                          jpriors=jpriors,
                          rosterlist=rosterlist,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'collect-sed',
                          'fermipy-collect-sed-sg',
                          ttype=ttype,
                          sim=sim_name,
                          enumbins=enumbins,
                          targetlist=targetlist,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'collect-limits',
                          'dmpipe-collect-limits-sg',
                          ttype=ttype,
                          sim=sim_name,
                          jpriors=jpriors,
                          targetlist=targetlist,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'collect-stacked-limits',
                          'dmpipe-collect-stacked-limits-sg',
                          ttype=ttype,
                          sim=sim_name,
                          jpriors=jpriors,
                          rosterlist=rosterlist,
                          seed=seed, nsims=nsims)

        config_plot_stacked_dm = get_plot_config(sim_plotting, 'plot-stacked-dm',
                                                 plot_channels_default, jpriors)
        if config_plot_stacked_dm is not None:
            insert_app_config(o_dict, link_prefix + 'plot-stacked-dm',
                              'dmpipe-plot-stacked-dm-sg',
                              ttype=ttype,
                              sim=sim_name,
                              rosterlist=rosterlist,
                              **config_plot_stacked_dm)

        config_plot_stacked_limits = get_plot_config(sim_plotting, 'plot-stacked-limits',
                                                     plot_channels_default, jpriors)
        if config_plot_stacked_limits is not None:
            insert_app_config(o_dict, link_prefix + 'plot-stacked-limits',
                              'dmpipe-plot-stacked-limits-sg',
                              ttype=ttype,
                              sim=sim_name,
                              rosterlist=rosterlist,
                              **config_plot_stacked_limits)

        return o_dict


class PipelineRandom(Chain):
    """Small class to chain together the steps of the dSphs pipeline
    """
    appname = 'dmpipe-pipeline-random'
    linkname_default = 'pipeline-random'
    usage = '%s [options]' % (appname)
    description = 'Data analysis pipeline for random directions'

    default_options = dict(config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(PipelineRandom, self).__init__(linkname, **init_dict)

    def _register_link_classes(self):
        from fermipy.jobs.factory import LinkFactory
        if self.appname in LinkFactory._class_dict:
            return

        from fermipy.jobs.target_analysis import register_classes as register_analysis
        from fermipy.jobs.target_collect import register_classes as register_collect
        from fermipy.jobs.target_sim import register_classes as register_sim
        from fermipy.jobs.target_plotting import register_classes as register_plotting
        from dmpipe.dm_plotting import register_classes as register_plotting_dm
        from dmpipe.dm_spectral import register_classes as register_spectral
        from dmpipe.dm_collect import register_classes as register_collect_dm
        register_analysis()
        register_collect()
        register_sim()
        register_plotting()
        register_plotting_dm()
        register_spectral()
        register_collect_dm()

    def _map_arguments(self, input_dict):
        """Map from the top-level arguments to the arguments provided to
        the indiviudal links """
        config_yaml = input_dict['config']
        o_dict = OrderedDict()
        config_dict = load_yaml(config_yaml)
        link_prefix = config_dict.get('link_prefix', '')

        ttype = config_dict.get('ttype')
        config_template = config_dict.get('config_template', None)
        config_localpath = config_dict.get('config_localpath', None)
        specfile = config_dict.get('specfile')
        rosterlist = config_dict.get('rosterlist')
        targetlist = config_dict.get('targetlist')
        jpriors = config_dict.get('jpriors')
        random = config_dict.get('random')
        enumbins = config_dict.get('enumbins', 12)

        rand_dirs = random.get('rand_dirs')
        sim_values = config_dict['sim_defaults']
        seed = sim_values.get('seed', 0)
        nsims = sim_values.get('nsims', 20)

        rand_plotting = config_dict.get('rand_plotting')
        plot_channels_default = config_dict.get('plot_channels', [])

        insert_app_config(o_dict, link_prefix + 'copy-base-roi',
                          'fermipy-copy-base-roi-sg',
                          ttype=ttype,
                          targetlist=targetlist,
                          rosterlist=rosterlist,
                          sim='random',
                          config=config_template)
        insert_app_config(o_dict, link_prefix + 'random-dir-gen',
                          'fermipy-random-dir-gen-sg',
                          ttype=ttype,
                          rand_config=rand_dirs,
                          targetlist=targetlist,
                          sim='random',
                          config=config_localpath)
        insert_app_config(o_dict, link_prefix + 'fermipy-sed',
                          'fermipy-analyze-sed-sg',
                          ttype=ttype,
                          sim='random',
                          skydirs='skydirs.yaml',
                          targetlist=targetlist,
                          config=config_localpath,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'convert-castro',
                          'dmpipe-convert-castro-sg',
                          ttype=ttype,
                          sim='random',
                          jpriors=jpriors,
                          targetlist=targetlist,
                          config=config_localpath,
                          specfile=specfile,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'stack-likelihood',
                          'dmpipe-stack-likelihood-sg',
                          ttype=ttype,
                          sim='random',
                          jpriors=jpriors,
                          rosterlist=rosterlist,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'collect-sed',
                          'fermipy-collect-sed-sg',
                          ttype=ttype,
                          sim='random',
                          enumbins=enumbins,
                          targetlist=targetlist,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'collect-limits',
                          'dmpipe-collect-limits-sg',
                          ttype=ttype,
                          sim='random',
                          jpriors=jpriors,
                          targetlist=targetlist,
                          seed=seed, nsims=nsims)
        insert_app_config(o_dict, link_prefix + 'collect-stacked-limits',
                          'dmpipe-collect-stacked-limits-sg',
                          ttype=ttype,
                          sim='random',
                          jpriors=jpriors,
                          rosterlist=rosterlist,
                          seed=seed, nsims=nsims)

        config_plot_stacked_dm = get_plot_config(rand_plotting, 'plot-stacked-dm',
                                                 plot_channels_default, jpriors)
        if config_plot_stacked_dm is not None:
            insert_app_config(o_dict, link_prefix + 'plot-stacked-dm',
                              'dmpipe-plot-stacked-dm-sg',
                              ttype=ttype,
                              sim='random',
                              rosterlist=rosterlist,
                              **config_plot_stacked_dm)

        config_plot_stacked_limits = get_plot_config(rand_plotting, 'plot-stacked-limits',
                                                     plot_channels_default, jpriors)
        if config_plot_stacked_limits is not None:
            insert_app_config(o_dict, link_prefix + 'plot-stacked-limits',
                              'dmpipe-plot-stacked-limits-sg',
                              ttype=ttype,
                              sim='random',
                              rosterlist=rosterlist,
                              **config_plot_stacked_limits)

        return o_dict


class Pipeline(Chain):
    """Small class to chain together the steps of the dSphs pipeline
    """
    appname = 'dmpipe-pipeline'
    linkname_default = 'pipeline'
    usage = '%s [options]' % (appname)
    description = 'Data analysis pipeline'

    default_options = dict(config=defaults.common['config'],
                           dry_run=defaults.common['dry_run'])

    def __init__(self, **kwargs):
        """C'tor
        """
        linkname, init_dict = self._init_dict(**kwargs)
        super(Pipeline, self).__init__(linkname, **init_dict)
        self._preconfigured = False

    def _register_link_classes(self):
        from dmpipe.dm_prepare import PrepareTargets
        from dmpipe.dm_spectral import SpecTable
        PrepareTargets.register_class()
        SpecTable.register_class()
        PipelineData.register_class()
        PipelineSim.register_class()
        PipelineRandom.register_class()

    def preconfigure(self, config_yaml):
        """ """
        o_dict = OrderedDict()
        config_dict = load_yaml(config_yaml)
        ttype = config_dict.get('ttype')
        self.link_prefix = "%s." % ttype
        config_template = config_dict.get('config_template', None)
        rosters = config_dict.get('rosters')
        sims = config_dict.get('sims', {})
        sim_names = []
        sim_names += sims.keys()
        if 'random' in config_dict:
            sim_names += ['random']

        insert_app_config(o_dict, 'prepare-targets',
                          'dmpipe-prepare-targets',
                          ttype=ttype,
                          rosters=rosters,
                          sims=sim_names,
                          config=config_template)
        self._arg_dict = o_dict
        self._load_arguments()
        link = self['prepare-targets']

        key = JobDetails.make_fullkey(link.full_linkname)
        if len(link.jobs) == 0:
            raise ValueError("No Jobs")
        link_status = link.check_job_status(key)
        if link_status == JobStatus.done:
            return
        elif link_status == JobStatus.failed:
            link.clean_jobs()
        link.run_with_log()

    def _map_arguments(self, input_dict):
        """Map from the top-level arguments to the arguments provided to
        the indiviudal links """

        config_yaml = input_dict['config']
        o_dict = OrderedDict()
        config_dict = load_yaml(config_yaml)
        ttype = config_dict.get('ttype')
        config_template = config_dict.get('config_template', None)
        rosters = config_dict.get('rosters')
        specfile = config_dict.get('specfile')
        sims = config_dict.get('sims', {})
        sim_names = []
        sim_names += sims.keys()
        if 'random' in config_dict:
            sim_names += ['random']
        dry_run = input_dict.get('dry_run', False)

        insert_app_config(o_dict, 'prepare-targets',
                          'dmpipe-prepare-targets',
                          ttype=ttype,
                          rosters=rosters,
                          sims=sim_names,
                          config=config_template)

        insert_app_config(o_dict, 'spec-table',
                          'dmpipe-spec-table',
                          ttype=ttype,
                          config=config_template,
                          specfile=specfile)

        insert_app_config(o_dict, 'data',
                          'dmpipe-pipeline-data',
                          linkname='data',
                          link_prefix='data.',
                          config=config_yaml,
                          dry_run=dry_run)

        for sim in sims.keys():
            linkname = 'sim_%s' % sim
            insert_app_config(o_dict, linkname,
                              'dmpipe-pipeline-sim',
                              link_prefix='%s.' % linkname,
                              linkname=linkname,
                              config=config_yaml,
                              sim=sim,
                              dry_run=dry_run)

        if 'random' in config_dict:
            insert_app_config(o_dict, 'random',
                              'dmpipe-pipeline-random',
                              link_prefix='random.',
                              linkname='random',
                              config=config_yaml,
                              dry_run=dry_run)

        return o_dict

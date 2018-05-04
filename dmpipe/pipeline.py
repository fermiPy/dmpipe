#!/usr/bin/env python

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Scripts to run the all-sky diffuse analysis
"""
from __future__ import absolute_import, division, print_function

from fermipy.utils import load_yaml

from fermipy.jobs.job_archive import JobStatus, JobDetails
from fermipy.jobs.chain import Chain, purge_dict
from fermipy.jobs.target_analysis import AnalyzeROI_SG, AnalyzeSED_SG
from fermipy.jobs.target_collect import CollectSED_SG
from fermipy.jobs.target_sim import CopyBaseROI_SG, SimulateROI_SG,\
    RandomDirGen_SG
from fermipy.jobs.target_plotting import PlotCastro_SG

from dmpipe.dm_spectral import ConvertCastro_SG, StackLikelihood_SG
from dmpipe.dm_collect import CollectLimits_SG, CollectStackedLimits_SG
from dmpipe.dm_plotting import PlotDM_SG, PlotLimits_SG,\
    PlotStackedDM_SG, PlotStackedLimits_SG
from dmpipe.dm_prepare import PrepareTargets
from dmpipe.dm_spectral import SpecTable

from dmpipe import defaults


def _get_plot_config(plotting_dict, key,
                     plot_channels_default=None, jpriors_default=None):
    """ Get the configuration for a type of plot """
    sub_dict = plotting_dict.get(key, None)
    if sub_dict is None:
        return None
    channels = sub_dict.get('plot_channels', plot_channels_default)
    jpriors = sub_dict.get('jpriors', jpriors_default)
    return purge_dict(dict(channels=channels,
                           jpriors=jpriors))


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


    def _map_arguments(self, input_dict):
        """Map from the top-level arguments to the arguments provided to
        the indiviudal links """
        config_yaml = input_dict['config']
        config_dict = load_yaml(config_yaml)
        ttype = config_dict.get('ttype')
        config_localpath = config_dict.get('config_localpath', None)
        specfile = config_dict.get('specfile')
        rosterlist = config_dict.get('rosterlist')
        targetlist = config_dict.get('targetlist')
        jpriors = config_dict.get('jpriors')

        data_plotting = config_dict.get('data_plotting')
        plot_channels_default = config_dict.get('plot_channels', [])

        self._load_link_args('analyze-roi',
                             AnalyzeROI_SG,
                             ttype=ttype,
                             targetlist=targetlist,
                             config=config_localpath)
        self._load_link_args('analyze-sed',
                             AnalyzeSED_SG,
                             ttype=ttype,
                             targetlist=targetlist,
                             config=config_localpath)
        self._load_link_args('convert-castro',
                             ConvertCastro_SG,
                             ttype=ttype,
                             jpriors=jpriors,
                             targetlist=targetlist,
                             config=config_localpath,
                             specfile=specfile)
        self._load_link_args('stack-likelihood',
                             StackLikelihood_SG,
                             ttype=ttype,
                             jpriors=jpriors,
                             rosterlist=rosterlist)

        config_plot_castro = _get_plot_config(data_plotting, 'plot-castro')
        if config_plot_castro is not None:
            self._load_link_args('plot-castro-sg',
                                 PlotCastro_SG,
                                 ttype=ttype,
                                 targetlist=targetlist,
                                 **config_plot_castro)

        config_plot_dm = _get_plot_config(data_plotting, 'plot-dm',
                                          plot_channels_default, jpriors)
        if config_plot_castro is not None:
            self._load_link_args('plot-dm-sg',
                                 PlotDM_SG,
                                 ttype=ttype,
                                 targetlist=targetlist,
                                 **config_plot_dm)

        config_plot_limits = _get_plot_config(data_plotting, 'plot-limits',
                                              plot_channels_default, jpriors)
        if config_plot_limits is not None:
            self._load_link_args('plot-limits-sg',
                                 PlotLimits_SG,
                                 ttype=ttype,
                                 targetlist=targetlist,
                                 **config_plot_limits)

        config_plot_stacked_dm = _get_plot_config(data_plotting, 'plot-stacked-dm',
                                                  plot_channels_default, jpriors)
        if config_plot_stacked_dm is not None:
            self._load_link_args('plot-stacked-dm-sg',
                                 PlotStackedDM_SG,
                                 ttype=ttype,
                                 rosterlist=rosterlist,
                                 **config_plot_stacked_dm)

        config_plot_stacked_limits = _get_plot_config(data_plotting, 'plot-stacked-limits',
                                                      plot_channels_default, jpriors)
        if config_plot_stacked_limits is not None:
            self._load_link_args('plot-stacked-limits-sg',
                                 PlotStackedLimits_SG,
                                 ttype=ttype,
                                 rosterlist=rosterlist,
                                 **config_plot_stacked_limits)


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

    def _map_arguments(self, input_dict):
        """Map from the top-level arguments to the arguments provided to
        the indiviudal links """
        config_yaml = input_dict['config']
        config_dict = load_yaml(config_yaml)

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

        self._load_link_args('copy-base-roi',
                             CopyBaseROI_SG,
                             ttype=ttype,
                             targetlist=targetlist,
                             rosterlist=rosterlist,
                             sim=sim_name,
                             config=config_template)
        self._load_link_args('simulate-roi',
                             SimulateROI_SG,
                             ttype=ttype,
                             sim=sim_name,
                             targetlist=targetlist,
                             config=config_localpath,
                             seed=seed, nsims=nsims)
        self._load_link_args('convert-castro',
                             ConvertCastro_SG,
                             ttype=ttype,
                             sim=sim_name,
                             jpriors=jpriors,
                             targetlist=targetlist,
                             config=config_localpath,
                             specfile=specfile,
                             seed=seed, nsims=nsims)
        self._load_link_args('stack-likelihood',
                             StackLikelihood_SG,
                             ttype=ttype,
                             sim=sim_name,
                             jpriors=jpriors,
                             rosterlist=rosterlist,
                             seed=seed, nsims=nsims)
        self._load_link_args('collect-sed',
                             CollectSED_SG,
                             ttype=ttype,
                             sim=sim_name,
                             enumbins=enumbins,
                             targetlist=targetlist,
                             seed=seed, nsims=nsims)
        self._load_link_args('collect-limits',
                             CollectLimits_SG,
                             ttype=ttype,
                             sim=sim_name,
                             jpriors=jpriors,
                             targetlist=targetlist,
                             seed=seed, nsims=nsims)
        self._load_link_args('collect-stacked-limits',
                             CollectStackedLimits_SG,
                             ttype=ttype,
                             sim=sim_name,
                             jpriors=jpriors,
                             rosterlist=rosterlist,
                             seed=seed, nsims=nsims)

        config_plot_stacked_dm = _get_plot_config(sim_plotting, 'plot-stacked-dm',
                                                  plot_channels_default, jpriors)
        if config_plot_stacked_dm is not None:
            self._load_link_args('plot-stacked-dm',
                                 PlotStackedDM_SG,
                                 ttype=ttype,
                                 sim=sim_name,
                                 rosterlist=rosterlist,
                                 **config_plot_stacked_dm)

        config_plot_stacked_limits = _get_plot_config(sim_plotting, 'plot-stacked-limits',
                                                      plot_channels_default, jpriors)
        if config_plot_stacked_limits is not None:
            self._load_link_args('plot-stacked-limits',
                                 PlotStackedLimits_SG,
                                 ttype=ttype,
                                 sim=sim_name,
                                 rosterlist=rosterlist,
                                 **config_plot_stacked_limits)


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


    def _map_arguments(self, input_dict):
        """Map from the top-level arguments to the arguments provided to
        the indiviudal links """
        config_yaml = input_dict['config']
        config_dict = load_yaml(config_yaml)

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

        self._load_link_args('copy-base-roi',
                             CopyBaseROI_SG,
                             ttype=ttype,
                             targetlist=targetlist,
                             rosterlist=rosterlist,
                             sim='random',
                             config=config_template)
        self._load_link_args('random-dir-gen',
                             RandomDirGen_SG,
                             ttype=ttype,
                             rand_config=rand_dirs,
                             targetlist=targetlist,
                             sim='random',
                             config=config_localpath)
        self._load_link_args('analyze-sed',
                             AnalyzeSED_SG,
                             ttype=ttype,
                             sim='random',
                             skydirs='skydirs.yaml',
                             targetlist=targetlist,
                             config=config_localpath,
                             seed=seed, nsims=nsims)
        self._load_link_args('convert-castro',
                             ConvertCastro_SG,
                             ttype=ttype,
                             sim='random',
                             jpriors=jpriors,
                             targetlist=targetlist,
                             config=config_localpath,
                             specfile=specfile,
                             seed=seed, nsims=nsims)
        self._load_link_args('stack-likelihood',
                             StackLikelihood_SG,
                             ttype=ttype,
                             sim='random',
                             jpriors=jpriors,
                             rosterlist=rosterlist,
                             seed=seed, nsims=nsims)
        self._load_link_args('collect-sed',
                             CollectSED_SG,
                             ttype=ttype,
                             sim='random',
                             enumbins=enumbins,
                             targetlist=targetlist,
                             seed=seed, nsims=nsims)
        self._load_link_args('collect-limits',
                             CollectLimits_SG,
                             ttype=ttype,
                             sim='random',
                             jpriors=jpriors,
                             targetlist=targetlist,
                             seed=seed, nsims=nsims)
        self._load_link_args('collect-stacked-limits',
                             CollectStackedLimits_SG,
                             ttype=ttype,
                             sim='random',
                             jpriors=jpriors,
                             rosterlist=rosterlist,
                             seed=seed, nsims=nsims)

        config_plot_stacked_dm = _get_plot_config(rand_plotting, 'plot-stacked-dm',
                                                  plot_channels_default, jpriors)
        if config_plot_stacked_dm is not None:
            self._load_link_args('plot-stacked-dm',
                                 PlotStackedDM_SG,
                                 ttype=ttype,
                                 sim='random',
                                 rosterlist=rosterlist,
                                 **config_plot_stacked_dm)

        config_plot_stacked_limits = _get_plot_config(rand_plotting, 'plot-stacked-limits',
                                                      plot_channels_default, jpriors)
        if config_plot_stacked_limits is not None:
            self._load_link_args('plot-stacked-limits',
                                 PlotStackedLimits_SG,
                                 ttype=ttype,
                                 sim='random',
                                 rosterlist=rosterlist,
                                 **config_plot_stacked_limits)



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


    def preconfigure(self, config_yaml):
        """ Run any links needed to build files
        that are used in _map_arguments """
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

        self._load_link_args('prepare-targets',
                             PrepareTargets,
                             ttype=ttype,
                             rosters=rosters,
                             sims=sim_names,
                             config=config_template)
        link = self['prepare-targets']

        key = JobDetails.make_fullkey(link.full_linkname)
        if not link.jobs:
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

        self._load_link_args('prepare-targets',
                             PrepareTargets,
                             ttype=ttype,
                             rosters=rosters,
                             sims=sim_names,
                             config=config_template)

        self._load_link_args('spec-table',
                             SpecTable,
                             ttype=ttype,
                             config=config_template,
                             specfile=specfile)

        self._load_link_args('data',
                             PipelineData,
                             link_prefix='data.',
                             config=config_yaml,
                             dry_run=dry_run)

        for sim in sims.keys():
            linkname = 'sim_%s' % sim
            self._load_link_args(linkname,
                                 PipelineSim,
                                 link_prefix='%s.' % linkname,
                                 config=config_yaml,
                                 sim=sim,
                                 dry_run=dry_run)

        if 'random' in config_dict:
            self._load_link_args('random',
                                 PipelineRandom,
                                 link_prefix='random.',
                                 config=config_yaml,
                                 dry_run=dry_run)

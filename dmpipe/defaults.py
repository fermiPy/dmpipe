# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Analysis defaults options for DM pipelien analysis
"""
from __future__ import absolute_import, division, print_function

# Options for diffuse analysis
generic = {
    'outfile': (None, 'Path to output file.', str),
    'infile': (None, 'Path to input file.', str),
    'limitfile': (None, 'Path to file with limits.', str),
    'summaryfile': (None, 'Path to file with results summaries.', str),
    }

common = {
    'ttype': (None, 'Type of target being analyzed.', str),
    'roster': (None, 'Name of a dmsky target roster.', str),
    'rosters': ([], 'Name of a dmsky target roster.', list),
    'target': (None, 'Name of analysis target.', str),
    'targetlist': (None, 'Path to the target list.', str),
    'rosterlist': (None, 'Path to the roster list.', str),
    'config': (None, 'Path to fermipy config file.', str),
    'specconfig': (None, 'Path to DM yaml file defining DM spectra of interest.', str),
    'specfile': (None, 'Path to DM spectrum file.', str),
    'sed_file': (None, 'Path to SED file.', str),
    'profile_file': (None, 'Path to yaml file with target profile', str),
    'j_value_file': (None, 'Path to yaml file with target j_value', str),
    'profiles' : ([], 'List of profiles to analyze', list),
    'jprior': (None, 'Types of Prior on J-factor', str),
    'jpriors': ([], 'Types of Prior on J-factor', list),
    'nsims': (-1, 'Number of simulations to run.', int),
    'channels': ([], 'DM annihilation channels', list),
    'chan': ('bb', 'DM annihilation channel', str),
    'mass': (100, 'DM particle mass', float),
    'spec_type' : ('eflux', 'Type of flux to consider', str),
    'dry_run' : (False, 'Print commands but do not run them.', bool),
    'clobber' : (False, 'Overwrite existing files.', bool),
    }

sims = {
    'sim': (None, 'Name of the simulation scenario.', str),
    'sims': ([], 'Names of the simulation scenario.', list),
    'nsims': (20, 'Number of simulations to run.', int),
    'seed': (0, 'Seed number for first simulation.', int),
    'skydirs': (None, 'Yaml file with blank sky directions.', str),
    'rand_config': (None, 'Path to config file for genaration random sky dirs', str),
   }

collect = {
    'bands': (None, 'Name of file with expected limit bands.', str),
    'write_full' : (False, 'Write file with full collected results', bool),
    'write_summary' : (False, 'Write file with summary of collected results', bool),
    }

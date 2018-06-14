# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Analysis defaults options for DM pipelien analysis
"""
from __future__ import absolute_import, division, print_function

from fermipy.jobs import defaults as base_defaults

generic = {
    'limitfile': (None, 'Path to file with limits.', str),
}
generic.update(base_defaults.generic)

common = {
    'roster': (None, 'Name of a dmsky target roster.', str),
    'rosters': ([], 'Name of a dmsky target roster.', list),
    'rosterlist': (None, 'Path to the roster list.', str),
    'alias_dict': (None, 'File to rename target version keys.', str),
    'specconfig': (None, 'Path to DM yaml file defining DM spectra of interest.', str),
    'specfile': (None, 'Path to DM spectrum file.', str),
    'astro_value_file': (None, 'Path to yaml file with target j_value or d_value', str),
    'astro_prior': (None, 'Types of Prior on J-factor or D-factor', str),
    'astro_priors': ([], 'Types of Prior on J-factor or D-factor', list),
    'spatial_models': ([], 'Types of spatial models to use', list),
    'channels': ([], 'DM annihilation channels', list),
    'chan': ('bb', 'DM annihilation channel', str),
    'mass': (100, 'DM particle mass', float),
    'spec_type': ('eflux', 'Type of flux to consider', str),
    'clobber': (False, 'Overwrite existing files.', bool),
}
common.update(base_defaults.common)

sims = {}
sims.update(base_defaults.sims)

collect = {
    'bands': (None, 'Name of file with expected limit bands.', str),
}
collect.update(base_defaults.collect)

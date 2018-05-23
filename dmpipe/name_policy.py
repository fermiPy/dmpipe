# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handle the naming conventions for DM pipeline analysis
"""
from __future__ import absolute_import, division, print_function

import sys

from fermipy.jobs.utils import is_null, is_not_null
from fermipy.jobs.name_policy import NameFactory as NameFactory_Base


class NameFactory(NameFactory_Base):
    """ Helper class to define file names and keys consistently. """

    # DM spectral file target configruation
    specconfig_format = 'config/dm_spectra_{target_type}.yaml'

    # DM spectral file
    specfile_format = 'dm_spectra_{target_type}.fits'

    # target keys, these are how we specify various files associated with
    # particular targets

    # Roster list file format
    rosterfile_format = '{target_type}/{rosterlist}'

    # Simulated rosterlist  file format
    sim_rosterfile_format = '{target_type}_sim/sim_{sim_name}/{rosterlist}'

    # Information about a particular target j factor
    j_valuefile_format = '{target_type}/{target_name}/j_val_{target_version}.yaml'

    # Information about a particular target j factor
    sim_j_valuefile_format = '{target_type}_sim/sim_{sim_name}/{target_name}/j_val_{target_version}.yaml'

    # DM likelilood file for a particular target (and j-factor prior)
    dmlikefile_format = '{target_type}/{target_name}/dmlike_{profile}_{jprior}.fits'

    # DM limits file for a particular target (and j-factor prior)
    dmlimitsfile_format = '{target_type}/{target_name}/dmlimits_{profile}_{jprior}.fits'

    # Stacked DM limits file for a particular roster (and j-factor prior)
    resultsfile_format = '{target_type}/stacked/results_{roster_name}_{jprior}.fits'

    # Stacked DM limits file for a particular roster (and j-factor prior)
    stackedlimitsfile_format = '{target_type}/stacked/limits_{roster_name}_{jprior}.fits'

    # Simulated DM likelilood file for a particular target (and j-factor prior)
    sim_dmlikefile_format = '{target_type}_sim/sim_{sim_name}/{target_name}/dmlike_{profile}_{jprior}_{seed}.fits'

    # Simulated DM limits file for a particular target (and j-factor prior)
    sim_dmlimitsfile_format = '{target_type}_sim/sim_{sim_name}/{target_name}/dmlimits_{profile}_{jprior}_{seed}.fits'

    # Simulated Stacked DM limits file for a particular roster (and j-factor
    # prior)
    sim_resultsfile_format = '{target_type}_sim/sim_{sim_name}/stacked/results_{roster_name}_{jprior}_{seed}.fits'

    # Stacked DM limits file for a particular roster (and j-factor prior)
    sim_stackedlimitsfile_format = '{target_type}_sim/sim_{sim_name}/stacked/limits_{roster_name}_{jprior}_{seed}.fits'

    def specconfig(self, **kwargs):
        """ return the name of the input configuration file
        """
        return self._format_from_dict(NameFactory.specconfig_format, **kwargs)

    def specfile(self, **kwargs):
        """ return the name of DM spectral file
        """
        return self._format_from_dict(NameFactory.specfile_format, **kwargs)

    def rosterfile(self, **kwargs):
        """ return the name for the Roster list file
        """
        return self._format_from_dict(NameFactory.rosterfile_format, **kwargs)

    def sim_rosterfile(self, **kwargs):
        """ return the name for the Roster list file for simulation
        """
        return self._format_from_dict(
            NameFactory.sim_rosterfile_format, **kwargs)

    def j_valuefile(self, **kwargs):
        """ return the name of the yaml file with information about a partiuclar target j factor
        """
        return self._format_from_dict(NameFactory.j_valuefile_format, **kwargs)

    def sim_j_valuefile(self, **kwargs):
        """ return the name of the yaml file with information about a partiuclar target j factor
        """
        return self._format_from_dict(
            NameFactory.sim_j_valuefile_format, **kwargs)

    def dmlikefile(self, **kwargs):
        """ return the name for the DM likelilood file for a particular target
        """
        return self._format_from_dict(NameFactory.dmlikefile_format, **kwargs)

    def dmlimitsfile(self, **kwargs):
        """ return the name for the DM limits file for a particular target
        """
        return self._format_from_dict(
            NameFactory.dmlimitsfile_format, **kwargs)

    def resultsfile(self, **kwargs):
        """ return the name for the stacked results file for a particular roster
        """
        return self._format_from_dict(NameFactory.resultsfile_format, **kwargs)

    def stackedlimitsfile(self, **kwargs):
        """ return the name for the stacked limits file for a particular roster
        """
        return self._format_from_dict(
            NameFactory.stackedlimitsfile_format, **kwargs)

    def sim_dmlikefile(self, **kwargs):
        """ return the name for the simulated DM likelilood file for a particular target
        """
        if 'seed' not in kwargs:
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(
            NameFactory.sim_dmlikefile_format, **kwargs)

    def sim_dmlimitsfile(self, **kwargs):
        """ return the name for the simulated DM limits file for a particular target
        """
        if 'seed' not in kwargs:
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(
            NameFactory.sim_dmlimitsfile_format, **kwargs)

    def sim_resultsfile(self, **kwargs):
        """ return the name for the stacked results file for a particular roster
        """
        if 'seed' not in kwargs:
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(
            NameFactory.sim_resultsfile_format, **kwargs)

    def sim_stackedlimitsfile(self, **kwargs):
        """ return the name for the stacked limits file for a particular roster
        """
        if 'seed' not in kwargs:
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(
            NameFactory.sim_stackedlimitsfile_format, **kwargs)

    def resolve_rosterfile(self, args, require_sim_name=False):
        """Get the name of the roster based on the job arguments"""
        ttype = args.get('ttype')
        if is_null(ttype):
            sys.stderr.write('Target type must be specified')
            return (None, None)

        sim = args.get('sim')
        if is_null(sim):
            if require_sim_name:
                sys.stderr.write('Simulation scenario must be specified')
                return (None, None)
            else:
                sim = None

        name_keys = dict(target_type=ttype,
                         rosterlist='roster_list.yaml',
                         sim_name=sim,
                         fullpath=True)
        if sim is None:
            rosterfile = self.rosterfile(**name_keys)
        else:
            rosterfile = self.sim_rosterfile(**name_keys)

        roster_override = args.get('rosterfile')
        if is_not_null(roster_override):
            rosterfile = roster_override

        return (rosterfile, sim)

    def resolve_specfile(self, args):
        """Get the name of the specturm file based on the job arguments"""
        ttype = args.get('ttype')
        if is_null(ttype):
            sys.stderr.write('Target type must be specified')
            return None
        name_keys = dict(target_type=ttype,
                         fullpath=True)
        specfile = self.specfile(**name_keys)
        spec_override = args.get('specfile')
        if is_not_null(spec_override):
            specfile = spec_override
        return specfile

    def resolve_specconfig(self, args):
        """Get the name of the specturm file based on the job arguments"""
        ttype = args.get('ttype')
        if is_null(ttype):
            sys.stderr.write('Target type must be specified')
            return None
        name_keys = dict(target_type=ttype,
                         fullpath=True)
        specconfig = self.specconfig(**name_keys)
        spec_override = args.get('specconfig')
        if is_not_null(spec_override):
            specconfig = spec_override
        return specconfig

# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Handle the naming conventions for DM pipeline analysis
"""
from __future__ import absolute_import, division, print_function

import yaml


class NameFactory(object):
    """ Helper class to define file names and keys consistently. """

    # Input configuration file
    ttypeconfig_format = 'config/config_{target_type}.yaml'

    # DM spectral file target configruation
    specconfig_format = 'config/dm_spectra_{target_type}.yaml'

    # DM spectral file
    specfile_format = 'dm_spectra_{target_type}.fits'

    # target keys, these are how we specify various files associated with particular targets

    # Directory for a particular target
    targetdir_format = '{target_type}/{target_name}'

    # Directory for simulations for a particular target
    sim_targetdir_format = '{target_type}_sim/sim_{sim_name}/{target_name}'

    # Roster list file format
    rosterfile_format = '{target_type}/{rosterlist}'

    # Simulated rosterlist  file format
    sim_rosterfile_format = '{target_type}_sim/sim_{sim_name}/{rosterlist}'

    # Targetlist file format
    targetfile_format = '{target_type}/{targetlist}'

    # Roster file format
    sim_targetfile_format = '{target_type}_sim/sim_{sim_name}/{targetlist}'

    # Information about a particular target profile
    profilefile_format = '{target_type}/{target_name}/profile_{profile}.yaml'

    # Information about a particular target profile
    sim_profilefile_format = '{target_type}_sim/sim_{sim_name}/{target_name}/profile_{profile}.yaml'

    # SED file for a particular target
    sedfile_format = '{target_type}/{target_name}/sed_{profile}.fits'

    # DM likelilood file for a particular target (and j-factor prior)
    dmlikefile_format = '{target_type}/{target_name}/dmlike_{profile}_{jprior}.fits'

    # DM limits file for a particular target (and j-factor prior)
    dmlimitsfile_format = '{target_type}/{target_name}/dmlimits_{profile}_{jprior}.fits'

    # Stacked DM limits file for a particular roster (and j-factor prior)
    resultsfile_format = '{target_type}/stacked/results_{roster_name}_{jprior}.fits'

    # Stacked DM limits file for a particular roster (and j-factor prior)
    stackedlimitsfile_format = '{target_type}/stacked/limits_{roster_name}_{jprior}.fits'

    # Simulated SED file for a particular target
    sim_sedfile_format = '{target_type}_sim/sim_{sim_name}/{target_name}/sed_{profile}_{seed}.fits'

    # Simulated DM likelilood file for a particular target (and j-factor prior)
    sim_dmlikefile_format = '{target_type}_sim/sim_{sim_name}/{target_name}/dmlike_{profile}_{jprior}_{seed}.fits'

    # Simulated DM limits file for a particular target (and j-factor prior)
    sim_dmlimitsfile_format = '{target_type}_sim/sim_{sim_name}/{target_name}/dmlimits_{profile}_{jprior}_{seed}.fits'

    # Simulated Stacked DM limits file for a particular roster (and j-factor prior)
    sim_resultsfile_format = '{target_type}_sim/sim_{sim_name}/stacked/results_{roster_name}_{jprior}_{seed}.fits'

    # Stacked DM limits file for a particular roster (and j-factor prior)
    sim_stackedlimitsfile_format = '{target_type}_sim/sim_{sim_name}/stacked/limits_{roster_name}_{jprior}_{seed}.fits'

    # Stamp files from scatter gather jobs
    stamp_format = 'stamps/{linkname}.stamp'

    # Full filepath
    fullpath_format = '{basedir}/{localpath}'

    def __init__(self, **kwargs):
        """ C'tor.  Set baseline dictionary used to resolve names
        """
        self.base_dict = kwargs.copy()

    def update_base_dict(self, yamlfile):
        """ Update the values in baseline dictionary used to resolve names
        """
        self.base_dict.update(**yaml.safe_load(open(yamlfile)))

    def _format_from_dict(self, format_string, **kwargs):
        """ Return a formatted file name dictionary components """
        kwargs_copy = self.base_dict.copy()
        kwargs_copy.update(**kwargs)
        localpath = format_string.format(**kwargs_copy)
        if kwargs.get('fullpath', False):
            return self.fullpath(localpath=localpath)
        else:
            return localpath

    def ttypeconfig(self, **kwargs):
        """ return the name of the input configuration file
        """
        return self._format_from_dict(NameFactory.ttypeconfig_format, **kwargs)

    def specconfig(self, **kwargs):
        """ return the name of the input configuration file
        """
        return self._format_from_dict(NameFactory.specconfig_format, **kwargs)
    
    def specfile(self, **kwargs):
        """ return the name of DM spectral file
        """
        return self._format_from_dict(NameFactory.specfile_format, **kwargs)
      
    def targetdir(self, **kwargs):
        """ return the name for the directory for a particular target
        """
        return self._format_from_dict(NameFactory.targetdir_format, **kwargs)
    
    def sim_targetdir(self, **kwargs):
        """ return the name for the directory for a particular target
        """
        return self._format_from_dict(NameFactory.sim_targetdir_format, **kwargs)

    def rosterfile(self, **kwargs):
        """ return the name for the Roster list file
        """
        return self._format_from_dict(NameFactory.rosterfile_format, **kwargs)
    
    def sim_rosterfile(self, **kwargs):
        """ return the name for the Roster list file for simulation
        """
        return self._format_from_dict(NameFactory.sim_rosterfile_format, **kwargs)
    
    def targetfile(self, **kwargs):
        """ return the name for the Target list file
        """
        return self._format_from_dict(NameFactory.targetfile_format, **kwargs)
    
    def sim_targetfile(self, **kwargs):
        """ return the name for the Target list file for simulation
        """
        return self._format_from_dict(NameFactory.sim_targetfile_format, **kwargs)
    
    def profilefile(self, **kwargs):
        """ return the name of the yaml file with information about a partiuclar profile
        """
        return self._format_from_dict(NameFactory.profilefile_format, **kwargs)

    def sim_profilefile(self, **kwargs):
        """ return the name of the yaml file with information about a partiuclar profile
        """
        return self._format_from_dict(NameFactory.sim_profilefile_format, **kwargs)

    def sedfile(self, **kwargs):
        """ return the name for the SED file for a particular target
        """
        return self._format_from_dict(NameFactory.sedfile_format, **kwargs)
    
    def dmlikefile(self, **kwargs):
        """ return the name for the DM likelilood file for a particular target
        """
        return self._format_from_dict(NameFactory.dmlikefile_format, **kwargs)
    
    def dmlimitsfile(self, **kwargs):
        """ return the name for the DM limits file for a particular target
        """
        return self._format_from_dict(NameFactory.dmlimitsfile_format, **kwargs)
 
    def resultsfile(self, **kwargs):
        """ return the name for the stacked results file for a particular roster
        """
        return self._format_from_dict(NameFactory.resultsfile_format, **kwargs)

    def stackedlimitsfile(self, **kwargs):
        """ return the name for the stacked limits file for a particular roster
        """
        return self._format_from_dict(NameFactory.stackedlimitsfile_format, **kwargs)
    
    def sim_sedfile(self, **kwargs):
        """ return the name for the simulated SED file for a particular target
        """
        if not kwargs.has_key('seed'):
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(NameFactory.sim_sedfile_format, **kwargs)
    
    def sim_dmlikefile(self, **kwargs):
        """ return the name for the simulated DM likelilood file for a particular target
        """
        if not kwargs.has_key('seed'):
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(NameFactory.sim_dmlikefile_format, **kwargs)
    
    def sim_dmlimitsfile(self, **kwargs):
        """ return the name for the simulated DM limits file for a particular target
        """
        if not kwargs.has_key('seed'):
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(NameFactory.sim_dmlimitsfile_format, **kwargs)

    def sim_resultsfile(self, **kwargs):
        """ return the name for the stacked results file for a particular roster
        """
        if not kwargs.has_key('seed'):
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(NameFactory.sim_resultsfile_format, **kwargs)
 
    def sim_stackedlimitsfile(self, **kwargs):
        """ return the name for the stacked limits file for a particular roster
        """
        if not kwargs.has_key('seed'):
            kwargs['seed'] = 'SEED'
        return self._format_from_dict(NameFactory.sim_stackedlimitsfile_format, **kwargs)
   
    def stamp(self, **kwargs):
        """Return the path for a stamp file for a scatter gather job"""
        kwargs_copy = self.base_dict.copy()
        kwargs_copy.update(**kwargs)
        return NameFactory.stamp_format.format(**kwargs_copy)

    def fullpath(self, **kwargs):
        """Return a full path name for a given file
        """
        kwargs_copy = self.base_dict.copy()
        kwargs_copy.update(**kwargs)
        return NameFactory.fullpath_format.format(**kwargs_copy)

    
    def resolve_targetfile(self, args, require_sim_name=False):
        """Get the name of the targetfile based on the job arguments"""
        ttype = args.get('ttype')
        if ttype in [None, 'none', 'None']:
            sys.stderr.write('Target type must be specified')
            return (None, None)
      
        sim = args.get('sim')
        if sim in [None, 'none', 'None']:
            if require_sim_name:
                sys.stderr.write('Simulation scenario must be specified')
                return (None, None)
            else:
                sim = None

        name_keys = dict(target_type=ttype,
                         targetlist='target_list.yaml',
                         sim_name=sim,
                         fullpath=True)
        if sim is None:
            targetfile = self.targetfile(**name_keys)
        else:
            targetfile = self.sim_targetfile(**name_keys)

        targets_override = args.get('targetfile')
        if targets_override not in [None, 'none', 'None']:
            targetfile = targets_override
        
        return (targetfile, sim)
  
    def resolve_rosterfile(self, args, require_sim_name=False):
        """Get the name of the roster based on the job arguments"""
        ttype = args.get('ttype')
        if ttype in [None, 'none', 'None']:
            sys.stderr.write('Target type must be specified')
            return (None, None)
      
        sim = args.get('sim')
        if sim in [None, 'none', 'None']:
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
        if roster_override not in [None, 'none', 'None']:
           rosterfile  = roster_override
        
        return (rosterfile, sim)
  
    def resolve_specfile(self, args, require_sim_name=False):
        """Get the name of the specturm file based on the job arguments"""
        ttype = args.get('ttype')
        if ttype in [None, 'none', 'None']:
            sys.stderr.write('Target type must be specified')
            return None
        name_keys = dict(target_type=ttype,
                         fullpath=True)
        specfile = self.specfile(**name_keys)
        spec_override = args.get('specfile')
        if spec_override not in [None, 'none', 'None']:
            specfile = spec_override
        return specfile

    def resolve_specconfig(self, args):
        """Get the name of the specturm file based on the job arguments"""
        ttype = args.get('ttype')
        if ttype in [None, 'none', 'None']:
            sys.stderr.write('Target type must be specified')
            return None
        name_keys = dict(target_type=ttype,
                         fullpath=True)
        specconfig = self.specconfig(**name_keys)
        spec_override = args.get('specconfig')
        if spec_override not in [None, 'none', 'None']:
            specconfig = spec_override
        return specconfig

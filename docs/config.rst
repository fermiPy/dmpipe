.. _config:

Configuration
=============

This page describes the configuration management scheme used within
the dmpipe package and documents the configuration parameters
that can be set in the configuration file.


##################################
Class Configuration
##################################

Analysis classes in the dmpipe package all inherit from the `fermipy.jobs.Link`
class, which allow user to invoke the class either interactively within python
or from the unix command line.

From the command line

.. code-block:: bash

   $ dmpipe-plot-dm --infile dSphs/segue_1/dmlike_ack2016_point_none.fits --chan bb --outfile dSphs/segue_1/dmlike_ack2016_point_none_bb.fits


From python there are a number of ways to do it, we recommend this:

.. code-block:: python

   from dmipe.dm_plotting import PlotDM
   link = PlotDM()
   link.update_args(dict(infile='dSphs/segue_1/dmlike_ack2016_point_none.fits',
                                     chan='bb', outfile='dSphs/segue_1/dmlike_ack2016_point_none_bb.fits'))
   link.run()				     

   

##################################
Configuration File
##################################

dmpipe uses `YAML <http://yaml.org/>`_ files to read and write its
configuration in a persistent format.  The configuration file has a
hierarchical structure that groups parameters into dictionaries that
are keyed to a section name (*data*, *binning*, etc.).

.. code-block:: yaml
   :caption: Sample Configuration

   ttype: dSphs
   rosters : ['test']
   jpriors : ['none', 'lgauss']
   spatial_models: ['point']
   alias_dict : 'config/aliases_dSphs.yaml'
   
   sim_defaults:
       seed : 0
       nsims : 20
       profile : ack2016

   sims:
       'null' : {}
       '100GeV_bb_1e25' : {}
       '100GeV_bb_1e26' : {}

   random: {}

   plot_channels = ['bb', 'tautau']

   data_plotting:
       plot-castro : {}
       plot-dm : {}
       plot-limits : {}
       plot-stacked-dm : {}
       plot-stacked-limits : {}

   sim_plotting:
       plot-stacked-dm : {}
       plot-stacked-limits : {}
       plot-control-limits : {}

   rand_plotting:
       plot-stacked-dm : {}
       plot-stacked-limits : {}


     
.. _config_top:

Top level configuration
-----------------------

Options at the top level apply to all parts of the analysis pipeline

.. code-block:: yaml
   :caption: Sample *top level* Configuration
                
     # Top level
     ttype : 'dSphs'     
     rosters : ['test']  
     jpriors : ['none', 'lgauss'] 
     spatial_models: ['point']     
     alias_dict : 'config/aliases_dSphs.yaml'

* ttype: str
  Target tpye.  This is used for bookkeeping mainly, to give the naem of the top-level directory, and to
  call out specfic configuration files.

* rosters: list
  List of `dmsky` rosters to analyze.   Each roster represents a self-consistent set of targets and DM models for each target.

* jpriors : list
  List of types of J-factor prior to use.

* spatial_models: : list
  List of types of spatial model to use when fitting the DM.  Options are
  * point : A point source
  * map: A spatial map (in a FITS file)
  * radial: A radial profile (in at text file) and central direction

* alias_dict : Filename [Optional]
  Path to a file that give short names for the DM model to use with each target.
  
.. note::  
  If multiple rosters include the same target and DM model, that target will only be analyzed once,
  and those results will be re-used when combining each roster.   


  
.. _config_sims:

Simulation configuration
------------------------

The *sim_defaults*, *sims* and *random* sections can be used to define
analysis configurations for control studies with simulations and
random sky directions.

.. code-block:: yaml
   :caption: Sample *simulation* Configuration

   sim_defaults:
       seed : 0
       nsims : 20
       profile : ack2016

   sims:
       'null' : {}
       '100GeV_bb_1e25' : {}
       '100GeV_bb_1e26' : {}

   random: {}

   
* sim_defaults : dict
  This is a dictionary of the parameters to use for simulations.
  This can be overridden for specific type of simulation.

  * seed : int
     Random number seed to use for the first simulation

  * nsims : int
     Number of simulations

  * profile : str
     Name of the DM spatial profile to use for simulations.  This must match a profile defined in the roster for each target.
     The 'alias_dict' file can be used to remap longer profile names, or to define a common name for all the profiles in a roster.

     
* sims : dict
  This is a dictionary of the simulation scenarious to consider, and
  of any option overrides for some of those scenarios.

  Each defined simulation needs a 'config/sim_{sim_name}.yaml' to define the injected source to use for that simulation.

* random: dict
  This is a dictionary of the options to use for random sky direction control studies.


_config_plotting

Plotting configuration
----------------------

.. code-block:: yaml
   :caption: Sample *plotting* Configuration


   plot_channels = ['bb', 'tautau']

   data_plotting:
       plot-castro : {}
       plot-dm : {}
       plot-limits : {}
       plot-stacked-dm : {}
       plot-stacked-limits : {}

   sim_plotting:
       plot-stacked-dm : {}
       plot-stacked-limits : {}
       plot-control-limits : {}

   rand_plotting:
       plot-stacked-dm : {}
       plot-stacked-limits : {}

  

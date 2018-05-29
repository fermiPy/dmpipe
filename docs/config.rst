.. _config:

Configuration
=============

This page describes the configuration management scheme used within
the dmpipe package and documents the configuration parameters
that can be set in the configuration file.

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

   


Master Configuration File
=========================

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

  


Additional Configuration files
==============================

In addition to the master configuration file, the pipeline needs a few additional files.


Fermipy Analysis Configuration Yaml
-----------------------------------

This is simply the a template of the `fermipy` configuration file to be used for the baseline analysis and SED fitting
in each ROI.  Details of the syntax and options are `here <https://fermipy.readthedocs.io/en/latest/config.html>` _
The actually direction and name of the target source in this file will be over written for each target.


Dark Matter Spectral Configuration Yaml
---------------------------------------

This file specifies the masses and channels to analyze the DM spectra for.  Here is an example of
this file:

.. code-block yaml

  # This is the list of channels we will analyze.
  # These must match channel names in the DMFitFunction class
  channels : ['ee', 'mumu', 'tautau', 'bb', 'tt', 'gg', 'ww', 'zz', 'cc', 'uu', 'dd', 'ss']

  # This defines the array of mass points we use (in GeV)
  # The points are sampled in log-space
  masses : 
    mass_min : 10.
    mass_max : 10000.
    mass_nstep : 13


Simulation Scenario Configuration Yaml
--------------------------------------

This file specifies the DM signal to inject in the analysis (if any).  Here is a example, note
that everything inside the 'injected_source' tag is in the format that `fermipy` expects to see
source defintions.

  .. code-block yaml

    # For positive control tests we with injected source.
    # In this case it is a DM annihilation spectrum.
    injected_source:
      name : dm
      source_model :
        SpatialModel : PointSource
        SpectrumType : DMFitFunction
        norm : 
          value : nan # This is the J-factor and depend on the target
        sigmav : 
          value: 1.0E-25 # cm^3 s^-1.  (i.e., very large cross section)
        mass : 
          value: 100.0 # GeV
        channel0 : 
          value : 4 # annihilation to b-quarks


For null simulations, you should include the 'injected_source' tag, but leave it blank
        
  .. code-block yaml

  # For positive control tests we with injected source.
  # In this case it is a DM annihilation spectrum.
  injected_source:


  
Profile Alias Configuration Yaml
--------------------------------

This is a small file that remaps the target profile names used by dmsky to shorter names (without
underscores in them).  Removing the underscores helps keep the file name fields more logical, and
dmpipe generally uses underscores as a field seperator.  This also keeps file names shorter, and allow
us to use roster with a mixed set of profile version to do simulations.  Here is an example:

  .. code-block yaml

  ackermann2016_photoj_0.6_nfw : ack2016
  geringer-sameth2015_nfw : gs2015


  
Random Direction Control Sample Configuration Yaml
--------------------------------------------------

The file define how we select random directions for the random direction control studies.  Here is an example:

  .. code-block yaml

    # These are the parameters for the random direction selection
    # The algorithm picks points on a grid 

    # File key for the first direction
    seed : 0
    # Number of directions to select
    nsims : 20

    # Step size between grid points (in deg)
    step_x : 1.0
    step_y : 1.0
    # Max distance from ROI center (in deg)
    max_x : 3.0
    max_y : 3.0





       

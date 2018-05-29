.. _output:

Output Files
============

Fermipy ROI Snapshots
---------------------

For each target, the pipeline will perform a baseline fit to the region
of interest (ROI) and produce a snapshot in the 
standard `fermipy ROI FITs file format <https://fermipy.readthedocs.io/en/latest/output.html>`_
These are producted by the `AnalyzeROI` link for data, and copied to directories used for simulations.



Fermipy SED Files
-----------------

For each target, and for each spatial profile used to model the target the pipeline will produce a
standard `fermipy SED FITs file format <https://fermipy.readthedocs.io/en/latest/advanced/sed.html#sed-fits-file>`_
These are producted by the `AnalyzeSED` link for data, or the `SimulateROI` link for simulations.


Dark Matter Likeihood 'Castro' Files
------------------------------------

For each target, spatial profile and J-factor prior combination, the pipeline will produce `DMCastroData` fits
file with the likelhoods as a function of DM interaction rate.   These are produced by the `ConvertCastro` link for
the individual targets, and by the `StackLikelihood` link for the stacked roster results.


Dark Matter Limits Files
------------------------

For each target, spatial profile and J-factor prior combination, the pipeline will produce `DMCastroData` fits
file with the upper limits on the DM interaction rate.   These will also be produced for the stacked results from
each roster and J-factor prior combination.


Expectation Band Files
----------------------

For every target, simulation scenario and J-factor prior combination, the pipeline will also produce FITs files with summeries of the limits to capture the expected limits bands.  




Bookkeeping and Generated Configuration Files
=============================================

Several files needed for bookeeping are created by the `PrepareTargets` script.


* Target List Yaml Files

  These are dictionaries of all the targets and all the profiles to consider for each target.  Here is an example from the
  simple test analysis:

  .. code-block:: yaml

    draco: [ack2016_point]
    segue_1: [ack2016_point]


* Roster List Yaml Files

   These are dictionaries of all the targets and target version that define each roster.  Here is an example from teh
   simple test analysis:

   test_point: ['segue_1:ack2016_point', 'draco:ack2016_point']
   

* ROI configuration Yaml Files

This is simply the `fermipy` configuration file to be used for the baseline analysis and SED fitting
in each ROI.  Details of the syntax and options are `here <https://fermipy.readthedocs.io/en/latest/config.html>` _
These are copied from the template version to each of the analysis directories and updated to include the
target name and direction.


* Spatial Profile Yaml Files

These file define the various spatial profiles used to fit each target.  The syntax is basically what
`fermipy` needs to create a new source.

  .. code-block:: yaml

    name: ack2016_point
    source_model: {DEC: 57.91528, RA: 260.05167, SpatialModel: PointSource, SpectrumType: PowerLaw}


* J-value Yaml Files

These file define the values of the J-factor for different profiles.  They are needed to convert the analysis results to DM annihilation rate.  Here is an example:

  .. code-block:: yaml

     {j_integ: 2.188e+18, j_sigma: 0.6, type: NFW}

* Simulation Input Yaml Files

These file define the various spatial profiles used to fit each target.  The syntax inside the 'injected_source' tag is
exactly what `fermipy` needs to create a new source.

  .. code-block:: yaml

    injected_source:
      name: dm
      source_model:
        SpatialModel: PointSource
        SpectrumType: DMFitFunction
        channel0: {value: 4}
        mass: {value: 100.0}
        norm: {value: 2.188e+18}
        sigmav: {value: 3.0e-26}


* Simulated Source Spectrum Yaml Files

  These files are created by the `SimulateROI` task, and contain some information about the simulated sources.
  

* Source Correlation Yaml Files

  These file are created by the `AnalyzeSED` (for data) or `SimulateROI` (for simulations) tasks, and contain the
  correlation factors between the target source and any other source in the ROI above the threshold for special
  treatment (typically 0.25).


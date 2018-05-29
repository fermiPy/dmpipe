.. dmpipe documentation master file, created by
   sphinx-quickstart on Fri Mar 20 15:40:47 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dmpipe's documentation!
===================================

##################################
Introduction
##################################

This is the dmpipe documentation page. 

dmpipe is a python package that implements
an analysis pipeline to search for Dark Matter (DM) signals
in the data from the Large Area Telescope (LAT).

For more information about the Fermi mission and the LAT instrument please
refer to the `Fermi Science Support Center <http://fermi.gsfc.nasa.gov/ssc/>`_.

The dmpipe package is built on the pyLikelihood interface of the
Fermi Science Tools and the `Fermipy <https://fermipy.readthedocs.io/en/latest/>`_
and `dmksy <https://dmsky.readthedocs.io/en/latest/>`_
software packages.


dmpipe implements a standard analysis pipeline that:

* Uses the dmsky package to define lists of analysis targets.

* Models the DM gamma-ray spectra in a standard set of annihilation channels.
    
* Performs standard source analyses in regions of interest (ROI) around each target.

* Extracts Spectal Energy Denisty (SED) likelihood information for each target,
  and possibly for multiple dark matter spatial profiles for each target.

* Converts the SED likelihood information to likelihood information in
  on the DM interaction rate (i.e., the thermally averaged cross
  section), accounting for uncertainties on the DM spatial profile of
  the each target.

* Combines the results for multiple targets by performing likelihood statcking. 

* Implements control versions of the analysis pipeline, using both
  simulated data, and randomly selected control directiions within
  each target's ROI.

dmpipe uses a configuration-file driven workflow in which the
analysis parameters (data selection, IRFs, and ROI model) are defined
in a small set of YAML configuration files.  Analysis is executed through a python
script that dispatchs small analysis jobs to the computer batch farm.

For instructions on installing dmpipe see the :ref:`install` page.
For a short introduction to using dmipe see the :ref:`quickstart`.

Getting Help
------------

If you have questions about using dmpipe please open a `GitHub Issue
<https://github.com/fermiPy/dmpipe/issues>`_.


Documentation Contents
----------------------

.. toctree::
   :includehidden:
   :maxdepth: 3

   install
   quickstart
   config
   output
   dmpipe

Indices and tables
==================

.. * :ref:`genindex`
.. * :ref:`modindex`
.. * :ref:`search`


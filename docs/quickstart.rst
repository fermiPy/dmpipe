.. _quickstart:

Quickstart Guide
================

This page walks through the steps to setup and perform a basic
analysis of two dark matter targets, the "Segue 1" and "Draco" dwarf spheroidal galaxies.

Getting input data
------------------

First, you will need get Fermi-LAT data to analyze.  In particular
you will need:

* An event list, a so-called 'FT1' file with event data.
* The correspond spacecraft pointing history file: the so-called 'FT2' file.
* A 'livetime' htype cube of the amout of time each direction in the sky was at a particular direction with respect to the LAT instrument boresight.




Running the example notebook
----------------------------

First, download the configuration files and python notebook for this analysis

.. code-block:: bash

   $ curl -OL https://raw.githubusercontent.com/fermiPy/fermipy-extras/master/data/dmpipe_example.tar.gz
   $ tar zxvf dmpipe_example.tar.gz
   $ cd dmpipe_example
   $ curl -OL https://raw.githubusercontent.com/fermiPy/fermipy-extras/master/notebooks/dSphs.ipynb

Now you need to do two things to set up to run the example notebook

* Point the dmksy package at the target "Roster" you just downloaded.
  .. code-block:: bash

    $ export DMSKY_PATH=<current_dir>

* Edit the 'config/config_dSphs.yaml' file so that the 'evfile', 'ltcube', and 'scfile' lines refer to the input data you set up above.

Now you can run the example notebook.

.. code-block:: bash

   $ jupyter-notebook dSphs.ipynb

   
Running the analysis
--------------------

Every step of the analysis, including the top-level script that runs the entire analysis,
can be invoked directly from the UNIX command line:

.. code-block:: bash

   $ dmpipe-pipeline --config config/master_dSphs.yaml


Extracting Analysis Results
---------------------------

The analysis pipeline produces a number of different outputs,
including:

* Combined plots of the DM interaction limits for the stacked analysis,
  as well as plots showing the results of the control tests.  These are
  in the dSphs/results directory.

* Indvidual plots of the SEDs and the DM interaction limits for each target
  in the analysis.

* Detailed intermediate results allowing the user to preproduce any plot
  or to refit any ROI or reproduce any other step of the analysis chain in isolation.


Loading and running interactively
---------------------------------

One can load the pipeline interactively in python, and see the current status of the analysis.

.. code-block:: python
   
   from dmpipe import Pipeline
   configfile = 'config/master_dSphs.yaml'
   pipe = Pipeline(linkname='dSphs')
   pipe.preconfigure(configfile)
   pipe.update_args(dict(config=configfile))

   # look at the current state
   pipe.print_status()
   
   # Continue running analysis starting from the previously saved
   # state 
   pipe.run()

Many other commands are demonstrated in the jupyter notebook example.

   

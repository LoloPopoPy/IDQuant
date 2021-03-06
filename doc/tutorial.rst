.. _tutorial:

Tutorial
========

Introduction
------------

Welcome to the IDQuant usage tutorial.

As input, IDQuant takes two different files:

    * A data file containing the MS experiment data that must be integrated
      using **Emzed C13 Profiler** software.
    * A calibration file that takes the form of a table.

As output, IDQuant will generate a directory where the notebook is situated using the
date & time as name. Different files will be put in this directory:

    * A pdf file containing the linear and regression curve plots and the associated
      residual plots for each metabolite
    * An excel file named *calibration_data* containing the clean data used for the
      calculation of the polynomials and the associated residual values
    * An excel file named *calculated_datas* containing the results of concentration
      predictions by **linear** and **quadratic** regression.
    * A text file with the run name + the extension *.log* in which the log of the
      run is recorded (equations, errors, warnings, etc...)
    * An excel file named *missing_values* containing the metabolites for which
      missing data was detected.

Modifying the calibration file
------------------------------

The concentrations used to create the calibration curve for each metabolite are
kept in the **calibration.csv** file.

To add metabolites and concentrations to the file, it is important to keep the structure
intact. If more calibration points need to be added, add their number at the top of the
table and the associated concentration value in the squares beneath for each metabolite.

.. warning:: The metabolite names given in the table must be exactly the same
             as the names used in the MS experiment.

Input data formatting
---------------------

The input data format is defined by Emzed C13 Profiler output. But IDQuant does not
use all the columns created by C13 Profiler, and so it is possible to create a custom
table to feed into the software. This custom table must have the minimum required columns
in a certain format (see example table below).

======== === ==== ======
compound mi  area source
======== === ==== ======
ATP       1   xx   expe
======== === ==== ======

Description of headers:
    * **Compound:** Names of the different metabolites to be analyzed
    * **mi:** Number of the isotopologue (see: :ref:`Isotopologue<Isotopologue>`).
      The M0 and the highest Mn must be referenced for the software to function
      correctly
    * **area:** Integration area
    * **source:** Sample name

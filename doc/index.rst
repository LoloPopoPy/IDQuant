.. IDQuant documentation master file, created by
   sphinx-quickstart on Fri Dec  4 10:45:10 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to IDQuant's documentation!
===================================

**IDQuant (Isotopic Dilution Quantification) is a software
that calculates sample concentrations of metabolites from C12
and C13 Mass Spectrometry Integrated Data**

It uses linear and quadratic regression equations to predict
concentrations using calibration data from (unlabeled and fully
labeled) Standard Molecules. It takes as input integrated C12/C13
MS data from **C13Profiler**.

It is one of the routine tools used by
`MetaToul platform <https://www6.toulouse.inrae.fr/metatoul>`_.

The code is open-source, and available on GitHub under a GPLv3 license.

.. rubric:: KEY FEATURES
	
	* **Calculation of metabolite concentrations using linear and quadratic regression equations **
	* **Automatic Regression and Residual plot building**
	* **Open-source, free and easy to install** (Python3 required)
	* **Simple, easy-to-use graphical interface built with ipywidgets and implemented through a Jupyter Notebook**


You can find information on the usage of IDQuant in the sections listed below.


.. toctree::
   :maxdepth: 2
   :caption: Usage Information
   
   quickstart.rst
   tutorial.rst

For more in-depth details on the inner workings of IDQuant, check the sections below.
  
.. toctree::
   :maxdepth: 2
   :caption: Miscellaneous
   
   definitions.rst
   librairy_doc.rst


For any requests or suggestions, please contact the author at legregam@insa-toulouse.fr.


* :ref:`search`



.. CyPDM documentation master file, created by
   sphinx-quickstart on Fri Dec 28 13:56:34 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

CyPDM
=================================
A fast package to apply the phase disperion minimization (PDM) algorithm, based on PyAstronomy module `PyPDM <https://github.com/sczesla/PyAstronomy>`_.

The package is written in C and wrapped with Cython, so that it can be called in both C and Python code.

The PDM alogrithm refers to the reference: `Stellingwerf (1978) <http://adsabs.harvard.edu/abs/1978ApJ...224..953S>`_.

Documentation
-------------

.. toctree::
   :maxdepth: 2
   
   getting_started.rst
   api.rst

Attribution
-----------

If you make use of this code, please cite it as:

.. code-block:: tex

  @misc{yan_rong_li_2018_2527266,
    author       = {Yan-Rong Li},
    title        = {{CyPDM: A fast package to apply the phase disperion 
                   minimization algorithm}},
    month        = dec,
    year         = 2018,
    doi          = {10.5281/zenodo.2527266},
    url          = {https://doi.org/10.5281/zenodo.2527266}
  }

Authors & License
-----------------

Copyright 2018 Yan-Rong Li.

Licensed under the `GNU General Public License v3.0 or later <https://opensource.org/licenses/GPL-3.0>`_.

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

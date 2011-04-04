nmrglue.sparky
==============

.. automodule:: nmrglue.analysis.analysisbase

This modules is imported as nmrglue.analysisbase and can be called as such.

Low-Level Functions
-------------------

These functions are typically not used directly by users.  They are called by
high level functions.

.. autofunction:: recnames
.. autofunction:: find_limits
.. autofunction:: limits2slice
.. autofunction:: slice2limits
.. autofunction:: squish
.. autofunction:: pick2linesh
.. autofunction:: linesh2pick
.. autofunction:: ls_str2class
.. autofunction:: center_fwhm 
.. autofunction:: center_fwhm_bymoments

Low-Level Classes
-----------------

.. autoclass:: gauss1D
.. autoclass:: peak1D 
.. autoclass:: lorentz1D 
.. autoclass:: scale1D
.. autoclass:: ndwindow 
.. autoclass:: ndwindow_index
.. autoclass:: ndwindow_inside
.. autoclass:: ndwindow_inside_index

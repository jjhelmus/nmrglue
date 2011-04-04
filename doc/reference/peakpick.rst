nmrglue.sparky
==============

.. automodule:: nmrglue.analysis.peakpick

This modules is imported as nmrglue.peakpick and can be called as such.

High-Level Functions
---------------------


These are the functions most users will use from the peakpick module.

.. autofunction:: pick 

Low-Level Functions
-------------------

These functions are typically not used directly by users.  They are called by
high level functions.  Developers and user who want fine control over peak
picking will be interested in these functions.

.. autofunction:: find_regions
.. autofunction:: regions2recarray
.. autofunction:: pad_regions
.. autofunction:: resolve_region_overlap
.. autofunction:: combine_regions
.. autofunction:: combine_2regions
.. autofunction:: limits_in_limits
.. autofunction:: is_overlappedND
.. autofunction:: is_overlapped1D
.. autofunction:: region2linesh
.. autofunction:: linesh2region
.. autofunction:: filter_by_distance
.. autofunction:: pts_in_limits
.. autofunction:: in_limits
.. autofunction:: guess_params_center
.. autofunction:: guess_params_segment
.. autofunction:: pick_connected
.. autofunction:: pick_downward
.. autofunction:: pick_thres
.. autofunction:: pick_thres_fast


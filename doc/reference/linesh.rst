nmrglue.linesh
==============

.. automodule:: nmrglue.analysis.linesh

This module is imported as nmrglue.linesh and can be called as such

High-Level Functions
--------------------

These are the functions the majority of users will use from the linesh module.

.. autofunction:: fit_spectrum
.. autofunction:: fit_NDregion
.. autofunction:: sim_NDregion
.. autofunction:: add_to_table
.. autofunction:: pack_table
.. autofunction:: unpack_table
.. autofunction:: estimate_scales


Low-Level Functions
-------------------

These functions are typically not used directly by users. They are called by
high level functions.  Developers and user who want fine control over lineshape
fitting may be interested in these functions.

.. autofunction:: f_NDregion
.. autofunction:: s_NDregion
.. autofunction:: s_single_NDregion
.. autofunction:: err_NDregion
.. autofunction:: make_slist
.. autofunction:: split_list


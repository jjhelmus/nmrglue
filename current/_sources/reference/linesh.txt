nmrglue.linesh
==============

.. automodule:: nmrglue.analysis.linesh

This module is imported as nmrglue.linesh and can be called as such

User Functions
--------------

These are the functions the majority of users will use from the linesh module.

.. autosummary:: 
    :toctree: generated/

    fit_spectrum
    fit_NDregion
    sim_NDregion
    
    add_to_table
    pack_table
    unpack_table
    estimate_scales


Developer Functions
-------------------

These functions are typically not used directly by users.  Developers who
want fine control over lineshape fitting may be interested in these 
functions.

.. autosummary::
    :toctree: generated/

    make_slist
    split_list
    calc_errors
    s_NDregion
    s_single_NDregion
    err_NDregion
    f_NDregion


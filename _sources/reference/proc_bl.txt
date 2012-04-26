.. _proc_bl:

nmrglue.proc_bl
===============

.. automodule:: nmrglue.process.proc_bl

This module is imported as nmrglue.proc_bl and can be called as such.

High-Level Functions
--------------------

.. autosummary::
    :toctree: generated/

    base
    cbf
    cbf_explicit
    med
    sol_general
    sol_boxcar
    sol_sine
    sol_sine2
    sol_gaussian

These functions are called by high-level function are and most uses will not
use them in common processing scripts.

Low-Level Functions
-------------------

.. autosummary::
    :toctree: generated/

    calc_bl_linear
    calc_bl_med

These functions are not implemented (yet) in nmrglue.  There are here as
shells for later implementation.

NotImplemented
--------------

.. autosummary::
    :toctree: generated/

    poly_td
    poly_fd

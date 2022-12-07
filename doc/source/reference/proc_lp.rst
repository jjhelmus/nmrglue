.. _proc_lp:

nmrglue.proc_lp
===============

.. automodule:: nmrglue.process.proc_lp

This module is imported as nmrglue.proc_lp and can be called as such.

User Functions
--------------

.. autosummary::
    :toctree: generated/

    lp
    lp_svd
    lp_qr
    lp_cho
    lp_tls
    lp2d
    cadzow
    lp_model



Developer Functions
-------------------

.. include:: ../../../nmrglue/process/proc_lp.py
    :start-line: 5
    :end-line: 20


These functions are called by high-level function are and most users will not
use them in common processing scripts.  Developers may be interested in them.

.. autosummary::
    :toctree: generated/

    lp_1d
    extrapolate_2d
    make_lp2d_Dd
    cadzow_single
    root2damp
    root2freq
    cof2amp
    cof2phase
    make_D
    make_little_d
    make_Dd
    make_mirror
    find_lpc
    find_lpc_svd
    pinv_diagsvd
    find_lpc_qr
    find_lpc_cholesky
    find_lpc_tls
    find_lpc_fb
    find_lpc_bf
    find_lproots_hsvd
    find_roots
    find_coeff
    reverse_filter
    fix_roots
    extrapolate

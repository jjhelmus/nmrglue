.. _proc_base:

nmrglue.proc_base
=================

.. automodule:: nmrglue.process.proc_base

This module is imported as nmrglue.proc_base and can be called as such.

Apodization
-----------

.. autosummary::
    :toctree: generated/

    em
    gm
    gmb
    jmod
    sp
    sine
    tm
    tri

Shifts
------

.. autosummary::
    :toctree: generated/

    rs
    ls
    cs
    roll
    fsh


Transforms
----------

.. autosummary::
    :toctree: generated/

    rft
    irft
    fft
    fft_norm
    fft_positive
    ifft
    ifft_norm
    ifft_positive
    ha
    ht

Standard NMR
------------

.. autosummary::
    :toctree: generated/

    di
    ps
    ps_exp
    tp
    ytp
    xy2yx
    tp_hyper
    zf_inter
    zf_pad
    zf
    zf_double
    zf_size
    zf_auto

Basic Utilities
---------------

.. autosummary::
    :toctree: generated/

    add
    add_ri
    dx
    ext
    ext_left
    ext_right
    ext_mid
    integ
    mc
    mc_pow
    mir_left
    mir_right
    mir_center
    mir_center_onepoint
    mult
    rev
    set
    set_complex
    set_real
    set_imag
    ri2c
    interleave_complex
    unpack_complex
    c2ri
    seperate_interleaved
    pack_complex
    decode_States
    ri2rr
    append_imag
    rr2ri
    unappend_imag
    exlr
    exchange_lr
    rolr
    rotate_lr
    swap
    swap_ri
    bswap
    byte_swap
    neg_left
    neg_right
    neg_middle
    neg_edges
    neg_all
    neg_real
    neg_imag
    neg_even
    neg_odd
    neg_alt
    abs
    sign

Misc
----

.. autosummary::
    :toctree: generated/

    coadd
    coad
    thres
    conv
    convolute
    corr
    correlate
    filter_median
    filter_min
    filter_max
    filter_percentile
    filter_rank
    filter_amin
    filter_amax
    filter_range
    filter_avg
    filter_dev
    filter_sum
    filter_generic
    nmr_reorder
    qart
    qart_auto
    gram_schmidt
    qmix
    smo
    center
    zd
    zd_boxcar
    zd_triangle
    zd_sinebell
    zd_gaussian

Low-Level Functions
-------------------
The following are functions called by other processing functions.  They
are included here for completeness.

.. autosummary::
    :toctree: generated/

    int2bin
    bin2int
    gray
    amin_flt
    amax_flt
    range_flt
    avg_flt
    std_flt
    sum_flt
    largest_power_of_2


Other processing functions
--------------------------

.. autosummary::
    :toctree: generated/

    expand_nus
nmrglue.rnmrtk
==============

.. automodule:: nmrglue.fileio.rnmrtk

This modules is imported as nmrglue.rnmrtk and can be called as such.

User Information
----------------

User Functions
^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    read
    write
    write_lowmem
    read_lowmem
    read_sec
    write_sec
    read_par
    write_par
    make_uc
    guess_udic
    create_dic


Developer Infomation
--------------------

Developer Functions
^^^^^^^^^^^^^^^^^^^

These functions are typically not used directly by users.  Developers who
want fine control over Rowland NMR Toolkit files will be interested in these
functions.

.. autosummary::
    :toctree: generated/

    create_data
    get_data
    get_trace
    put_trace
    uninterleave_data
    interleave_data
    find_shape
    permute_dic
    parse_par_line
    make_empty_dic

Developer Classes
^^^^^^^^^^^^^^^^^

.. autosummary:: 
    :toctree: generated/

    rnmrtk_nd 

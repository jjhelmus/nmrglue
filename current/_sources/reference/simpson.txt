.. _simpson_module:

nmrglue.simpson
===============

.. automodule:: nmrglue.fileio.simpson

This modules is imported as nmrglue.simpson and can be called as such.

User Information
----------------

User Functions
^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    read
    read_text
    read_xreim
    read_xyreim
    read_raw_bin_1d
    read_raw_bin_2d
    read_binary

Developer Infomation
--------------------

Developer Functions
^^^^^^^^^^^^^^^^^^^

These functions are typically not used directly by users.  Developers who
want fine control over SIMPSON files will be interested in these 
functions.

.. autosummary::
    :toctree: generated/

    guess_ftype
    unappend_data
    chars2bytes
    bytes2float
    bytes2float_bitarray

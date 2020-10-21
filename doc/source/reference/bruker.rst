nmrglue.bruker
==============

.. automodule:: nmrglue.fileio.bruker

This modules is imported as nmrglue.bruker and can be called as such.

User Information
----------------

User Functions
^^^^^^^^^^^^^^

These are functions which are targetted for users of nmrglue.

.. autosummary::
    :toctree: generated/

    read
    write
    read_pdata
    write_pdata
    remove_digital_filter
    read_lowmem
    write_lowmem
    read_binary
    write_binary
    read_pdata_binary
    scale_pdata
    read_binary_lowmem
    write_binary_lowmem
    read_jcamp
    write_jcamp
    read_pprog
    write_pprog
    guess_udic
    create_dic



Developer Infomation
--------------------

.. include:: ../../../nmrglue/fileio/bruker.py
    :start-line: 11
    :end-line: 29

Developer Functions
^^^^^^^^^^^^^^^^^^^

These functions are typically not used directly by users.  Developers who want
fine control over Bruker files may be interested in these functions.

.. autosummary::
    :toctree: generated/

    create_data
    add_axis_to_udic
    create_acqus_dic
    guess_shape
    guess_shape_and_submatrix_shape
    get_data
    put_data
    complexify_data
    uncomplexify_data
    reorder_submatrix
    rm_dig_filter
    parse_jcamp_line
    parse_jcamp_value
    write_jcamp_pair
    read_acqus_file
    read_procs_file


Developer Classes
^^^^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    bruker_nd

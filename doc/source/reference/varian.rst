.. _varian_module:

nmrglue.varian
==============

.. automodule:: nmrglue.fileio.varian

This modules is imported as nmrglue.varian and can be called as such.  These
functions and classes can also be access from nmrglue.agilent.

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
    read_fid
    write_fid
    read_fid_lowmem
    write_fid_lowmem 
    read_fid_ntraces
    read_procpar
    write_procpar
    guess_udic
    create_dic


Developer Infomation
--------------------

.. include:: ../../../nmrglue/fileio/varian.py
    :start-line: 8
    :end-line: 18

Developer Functions
^^^^^^^^^^^^^^^^^^^

These functions are typically not used directly by users.  Developers who
want fine control over Agilent/Varian files will be interested in these 
functions.

.. autosummary::
    :toctree: generated/

    create_data
    create_pdic_param
    find_torder
    torder2i2t
    torder2t2i
    reorder_data
    order_data
    get_nblocks
    get_block
    get_nblocks_ntraces
    get_block_ntraces
    get_trace
    get_fileheader
    get_blockheader
    skip_blockheader
    get_hyperheader
    put_block
    put_trace
    put_fileheader
    put_blockheader
    put_hyperheader
    hyperheader2dic
    repack_hyperheader
    dic2hyperheader
    make_blockheader
    blockheader2dic
    repack_blockheader
    dic2blockheader
    fileheader2dic
    repack_fileheader
    dic2fileheader
    find_shape
    find_cdtype
    find_dtype
    uninterleave_data 
    interleave_data
    get_parameter


Developer Classes
^^^^^^^^^^^^^^^^^

.. autosummary:: 
    :toctree: generated/

    fid_nd 

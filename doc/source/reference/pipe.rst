nmrglue.pipe
==============

.. automodule:: nmrglue.fileio.pipe

This modules is imported as nmrglue.pipe and can be called as such.

User Information
----------------

User Functions
^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    read
    write 
    read_lowmem
    write_lowmem
    read_table
    write_table
    make_uc
    guess_udic
    create_dic
    datetime2dic
    dic2datetime

User Classes
^^^^^^^^^^^^

.. _pipe_iter3D:

.. autosummary::
    :toctree: generated/

    iter3D

Developer Information
---------------------

.. include:: ../../../nmrglue/fileio/pipe.py
    :start-line: 7
    :end-line: 8

Developer Functions
^^^^^^^^^^^^^^^^^^^

These functions are typically not used directly by users.  Developers who want
fine control over NMRPipe files may be interested in these functions.

.. autosummary::
    :toctree: generated/

    create_data
    add_axis_to_dic
    create_empty_dic
    read_1D
    read_2D 
    read_lowmem_2D
    read_stream
    read_lowmem_stream
    read_3D
    read_lowmem_3D
    read_4D
    read_lowmem_4D
    write_single
    write_3D
    write_4D
    write_lowmem_2D
    write_lowmem_3D
    write_lowmem_3Ds 
    write_lowmem_4D
    write_lowmem_4Ds
    put_fdata
    put_trace 
    put_data
    write_slice_3D
    pack_complex
    transpose_3D 
    find_shape
    reshape_data
    unshape_data 
    unappend_data
    append_data
    fdata2dic
    dic2fdata 
    get_fdata
    get_data
    get_fdata_data
    get_trace


Developer Classes
^^^^^^^^^^^^^^^^^
.. autosummary::
    :toctree: generated/

    pipe_2d
    pipe_3d
    pipestream_3d
    pipe_4d
    pipestream_4d

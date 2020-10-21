nmrglue.sparky
==============

.. automodule:: nmrglue.fileio.sparky

This modules is imported as nmrglue.sparky and can be called as such.

User Information
----------------

User Functions
^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/

    read
    write
    read_savefile
    read_lowmem
    write_lowmem
    make_uc
    guess_udic
    create_dic
    datetime2dic 
    dic2datetime

Developer Information
---------------------

.. include:: ../../../nmrglue/fileio/sparky.py
    :start-line: 7
    :end-line: 13

Developer Functions
^^^^^^^^^^^^^^^^^^^

These functions are typically not used directly by users.  Developers who 
want fine control over Sparky files may be interested in these functions.

.. autosummary::
    :toctree: generated/

    create_data
    create_axisdic
    calc_tshape
    read_2D
    write_2D
    read_3D
    write_3D
    read_lowmem_2D
    read_lowmem_3D
    get_tilen
    get_tile
    put_tile
    get_data
    put_data
    find_tilen_2d
    tile_data2d
    untile_data2D
    find_tilen_3d
    tile_data3d
    untile_data3D
    get_fileheader 
    put_fileheader 
    fileheader2dic
    dic2fileheader
    get_axisheader
    put_axisheader
    axisheader2dic 
    dic2axisheader

Developer Classes
^^^^^^^^^^^^^^^^^

.. autosummary::
    :toctree: generated/
    
    sparky_2d
    sparky_3d

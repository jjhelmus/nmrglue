nmrglue.sparky
==============

.. automodule:: nmrglue.fileio.sparky

This modules is imported as nmrglue.sparky and can be called as such.

High-Level Functions
---------------------


These are the functions most users will use from the sparky module.

.. autofunction:: read
.. autofunction:: read_lowmem
.. autofunction:: write
.. autofunction:: write_lowmem
.. autofunction:: make_uc
.. autofunction:: guess_udic
.. autofunction:: create_dic
.. autofunction:: datetime2dic 
.. autofunction:: dic2datetime

Low-Level Functions
-------------------

These functions are typically not used directly by users.  They are called by
high level functions.  Developers and user who want fine control over sparky
files will be interested in these functions.

.. autofunction:: create_data
.. autofunction:: create_axisdic
.. autofunction:: calc_tshape
.. autofunction:: read_2D
.. autofunction:: write_2D
.. autofunction:: read_3D
.. autofunction:: write_3D
.. autofunction:: read_lowmem_2D
.. autofunction:: read_lowmem_3D
.. autofunction:: get_tilen
.. autofunction:: get_tile
.. autofunction:: put_tile
.. autofunction:: get_data
.. autofunction:: put_data
.. autofunction:: find_tilen_2d
.. autofunction:: tile_data2d
.. autofunction:: untile_data2D
.. autofunction:: find_tilen_3d
.. autofunction:: tile_data3d
.. autofunction:: untile_data3D
.. autofunction:: get_fileheader 
.. autofunction:: put_fileheader 
.. autofunction:: fileheader2dic
.. autofunction:: dic2fileheader
.. autofunction:: get_axisheader
.. autofunction:: put_axisheader
.. autofunction:: axisheader2dic 
.. autofunction:: dic2axisheader

Low-Level Classes
-----------------

.. autoclass:: sparky_2d
.. autoclass:: sparky_3d

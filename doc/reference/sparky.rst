nmrglue.sparky
==============

.. automodule:: nmrglue.fileio.sparky

This modules is imported as nmrglue.sparky and can be called as such.

High-Level Functions
---------------------


These are the functions the majority of users will use from the sparky module.

.. autofunction:: read
.. autofunction:: read_2D
.. autofunction:: read_3D
.. autofunction:: read_lowmem
.. autofunction:: read_lowmem_2D
.. autofunction:: read_lowmem_3D
.. autofunction:: write
.. autofunction:: write_2D
.. autofunction:: write_3D
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
.. autofunction:: get_tilen
.. autofunction:: get_tile
.. autofunction:: get_data
.. autofunction:: get_fileheader
.. autofunction:: get_axisheader
.. autofunction:: put_tile
.. autofunction:: put_data
.. autofunction:: put_fileheader
.. autofunction:: put_axisheader
.. autofunction:: find_tilen_2d
.. autofunction:: find_tilen_3d
.. autofunction:: tile_data2d
.. autofunction:: tile_data3d
.. autofunction:: untile_data2D
.. autofunction:: untile_data3D
.. autofunction:: calc_tshape
.. autofunction:: fileheader2dic
.. autofunction:: dic2fileheader
.. autofunction:: axisheader2dic
.. autofunction:: dic2axisheader

Low-Level Classes
-----------------

.. autoclass:: sparky_2d
.. autoclass:: sparky_3d

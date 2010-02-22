nmrglue.varian
==============

.. automodule:: nmrglue.fileio.varian

This modules is imported as nmrglue.varian and can be called as such.

High-Level Functions
---------------------

These are the functions the majority of users will use from the varian module.

.. autofunction:: read_fid
.. autofunction:: read_fid_lowmem
.. autofunction:: read_fid_lowmem_2D
.. autofunction:: read_fid_lowmem_3D
.. autofunction:: write_fid
.. autofunction:: write_fid_1D
.. autofunction:: write_fid_2D
.. autofunction:: write_fid_3D
.. autofunction:: read_procpar
.. autofunction:: write_procpar
.. autofunction:: guess_udic
.. autofunction:: create_dic
.. autofunction:: sign_adj_2Ddata


Low-Level Functions
-------------------


These functions are typically not used directly by users.  They are called by
high level functions.  Developers and user who want fine control over varian
files will be interested in these functions.

.. autofunction:: create_data
.. autofunction:: get_parameter
.. autofunction:: find_cdtype
.. autofunction:: find_dtype
.. autofunction:: uninterleave_data
.. autofunction:: interleave_data
.. autofunction:: get_block
.. autofunction:: get_nblocks
.. autofunction:: get_fileheader
.. autofunction:: get_blockheader
.. autofunction:: get_hyperheader
.. autofunction:: get_trace
.. autofunction:: skip_blockheader
.. autofunction:: put_block
.. autofunction:: put_nblocks
.. autofunction:: put_fileheader
.. autofunction:: put_blockheader
.. autofunction:: put_hyperheader
.. autofunction:: put_trace
.. autofunction:: hyperheader2dic
.. autofunction:: blockheader2dic
.. autofunction:: fileheader2dic
.. autofunction:: dic2hyperheader
.. autofunction:: dic2blockheader
.. autofunction:: dic2fileheader
.. autofunction:: repack_hyperheader
.. autofunction:: repack_blockheader
.. autofunction:: repack_fileheader
.. autofunction:: make_blockheader


Low-Level Classes
-----------------


.. autoclass:: fid_2d
.. autoclass:: fid_3d

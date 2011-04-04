nmrglue.varian
==============

.. automodule:: nmrglue.fileio.varian

This modules is imported as nmrglue.varian and can be called as such.

High-Level Functions
---------------------

These are the functions most users will use from the varian module.

.. autofunction:: read
.. autofunction:: read_lowmem
.. autofunction:: write
.. autofunction:: write_lowmem
.. autofunction:: read_fid
.. autofunction:: read_fid_lowmem
.. autofunction:: read_fid_ntraces
.. autofunction:: write_fid
.. autofunction:: write_fid_lowmem 
.. autofunction:: read_procpar
.. autofunction:: write_procpar
.. autofunction:: guess_udic
.. autofunction:: create_dic


Low-Level Functions
-------------------


These functions are typically not used directly by users.  They are called by
high level functions.  Developers and user who want fine control over varian
files will be interested in these functions.

.. autofunction:: create_data
.. autofunction:: create_pdic_param
.. autofunction:: find_torder
.. autofunction:: torder2i2t
.. autofunction:: torder2t2i
.. autofunction:: reorder_data
.. autofunction:: order_data
.. autofunction:: get_nblocks
.. autofunction:: get_block
.. autofunction:: get_nblocks_ntraces
.. autofunction:: get_block_ntraces
.. autofunction:: get_trace
.. autofunction:: get_fileheader
.. autofunction:: get_blockheader
.. autofunction:: skip_blockheader
.. autofunction:: get_hyperheader
.. autofunction:: put_block
.. autofunction:: put_trace
.. autofunction:: put_fileheader
.. autofunction:: put_blockheader
.. autofunction:: put_hyperheader
.. autofunction:: hyperheader2dic
.. autofunction:: repack_hyperheader
.. autofunction:: dic2hyperheader
.. autofunction:: make_blockheader
.. autofunction:: blockheader2dic
.. autofunction:: repack_blockheader
.. autofunction:: dic2blockheader
.. autofunction:: fileheader2dic
.. autofunction:: repack_fileheader
.. autofunction:: dic2fileheader
.. autofunction:: find_shape
.. autofunction:: find_cdtype
.. autofunction:: find_dtype
.. autofunction:: uninterleave_data 
.. autofunction:: interleave_data
.. autofunction:: get_parameter


Low-Level Classes
-----------------


.. autoclass:: fid_nd 

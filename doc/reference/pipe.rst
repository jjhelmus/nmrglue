nmrglue.pipe
==============

.. automodule:: nmrglue.fileio.pipe

This modules is imported as nmrglue.pipe and can be called as such.

High-Level Functions
---------------------

These are the functions most users will use from the pipe module.

.. autofunction:: read
.. autofunction:: read_lowmem
.. autofunction:: write 
.. autofunction:: write_lowmem
.. autofunction:: read_table
.. autofunction:: write_table
.. autofunction:: make_uc
.. autofunction:: guess_udic
.. autofunction:: create_dic
.. autofunction:: datetime2dic
.. autofunction:: dic2datetime



High-Level Classes
-------------------


These are the classes most users will use from the pipe module.

.. _pipe_iter3D:

.. autoclass:: iter3D


Low-Level Functions
-------------------


These functions are typically not used directly by users.  They are called by
high level functions.  Developers and user who want fine control over pipe
files will be interested in these functions.

.. autofunction:: create_data
.. autofunction:: add_axis_to_dic
.. autofunction:: create_empty_dic
.. autofunction:: read_1D
.. autofunction:: read_2D 
.. autofunction:: read_lowmem_2D
.. autofunction:: read_stream
.. autofunction:: read_lowmem_stream
.. autofunction:: read_3D
.. autofunction:: read_lowmem_3D
.. autofunction:: read_4D
.. autofunction:: read_lowmem_4D
.. autofunction:: write_single
.. autofunction:: write_3D
.. autofunction:: write_4D
.. autofunction:: write_lowmem_2D
.. autofunction:: write_lowmem_3D
.. autofunction:: write_lowmem_3Ds 
.. autofunction:: write_lowmem_4D
.. autofunction:: write_lowmem_4Ds
.. autofunction:: put_fdata
.. autofunction:: put_trace 
.. autofunction:: put_data
.. autofunction:: write_slice_3D
.. autofunction:: pack_complex
.. autofunction:: transpose_3D 
.. autofunction:: find_shape
.. autofunction:: reshape_data
.. autofunction:: unshape_data 
.. autofunction:: unappend_data
.. autofunction:: append_data
.. autofunction:: fdata2dic
.. autofunction:: dic2fdata 
.. autofunction:: get_fdata
.. autofunction:: get_data
.. autofunction:: get_fdata_data
.. autofunction:: get_trace


Low-Level Classes
-----------------

.. autoclass:: pipe_2d
.. autoclass:: pipe_3d
.. autoclass:: pipestream_3d
.. autoclass:: pipe_4d
.. autoclass:: pipestream_4d

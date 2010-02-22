nmrglue.pipe
==============

.. automodule:: nmrglue.fileio.pipe

This modules is imported as nmrglue.pipe and can be called as such.

High-Level Functions
---------------------


These are the functions the majority of users will use from the pipe module.

.. autofunction:: read
.. autofunction:: read_1D
.. autofunction:: read_2D
.. autofunction:: read_3D
.. autofunction:: read_lowmem
.. autofunction:: read_lowmem_3D
.. autofunction:: write
.. autofunction:: write_1D
.. autofunction:: write_2D
.. autofunction:: write_3D
.. autofunction:: read_tab
.. autofunction:: write_tab
.. autofunction:: make_uc
.. autofunction:: guess_udic
.. autofunction:: create_dic
.. autofunction:: datetime2dic
.. autofunction:: dic2datetime



High-Level Classes
-------------------


These are the classes the majority of users will use from the pipe module.

.. autoclass:: iter3D


Low-Level Functions
-------------------


These functions are typically not used directly by users.  They are called by
high level functions.  Developers and user who want fine control over pipe
files will be interested in these functions.

.. autofunction:: create_data
.. autofunction:: create_empty_dic
.. autofunction:: add_axis_to_dic
.. autofunction:: flist_from_filemask
.. autofunction:: write_slice_3D
.. autofunction:: transpose_3D
.. autofunction:: find_shape
.. autofunction:: reshape_data
.. autofunction:: unshape_data
.. autofunction:: pack_complex
.. autofunction:: unappend_data
.. autofunction:: append_data
.. autofunction:: fdata2dic
.. autofunction:: dic2fdata
.. autofunction:: get_fdata
.. autofunction:: get_data
.. autofunction:: get_fdata_data
.. autofunction:: put_data


Low-Level Classes
-----------------


.. autoclass:: pipe_3d

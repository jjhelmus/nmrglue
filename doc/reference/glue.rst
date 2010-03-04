nmrglue.glue
==============

.. automodule:: nmrglue.fileio.glue

This modules is imported as nmrglue.glue and can be called as such.

High-Level Functions
---------------------


These are the functions the majority of users will use from the glue module.

.. autofunction:: read
.. autofunction:: read_lowmem
.. autofunction:: write
.. autofunction:: make_uc
.. autofunction:: guess_udic
.. autofunction:: create_dic


Low-Level Functions
-------------------


These functions are typically not used directly by users.  They are called by
high level functions.  Developers and user who want fine control over glue
files will be interested in these functions.

.. autofunction:: create_data
.. autofunction:: wrap_data
.. autofunction:: get_dic
.. autofunction:: put_dic


Low-Level Classes
-----------------

.. autoclass:: glue_2d
.. autoclass:: glue_3d

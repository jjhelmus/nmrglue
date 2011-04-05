nmrglue.bruker
==============

.. automodule:: nmrglue.fileio.bruker

This modules is imported as nmrglue.bruker and can be called as such.

High-Level Functions
---------------------

These are the functions that most users will use from the bruker module.

.. autofunction:: read
.. autofunction:: read_lowmem
.. autofunction:: write
.. autofunction:: write_lowmem
.. autofunction:: read_binary
.. autofunction:: read_binary_lowmem
.. autofunction:: write_binary
.. autofunction:: write_binary_lowmem
.. autofunction:: read_jcamp
.. autofunction:: write_jcamp
.. autofunction:: read_pprog
.. autofunction:: write_pprog
.. autofunction:: remove_digital_filter
.. autofunction:: guess_udic
.. autofunction:: create_dic


Low-Level Functions
-------------------


These functions are typically not used directly by users.  They are called by
high level functions.  Developers and user who want fine control over bruker
files will be interested in these functions.

.. autofunction:: create_data
.. autofunction:: create_acqus_dic
.. autofunction:: guess_shape
.. autofunction:: get_data
.. autofunction:: put_data
.. autofunction:: complexify_data
.. autofunction:: uncomplexify_data
.. autofunction:: rm_dig_filter
.. autofunction:: parse_jcamp_line
.. autofunction:: write_jcamp_pair


Low-Level Classes
-----------------

.. autoclass:: bruker_nd

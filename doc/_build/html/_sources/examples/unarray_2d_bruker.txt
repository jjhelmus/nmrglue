.. _unarray_2d_bruker:

unarray example: unarray_2d_bruker 
==================================

This example shows how to use nmrglue to seperate bruker data collected with
an innermost parameter interleaved.  The full experimental data in the 
``arrayed_data.dir`` directory is unpacked into a series of directories named 
``1`` , ``2`` , ``3`` , ... which can be converted with nmrglue or NMRPipe. 
The data shape, array size and additional files to copy to the new directories 
must be determined by the user.

[`source code <el/unarray/unarray_2d_bruker/unarray.py>`_]

.. literalinclude:: el/unarray/unarray_2d_bruker/unarray.py

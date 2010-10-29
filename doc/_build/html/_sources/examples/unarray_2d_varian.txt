.. _unarray_2d_varian:

unarray example: unarray_2d_varian 
==================================

This example shows how to use nmrglue to seperate varian data collected with
an innermost parameter interleaved.  The full experimental data in the 
``arrayed_data.dir`` directory is unpacked into a series of directories with 
names ``tXmix_*.fid`` which can be converted with nmrglue or NMRPipe. The name
and values of the interleaved parameter is determined from the ``procpar`` 
file in the ``arrayed_data.dir`` directory.

[`source code <el/unarray/unarray_2d_varian/unarray.py>`_]

.. literalinclude:: el/unarray/unarray_2d_varian/unarray.py

.. _unarray_1d_varian:

unarray example: unarray_1d_varian 
==================================

This example shows how to use nmrglue to seperate varian data collected with
an innermost parameter interleaved.  The full experimental data in the 
``arrayed_data.dir`` directory is unpacked into a series of directories with 
names ``tHX_*.fid`` which can be converted with nmrglue or NMRPipe. The name
and values of the interleaved parameter is determined from the ``procpar`` 
file in the ``arrayed_data.dir`` directory.

[:download:`source code <../../../examples/unarray/unarray_1d_varian/unarray.py>`]

.. literalinclude:: ../../../examples/unarray/unarray_1d_varian/unarray.py

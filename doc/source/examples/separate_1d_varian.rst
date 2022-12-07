.. _separate_1d_varian:

separate example: separate_1d_varian
====================================

This example shows how to use nmrglue to separate Agilent/Varian data
collected with an innermost parameter interleaved.  The full experimental
data in the ``arrayed_data.dir`` directory is unpacked into a series of
directories with names ``tHX_*.fid`` which can be converted with nmrglue or
NMRPipe. The name and values of the interleaved parameter is determined from
the ``procpar`` file in the ``arrayed_data.dir`` directory.

The data used in this example is available for
`download. <https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/nmrglue/example_separate_1d_varian.zip>`_

[:download:`source code <../../../examples/separate/separate_1d_varian/separate.py>`]

.. literalinclude:: ../../../examples/separate/separate_1d_varian/separate.py

.. _separate_2d_bruker:

separate example: separate_2d_bruker 
====================================

This example shows how to use nmrglue to separate Bruker data collected with
an innermost parameter interleaved.  The full experimental data in the 
``arrayed_data.dir`` directory is unpacked into a series of directories named 
``1`` , ``2`` , ``3`` , ... which can be converted with nmrglue or NMRPipe. 
The data shape, array size and additional files to copy to the new directories 
must be determined by the user.

The data used in this example is available for 
`download. <https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/nmrglue/example_separate_2d_bruker.zip>`_

[:download:`source code <../../../examples/separate/separate_2d_bruker/separate.py>`]

.. literalinclude:: ../../../examples/separate/separate_2d_bruker/separate.py

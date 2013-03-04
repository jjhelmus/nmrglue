.. _process_pipe_3d:

process example: process_pipe_3d
================================

This example shows how nmrglue can be used to process NMR data.  In this 
script a 3D time domain NMRPipe file is processing into a 3D NMRPipe 
frequency domain file. For 3D processing the iter3D object is used to loop 
over XY and ZX planes.  Detail on this object can be found in the
:ref:`varian_module` documentation.

The data used in this example is available for 
`download. <http://nmrglue.googlecode.com/files/example_process_pipe_3d.zip>`_

[:download:`source code <../../../examples/processing/process_pipe_3d.py>`]

.. literalinclude:: ../../../examples/processing/process_pipe_3d.py

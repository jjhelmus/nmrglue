.. _process_pipe_3d:

process example: process_pipe_3d
================================

This example shows how nmrglue can be used to process NMR data.  In this 
script a 3D time domain NMRPipe file is processing into a 3D NMRPipe 
frequency domain file. For 3D processing the XXX iter3D object is used to 
loop over XY and ZX planes.After processing the resulting file is compared to
one processed with NMRPipe.


[`source code <el/processing/process_pipe_3d.py>`_]

.. literalinclude:: el/processing/process_pipe_3d.py

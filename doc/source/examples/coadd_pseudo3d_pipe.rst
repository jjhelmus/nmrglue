.. _coadd_pseudo3d_pipe:

coadd example: coadd_pseudo3d_pipe
==================================

This example shows how to use nmrglue to coadd a number of collections of 2D
NMRPipe files which constitute a pseudo-3D data set.  The two (or more) 
psuedo-3D data sets are assumed to be in directories named  ``run*.dir`` with 
subdirectories named ``*.fid`` containing a ``test.fid`` file.  The directory 
``coadded_data.dir`` is created with the same subdirectory structure containing 
``test.fid`` files containing data created by coadding each pseudo-3D.

[:download:`source code <../../../examples/coadd/coadd_pseudo3d_pipe/coadd_pseudo3d.py>`]

.. literalinclude:: ../../../examples/coadd/coadd_pseudo3d_pipe/coadd_pseudo3d.py

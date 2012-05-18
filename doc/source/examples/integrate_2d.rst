.. _integrate_2d:

integration example: integrate_2d
==================================

This example shows how to use nmrglue to integrate a 2D NMRPipe spectra.  The
script reads in point limits from ``limits.in`` and takes a simple 
summation integral of all points in each box described.  The integrated 
volumes are writting to ``volumes.out``.  For a method to graphically examine
these limits see :ref:`plot_2d_boxes`.  Similarly to check the peak 
assignments see :ref:`plot_2d_assignments`.


[:download:`source code <../../../examples/integration/integrate_2d/integrate_2d.py>`]

.. literalinclude:: ../../../examples/integration/integrate_2d/integrate_2d.py


[:download:`input file <../../../examples/integration/integrate_2d/limits.in>`]

.. literalinclude:: ../../../examples/integration/integrate_2d/limits.in

Results:

[:download:`output file <../../../examples/integration/integrate_2d/volumes.out>`]

.. literalinclude:: ../../../examples/integration/integrate_2d/volumes.out

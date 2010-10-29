.. _integrate_1d:

integration example: integrate_1d
==================================

This example shows how to use nmrglue to integrate a 1D NMRPipe spectra.  The
script reads in ppm peak limits from ``limits.in`` and takes a simple 
summation integral of each peak using the spectra contained in ``1d_data.ft``.  The integration values are writting to ``area.out`` and a plot is make showing
the integration limits and values overlayed on the spectra to ``plot.png``.


[`source code <el/integration/integrate_1d/integrate_1d.py>`_]

.. literalinclude:: el/integration/integrate_1d/integrate_1d.py


[`input file <el/integration/integrate_1d/limits.in>`_]

.. literalinclude:: el/integration/integrate_1d/limits.in

Results:

[`output file <el/integration/integrate_1d/area.out>`_]

.. literalinclude:: el/integration/integrate_1d/area.out

[`figure <el/integration/integrate_1d/plot.png>`_]

.. image:: el/integration/integrate_1d/plot.png

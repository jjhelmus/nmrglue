.. _integrate_1d:

integration example: integrate_1d
==================================

This example shows how to use nmrglue to integrate a 1D NMRPipe spectra.  The
script reads in ppm peak limits from ``limits.in`` and takes a simple
summation integral of each peak using the spectra contained in ``1d_data.ft``.  The integration values are writing to ``area.out`` and a plot is make showing
the integration limits and values overlaid on the spectra to ``plot.png``.

The data used in this example is available for
`download. <https://storage.googleapis.com/google-code-archive-downloads/v2/code.google.com/nmrglue/example_integrate_1d.zip>`_

[:download:`source code <../../../examples/integration/integrate_1d/integrate_1d.py>`]

.. literalinclude:: ../../../examples/integration/integrate_1d/integrate_1d.py


[:download:`input file <../../../examples/integration/integrate_1d/limits.in>`]

.. literalinclude:: ../../../examples/integration/integrate_1d/limits.in

Results:

[:download:`output file <../../../examples/integration/integrate_1d/area.out>`]

.. literalinclude:: ../../../examples/integration/integrate_1d/area.out

[:download:`figure <../../../examples/integration/integrate_1d/plot.png>`]

.. image:: ../../../examples/integration/integrate_1d/plot.png

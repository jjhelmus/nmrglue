.. _apod_viewer_1win:

application example: apod_viewer_1win
=====================================

This is a sample GUI application showing how nmrglue can be used with
additional python modules like
`matplotlib <http://matplotlib.sourceforge.net/index.html>`_ and
`wxPython <http://www.wxpython.org/>`_ to create full fledged NMR applications.

In this application users can examine graphically the apodization windows
produced by the various window functions supported by NMRPipe.  In this example
the canvas in which the apodization windows are drawn and the location to input
the apodization parameter are contained in the same window.  The
:ref:`apod_viewer_2win` example has the canvas and input area in separate
windows.

[:download:`source code <../../../examples/sample_applications/apod_viewer_1win.py>`]

.. literalinclude:: ../../../examples/sample_applications/apod_viewer_1win.py

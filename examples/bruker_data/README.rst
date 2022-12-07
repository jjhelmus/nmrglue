.. _proc_bruker_1d:

Processing 1D Bruker Data
=========================

This example shows how nmrglue can be used to process and display one
dimensional Bruker data.

Raw Bruker data from modern spectrometers contains a *group delay artifact*
which must be removed during processing.  There has been much speculation as to
the origins of this artifact and many methods for removing the artifact have
been suggested [`1`_], [`2`_], [`3`_], [`4`_], [`5`_].

Nmrglue provides an algorithm for removing this artifact based on the protocol
presented in "DMX DIGITAL FILTERS AND NON-BRUKER OFFLINE PROCESSING III" by
W. M. Westler and F.  Abildgaard.
This method is available for use through
:py:func:`nmrglue.fileio.bruker.remove_digital_filter`.
Nmrglue users can use this included function to remove the artifact or
implement their own method if they are unsatisfied with the results.

.. _1: http://nmr-analysis.blogspot.com/2008/02/why-arent-bruker-fids-time-corrected.html
.. _2: http://nmr-analysis.blogspot.com/2010/05/bruker-smiles.html
.. _3: http://www.ebyte.it/stan/blog08a.html
.. _4: http://qa.nmrwiki.org/question/196/bruker-grpdly-parameter
.. _5: http://www.acornnmr.com/NutsHelp/rd_show.html

In this example a 1D NMR spectrum of 1,3 diaminopropane is processed and
plotted using nmrglue.  The results can be compared with the spectrum produced
from NMRPipe which provides a different artifact removal algorithm.  Note that
no apodization or baseline corrections are performed on these spectra.

Instructions
------------

* Download the 1D proton spectrum of 1,3 diaminopropane and unpack in this
  directory. This raw data is available from the
  `Madison Metabolomics Consortium Database`_ as `expnmr_00001_1.tar`_.

.. _Madison Metabolomics Consortium Database: http://mmcd.nmrfam.wisc.edu/
.. _expnmr_00001_1.tar: http://mmcd.nmrfam.wisc.edu/rawnmr/expnmr_00001_1.tar

* Execute ``process_and_plot_nmrglue.py`` to process and plot the 1D spectrum.
  This creates the file ``figure_nmrglue.png``.

* Optionally, the data can be processed with NMRPipe using the script
  ``nmrpipe_proc.com``.  Then ``plot_nmrpipe.py`` can be used to plot the
  resulting spectrum.  This creates the file ``figure_nmrpipe.png``.

``process_and_plot_nmrglue.py`` [:download:`source code <../../../examples/bruker_data/process_and_plot_nmrglue.py>`]

.. literalinclude:: ../../../examples/bruker_data/process_and_plot_nmrglue.py

Output:

[:download:`figure_nmrglue.png <../../../examples/bruker_data/figure_nmrglue.png>`]

.. image:: ../../../examples/bruker_data/figure_nmrglue.png

``nmrpipe_proc.com`` [:download:`source code <../../../examples/bruker_data/nmrpipe_proc.com>`]

.. literalinclude:: ../../../examples/bruker_data/nmrpipe_proc.com


``plot_nmrpipe.py`` [:download:`source code <../../../examples/bruker_data/plot_nmrpipe.py>`]

.. literalinclude:: ../../../examples/bruker_data/plot_nmrpipe.py

Output:

[:download:`figure_nmrpipe.png <../../../examples/bruker_data/figure_nmrpipe.png>`]

.. image:: ../../../examples/bruker_data/figure_nmrpipe.png

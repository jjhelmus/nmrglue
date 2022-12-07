.. _simulate_ucsf:

Simulate a ucsf (Sparky) file
=============================


This example shows how to use the
:py:func:`nmrglue.analysis.linesh.sim_NDregion` function to simulate a HSQC
spectrum with peak locations provided by a text file.  The simulated spectrum
is saved as a ucsf (Sparky) file named `test.ucsf`.

[:download:`make_ucsf.py <../../../examples/simulation/simulate_ucsf/make_ucsf.py>`]

.. literalinclude:: ../../../examples/simulation/simulate_ucsf/make_ucsf.py

[:download:`peaks.txt <../../../examples/simulation/simulate_ucsf/peaks.txt>`]

.. literalinclude:: ../../../examples/simulation/simulate_ucsf/peaks.txt

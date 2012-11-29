.. _testing:

=======
Testing
=======

Tests for verifying the functionality of nmrglue are available in the test
directory.  These tests use the nose_ testing infrastructure.

.. _nose: https://nose.readthedocs.org/en/latest/


Requirements
------------
To run these tests NumPy, SciPy, nmrglue, and nose must be installed and in the
Python search path.  NMRPipe must be installed to run the pipe_proc tests.

In addition, the location of the the test data sets must be specified in the 
``setup.py`` file in the test directory.  The `nmrglue test data`_ is available for download. 

.. _`nmrglue test data`: http://code.google.com/p/nmrglue/downloads/list?q=label:Test-Data

In order to run all nmrglue unit tests, the tests data sets must be 
downloaded, unpacked, and the all conversions scripts contained in the
archive must be run.  Many of these scripts require additional NMR software 
(NMRPipe, etc), see the ``README`` file in the test data achive for additional 
details.  A subset of the full test suite can be run without installing any 
additional software.


Running the unit tests
----------------------
After ensuring that all required packages are installed and ``setup.py`` 
correctly points to the location of the test data directory, the unit tests can
be run using the following::

    nosetest

Unit tests for a specific module can be run using::

    nosetest tests/test_pipe.py

Additional information on the `usage of the nosetest command
<https://nose.readthedocs.org/en/latest/usage.html>`_ is available.

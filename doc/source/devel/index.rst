.. _development-guide:

==================
Developement Guide
==================

This guide provides instructions for setting up an environment for developing
nmrglue and an overview of the project layout and contribution process.


Requirements
------------

To create an environment for developing nmrglue the following must be installed
and available.

* `Numpy <http://numpy.scipy.org>`_

* `SciPy <http://scipy.org>`_

* `nose <https://nose.readthedocs.org/en/latest/>`_

* `Sphinx <http://sphinx-doc.org/>`_

* `git <http://git-scm.com>`_

In addition the following Python packages are highly recommended. 

* `matplotlib <http://matplotlib.org/>`_

* `IPython <http://ipython.org/>`_

A easy way of obtaining and installing these packages is to use a Python 
distribution which provides these packages, such as 
`EPD <http://www.enthought.com/products/epd.php>`_.  Detailed information on
`installing a Scipy stack <http://scipy.github.com/install.html>`_ is 
available.

Finally, other NMR software packages must be installed to process and convert 
the test and example data.  These are not required for using nmrglue, but are
needed to verify its functionality and to run some of the examples.

* `NMRPipe <http://spin.niddk.nih.gov/NMRPipe/>`_

* `Sparky <http://www.cgl.ucsf.edu/home/sparky/>`_

* `The Rowland NMR Toolkit <http://rnmrtk.uchc.edu/rnmrtk/RNMRTK.html>`_

* `SIMPSON <http://bionmr.chem.au.dk/bionmr/software/simpson.php>`_

.. _source_code:

Source Code
-----------

nmrglue uses `github <http://github.com>`_ for source code hosting.  For access
to the source code, see the 
`nmrglue github <http://github.com/jjhelmus/nmrglue>`_ site.


To check out the latest version of nmrglue use git::
    
    git clone git://github.com/jjhelmus/nmrglue.git

nmrglue is a pure python module, the root directory can be included in your
PYTHONPATH directly, or a symbolic link can be added to the *site-packages*
directory of your Python install.  In this way any modifications to the nmrglue
source tree will be picked up when nmrglue is imported.  


Test and Example Data 
---------------------

nmrglue uses experimental and simulated NMR data for testing and in many
examples, this data is divided into two archives, the test data set and 
additional data needed for the examples.  

.. _test_data:

The nmrglue test data sets must be downloaded and unpacked into a directory
(a directory named ``data`` under the root directory is recommended but not
required).  The conversions scripts contained in the archive must be run to
convert and process the time domain NMR data.  Additional NMR software 
(NMRPipe, etc) are requires for this processing and conversion, see the 
``README`` file in the test data archive for details.  After installing this
test data edit the ``setup.py`` file in the ``test`` directory and the
``make_links.sh`` file in the ``examples`` directory to correctly point to the 
location of the test data directory. 

.. _`nmrglue test data`: http://code.google.com/p/nmrglue/downloads/list?q=label:Test-Data

.. _example_data:

Additional data required for the nmrglue examples can be downloaded as a 
`single archive 
<http://nmrglue.googlecode.com/files/all_none_test_example_data.zip>`_.  
Unpack this archive in the ``examples`` directory.  Run the
``make_links.sh`` shell script to make symbolic links to the test data which
reused in a number of example.  On operating systems which do not support
symbolic links (Windows), the data in the test data directory will need to be 
copied by hand into the appropiate locations.  


Project Layout
--------------

The directory layout of the nmrglue project is as follows.

* ``nmrglue`` : source code for the project.

* ``doc`` : contains the setup file and source code for building the
  nmrglue documentation using `Sphinx <http://sphinx-doc.org/>`__.   

* ``tests`` : unit tests which use the  
  `nose <https://nose.readthedocs.org/en/latest/>`_ framework to verify the
  functionality of nmrglue.  See the :ref:`testing` section for details.

* ``example`` : contains numerous examples in which nmrglue is used to solve
  many real world NMR problems.


Two additional directories can be created to aid in developments.  These are
not required but will be ignored by git using the default ``.gitignore`` file 

* ``data`` : Suggested location to hold the :ref:`test data <test_data>`.

* ``sandbox`` : Suggested location to store code, data, etc not yet ready to be
  include in nmrglue.  


Suggestions
-----------

When working with the nmrglue source code please consider the following when
preparing patches.  

* Coding Style : The nmrglue source code trys to follow the 
  `PEP8 <http://www.python.org/dev/peps/pep-0008/>`_ style guide.  Consider 
  using a tool, such as `pep8 <http://pypi.python.org/pypi/pep8>`__ or 
  `pylint <http://www.logilab.org/857>`_ to check your Python code against 
  these conventions.

* Documentation : All public functions and classes should have a docstring which
  follows the `NumPy/SciPy documentation standard 
  <https://github.com/numpy/numpy/blob/master/doc/HOWTO_DOCUMENT.rst.txt>`_.  
  Private functions and classes may have shorter dostrings.  The nmrglue 
  documentation is build using `Sphinx <http://sphinx.pocoo.org/>`__.  Sphinx 
  translates `reST <http://docutils.sourceforge.net/rst.html>`_ formatted 
  documents (including docstring) into html.  When adding new functions,
  classes or parameters to nmrglue please update the docstring and make any
  necessary changes to the Sphinx files in the doc directory.
  
* :ref:`testing` : Tests are available for verifying the functionality of
  nmrglue, please include a tests when adding new functionality to the package.

* Examples : Numerous examples showing real world use of nmrglue are provided in
  the ``examples`` directory.  Contributions of additional example are welcome
  and appreciated.  


.. _testing:

Testing
-------

Tests for verifying the functionality of nmrglue are available in the test
directory.  These tests use the nose_ testing infrastructure.

.. _nose: https://nose.readthedocs.org/en/latest/


Requirements
^^^^^^^^^^^^

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
^^^^^^^^^^^^^^^^^^^^^^

After ensuring that all required packages are installed and ``setup.py`` 
correctly points to the location of the test data directory, the unit tests can
be run using the following::

    nosetests

Unit tests for a specific module can be run using::

    nosetests tests/test_pipe.py

Additional information on the `usage of the nosetests command
<https://nose.readthedocs.org/en/latest/usage.html>`_ is available.


Reporting Bugs
--------------

The preferred location for submitting feature requests and reporting bugs
is the `github issue tracker <https://github.com/jjhelmus/nmrglue/issues>`_.
Reports are also welcomed on the 
`nmrglue mailing list <http://groups.google.com/group/nmrglue-discuss>`_ or by
contacting `Jonathan Helmus <http://nmrglue.com/jhelmus>`_ directly.

Contributions
-------------

Contribution of source code or examples to nmrglue is welcomed provided the
contents can be distributed under the 
`New BSD License <http://opensource.org/licenses/BSD-3-Clause>`_.  The 
preferred method for contributing is by creating a feature branch on a github
fork of nmrglue and submitting a pull request, although patches are also
accepted.  Refer to the Numpy/SciPy 
`git workflow <http://docs.scipy.org/doc/numpy/dev/gitwash/index.html>`_ for
details on how to prepare a patch or submit a pull request.  

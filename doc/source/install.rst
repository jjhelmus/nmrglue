==================
Installation Guide
==================

Where to get nmrglue
--------------------

Install files for all platforms are available for `download
<https://github.com/jjhelmus/nmrglue/releases>`_.
The .tar.gz file should be used on Linux and OS X and a binary, .exe file for
Windows.

Requirements
------------


nmrglue depends on the following open-source packages for scientific computing
in the Python ecosystem.

+------------+------------+---------------------------------------+
| Package    | Version    | Details                               |
+============+============+=======================================+
| Python     | 3.6.0+     |                                       |
+------------+------------+---------------------------------------+
| Numpy      | 1.16+      | Required for all basic data types     |
+------------+------------+---------------------------------------+
| Scipy      | 0.16+      | Required for processing functions     |
+------------+------------+---------------------------------------+
| Matplotlib | 2.2.3+     | Optional, required for some functions |
|            |            | such as interactive phase correction  |
+------------+------------+---------------------------------------+

Additionally, an interactive environment such as `IPython <http://ipython.org/>`_, (available via several distributions such as `Jupyterlab <https://jupyterlab.readthedocs.io/en/stable/>`_, `Spyder <https://www.spyder-ide.org/>`_, `Google Colaboratory <https://colab.research.google.com/>`_, etc.) is highly recommended. An easy way of obtaining and installing these packages is to use a Python distribution which provides these packages, such as `Anaconda/Miniconda <https://www.anaconda.com/>`_. Detailed information on `installing a Scipy stack <https://scipy.org/install.html>`_ is available. The nmrglue codebase no longer supports Python 2, the last version to support Python 2.7 was 0.8.


Platform Independent Installation
---------------------------------

nmrglue is available for installation via the Python Package Index. The latest
stable version can be installed using::

    $ python -m pip install nmrglue

The current development version can be installed directly from GitHub using::

    $ python -m pip install git+git://github.com/jjhelmus/nmrglue

This requires `git` to be installed and available.


Unix/OSX Installation
---------------------

After installing the above dependencies download and extract the source
distribution and run::

    $ python setup.py install

Windows Installation
--------------------

Installation of the Scipy stack via a distribution such as `Anaconda/Miniconda`_ is recommended. nmrglue can then be installed using pip, as described above. Alternately, you can download the binary installer and run it.

Installing from source code
---------------------------

nmrglue can also be installed from source code.  See the :ref:`source_code`
section of the :ref:`development-guide` for details on this process.

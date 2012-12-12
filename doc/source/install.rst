==================
Installation Guide
==================

Where to get nmrglue
--------------------

Downloads for all platforms are available at 
`<http://code.google.com/p/nmrglue/>`_. Tar files are available for UNIX-like
systems (Linux or OSX) and a binary installer for Windows.

Requirements
------------

nmrglue requires `NumPy <http://numpy.scipy.org>`_ and 
`SciPy <http://www.scipy.org>`_ to be installed. The 
`matplotlib <http://matplotlib.org/>`_ and `IPython <http://ipython.org/>`_
packages are highly recommended.  A easy way of obtaining 
and installing these packages is to use a Python distribution which provides 
these packages, such as `EPD <http://www.enthought.com/products/epd.php>`_.  
Detailed information on 
`installing a Scipy stack <http://scipy.github.com/install.html>`_ is available.


Unix/OSX Installation
---------------------

After installing the above dependencies download and extract the source 
distribution and run::

    $ python setup.py install

Windows Installation
--------------------

Download the binary installer and run it.

Installing from source code
---------------------------

nmrglue can also be installed from source code.  See the :ref:`source_code` 
section of the :ref:`development-guide` for details on this process.

#!/usr/bin/env python

# setup script for nmrglue
#
# Copyright (C) 2010 Jonathan Helmus
# http://code.google.com/p/nmrglue
# License: BSD (See LICENSE.txt for full license)
#
# $Date$

from setuptools import setup,find_packages

cls_txt = \
"""
Intended Audience :: Science/Research
Intended Audience :: Developers
License :: OSI Approved :: BSD License
Programming Language :: Python
Topic :: Scientific/Engineering
Operating System :: Unix
Operating System :: POSIX :: Linux
"""

short_desc = "Access and modify data from common NMR files from Python"

long_desc = \
"""
The nmrglue package provides method for reading, writing, converting between, 
and processing spectral data contained in a number of common NMR file formats.

"""

package_data = {'':['doc/*.rst','doc/*.txt','doc/*.py']}

setup(
    # basic informations
    name='nmrglue',
    version = '0.1',
    packages = find_packages(),
    
    
    # metadata for upload to PyPI
    author = 'Jonathan Helmus',
    author_email = 'jjhelmus at gmail.com',
    description = short_desc,
    long_description = long_desc,
    classifiers = [x for x in cls_txt.split("\n") if x],
    url = 'http://code.google.com/p/nmrglue',
    download_url = 'http://code.google.com/p/nmrglue/downloads/list',
    requires = ['numpy (>=1.3.0)','scipy (>=0.7.1)','h5py (>=1.2.1)']
)

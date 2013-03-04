#!/usr/bin/env python

# setup script for nmrglue
#
# Copyright (C) 2010 Jonathan Helmus
# http://code.google.com/p/nmrglue
# License: BSD (See LICENSE.txt for full license)
#
# $Date$

from distutils.core import setup

setup(
    name='nmrglue',
    version='0.4',
    author = 'Jonathan J. Helmus',
    author_email = 'jjhelmus@gmail.com',
    packages=['nmrglue','nmrglue.fileio','nmrglue.analysis','nmrglue.process',
              'nmrglue.util'],
    license = 'New BSD License',
    url = 'http://code.google.com/p/nmrglue', 
    description = 'A module for working with NMR data in Python',
    long_description = open('README.txt').read(),
    requires = ['numpy','scipy'],
    classifiers = [
        'Intended Audience :: Science/Research',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
        'Topic :: Scientific/Engineering',
        'Operating System :: MacOS :: MacOS X',
        'Operating System :: Microsoft :: Windows',
        'Operating System :: POSIX :: Linux'
    ]
)

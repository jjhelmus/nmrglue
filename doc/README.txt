nmrglue documentation
=====================

This is the top level build directory for the nmrglue documentation.
nmrglue's documentation is build using sphinx which creates beautiful
documentation from reStructuredText files.

This directory contains

* _build - sphinx puts generated documentation here.

* _static - used by sphinx build system.

* Makefile - sphinx Makefile for building documentation.

* conf.py - sphinx configuration file.

* index.rst - top level nmrglue documentation ReST file.

* examples - examples documentation.

* reference - reference guide documentation

* tutorial.rst - tutorial

* plot_1d.png - image used in tutorial

* screenshot.jpg - image used in tutorial

* files_to_move.txt - list of source files to move in order that they
   be accessible to the documentation.  Can be created with:
   find ./examples/el/ -name "*.py" > files_to_move.txt 

* move_files.py - python script to move python source files to _build


To build the HTML documentation, install sphinx (0.64 or greater) and type
"make html" in this directory.  Once the sphinx building process completes
run the move_files.py  script to move the necessary source code files to 
the _build directory.  The resulting documentation will be contained in 
_build/html.

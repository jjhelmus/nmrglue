==================================
Building documentation for nmrglue
==================================

nmrglue uses Sphinx for building it's documentation.  Sphinx version 2.1.0 has
been tested, earlier and later version may work.

Instructions
------------

1.  Create and activate a environment with nmrglue dependencies, sphinx,
    numpydoc and sphinx_rtd_theme installed.  For example using conda:

        conda create -n nmrglue_dev_38 python=3.8 numpy scipy pytest sphinx pip numpydoc sphinx_rtd_theme -c conda-forge
        conda activate nmrglue_dev_38

2.  Install nmrglue in editable mode using the following from the root
    directory.

    pip install -e .

3.  Run ``make html``. in this directory.

The html documentation should be stored in the build/html directory.  


ReadTheDocs
-----------

The documentation is automatically built for each commit by ReadTheDocs and
published to https://nmrglue.readthedocs.io

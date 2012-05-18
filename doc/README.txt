==================================
Building documentation for nmrglue
==================================

nmrglue uses Sphinx for building it's documentation.  Sphinx version 1.1.3 has
been tested, earlier and later version may work.

Instructions
------------

1. Run ``make html``.

The html documentation should be stored in the build/html directory.  


Upload to gh-pages
------------------

On my machine to publish the newly build documentation to github run the
following:

rm -r ../../nmrglue_gh-pages/*
cp -r build/html/* ../../nmrglue_gh-pages/
cd ../../nmrglue_gh-pages
git add *
git commit -am "Updated documentation

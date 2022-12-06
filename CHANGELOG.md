0.10 (Unreleased)
=================



0.9 (2022-04-20)
================

FileIO
------
* Fix a bug fix passing specific procs_file in read_fid or read_pdata (#134, #141)
* Add support for reading 4D UCSF (sparky) files (#142)
* Update NMRPipe metadata keys (#144)
* Support for reading Spinsolve files (#147, #151, #155, #161)
* The converter objects can output csdm-formatted data (#152)
* nmrglue.pipe.read can read from io.Bytes objects (#153)
* nmrglue.pipe can read from Path objects (#159)
* Fix a bug in reading sparky .save files (#164)


Documentation
-------------
* Update examples with python 3 syntax (#137, #138)
* Document support for reading spinsolve files (#151)
* Spelling fixes in documentation (#153, #165)
* Example of using the peakpick.pick function (#163)
* Various other small fixes to the documentation (#162)


Analysis
--------
* Update leastsqbound to support scipy 1.8 (#167)


Misc
----
* CI framework using GitHub Actions (#132)
* Fix various test warnings (#133)
* Drop Python 2.7 support (#131)

0.8 (2020-10-29)
================

File IO
-------
* Handle data processed with nmrPipe EXT in the z-dimension (#92)
* Guess spectral widths in multi-dimensional Bruker datasets (#98)
* Allow bracket-less text values in JCAMP files (#100)
* Allow different character encodings in JCAMP files (#101)
* Only consider relevant DATATYPEs when reading JCAMP-DX files (#120)
* Support reading of Sparky .save files (#118)

Processing
----------
* Add the util/xcpy.py module for interfacing nmrglue with TopSpin (#103, #105)
* Add additional options to the autops function (#108)
* Improvements to the autops function (#124)

Documentation
-------------
* Clarify tutorial documentation (#96)
* Docstring fixes (#107)
* Document the use of nmrglue on Google Colabs (#125)
* Documentation can be built with sphinx 2.1.0 (127)
* Documentation is update at ReadTheDocs
* Various updates to the documentation (#129)
* Added a CHANGELOG.md file (#130)


Misc
----
* Replace deprecated matplotlib axisbg keyword with facecolor (#93)
* Case to int after round (#95, #106)
* Fix FutureWarnings from recent NumPy releases (#99)
* Support numpy >=1.18 (#114)

nmrglue.misc.xcpy
=================

This module is intended for use with the Bruker Topspin software.
It runs external scripts via Jython (subprocess module) that ships with Topspin.
Currently, it allows only external CPython programs to run. By default, it passes
the current folder, expno and procno to the external CPython program (if available).
For an example of how this should look like in practice, see
`PR #103 <https://github.com/jjhelmus/nmrglue/pull/103>`_.

A test script `xcpy_test.py` is provided in the same folder. This can be used to test
the setup. After copying over xcpy.py and xcpy.test to the prescribed locations (see below),
the output should look something like this:

.. image:: https://user-images.githubusercontent.com/7735086/60869506-baff0480-a24c-11e9-92d2-63d0e7fec558.gif


Usage
-----
.. code-block::

  xcpy
  xcpy [OPTIONS]
  xcpy [OPTIONS] [SCRIPTNAME]


Installation
------------
1. Copy (or symlink) this file to the following directory:
<topspin_location>/exp/stan/nmr/py/user/
2. If you now type 'xcpy' on the command line within Topspin,
this documentation should pop up
3. A configuration file needs to be written out so that xcpy
knows the location of the CPython executable and a folder where
the .py scripts are located. This can be done by 'xcpy -s'. This
can be rewritten at any time point using the same command.


Description
-----------


Options
-------
1. `-h, --help`: Brings up this document. Also brings it up when no other
option is given.

2. `-s, --settings`: Opens a dialog to write out a configuration file. The
location of the Cpython executable and the folder where all scripts
are located can be given. Use full path names for this, i.e., use
'/usr/bin/python3' for \*nix instead of 'python3' or '~/folder/python3',
and 'C:\python.exe' for Windows (note the .exe) extension. If the
configuration file already exists, the entries will be verified and
the dialog will be populated with them if these are found to be correct.
Else, an error message with the configuration file will be returned and
the dialogs will be kept empty.

3. `-n, --name`: Opens a dialog box to give the path of a script to be run

4. `-c, --config`: Prints contents of the configuration file, if it exists.
If no configuration file is found, it prints 'None'

5. `-d, --dry-run`: Prints out the command that will be executed by the subprocess
module if this option is not given.

6. `--no-args`: Does not pass any arguments (current data folder, etc) to
the external program

7. `--use-shell`: Uses shell to run the subprocess command. By default, this is
not used for \*nix, but is used for Windows. (*Warning*: this is a known
security risk and users are advised to be careful with their input while using this
option)

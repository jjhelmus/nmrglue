nmrglue.spinsolve
=================

.. automodule:: nmrglue.fileio.spinsolve

This modules is imported as nmrglue.spinsolve and can be called as such.

User Information
----------------

User Functions
^^^^^^^^^^^^^^

These are functions which are targeted for users of nmrglue.

.. autosummary::
    :toctree: generated/

    read
    guess_udic



Developer Information
---------------------

.. include:: ../../../nmrglue/fileio/spinsolve.py
    :start-line: 14
    :end-line: 37

Developer Functions
^^^^^^^^^^^^^^^^^^^

These functions are typically not used directly by users.  Developers who want
fine control over Spinsolve files may be interested in these functions.

.. autosummary::
    :toctree: generated/

    parse_spinsolve_par_line
    get_udic_from_acqu_dict
    get_udic_from_jcamp_dict

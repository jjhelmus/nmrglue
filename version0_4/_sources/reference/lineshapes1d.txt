nmrglue.lineshapes1d
====================

.. automodule:: nmrglue.analysis.lineshapes1d

This modules is imported as nmrglue.lineshapes1d and can be called as such.

Developer Functions
-------------------

These functions are typically not used directly by users.  Developers who 
want to segment spectra will be interested in these functions.

.. autosummary::
    :toctree: generated/

    sim_gauss_fwhm
    sim_lorentz_fwhm
    sim_voigt_fwhm
    sim_pvoigt_fwhm
    sim_gauss_sigma
    sim_lorentz_gamma
    sim_voigt_sigmagamma
    ls_str2class
    center_fwhm
    center_fwhm_bymoments

Developer Classes
-----------------

.. autosummary::
    :toctree: generated/

    location_scale
    gauss_sigma
    gauss_fwhm
    lorentz_gamma
    lorentz_fwhm
    location_2params
    voigt_fwhm
    voigt_sigmagamma
    pvoigt_fwhm
    scale


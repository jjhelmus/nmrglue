"""
Misc. functions
"""

import numpy as np

# default tolerences
RTOL = 1.001e-01
ATOL = 1.001e-01
DTOL = 5.001e-01


def pair_similar(dic1, data1, dic2, data2, verb=False, atol=ATOL, rtol=RTOL,
                 dtol=DTOL, ignore_pipe_display=False):
    """
    Check a dic, data pair against a second dic, data pair for differences.

    Parameters
    ----------
    dic1 : dict
        First dictionary of NMR parameters.
    data1 : ndarray
        First array of NMR data
    dic2 : dict
        Second dictionary of NMR parameters
    data2 : ndarray
        Second array of NMR data
    verb : bool, optional
        Set True for verbose reporting.
    atol : float, optional
        The absolute tolerent parameter to pass to numpy.allclose.
    rtol : float, optional
        The relative tolenance parameter to pass to numpy.allclose.

    Returns
    -------
    r1 : bool
        True is data1 and data2 are similar, False if they differ.
    r2 : bool
        True is dic1 and dic2 are similar, False if they differ.

    """
    r1 = isdatasimilar(data1, data2, verb, atol, rtol)
    r2 = isdicsimilar(dict(dic1), dict(dic2), verb, dtol,
                      ignore_pipe_display=ignore_pipe_display)
    return r1, r2


def isdatasimilar(data1, data2, verb=False, atol=ATOL, rtol=RTOL):
    """
    Check that two sets of NMR data are equal within a tolerance.

    Parameters
    ----------
    data1 : ndarray
        First array of NMR data
    data2 : ndarray
        Second array of NMR data
    verb : bool, optional
        Set True for verbose reporting.
    atol : float, optional
        The absolute tolerent parameter to pass to numpy.allclose.
    rtol : float, optional
        The relative tolenance parameter to pass to numpy.allclose.

    Returns
    -------
    r1 : bool
        True is data1 and data2 are similar, False if they differ.

    """
    r = True
    if data1.dtype != data2.dtype:
        r = False
        if verb:
            print "Dtypes do not match:", data1.dtype, data2.dtype
    if data1.shape != data2.shape:
        r = False
        if verb:
            print "Shapes do not match:", data1.shape, data2.shape
    if np.allclose(data1, data2, rtol=rtol, atol=atol) is False:
        r = False
        if verb:
            print "Data does not match"
    return r


def isitemsimilar(v1, v2, verb=False, dtol=DTOL):
    """
    Compare two values for differences

    See :py:func:`isdicsimilar` for Parameters.

    """

    r = True

    # type checking
    if type(v1) != type(v2):
        r = False
        if verb:
            print "Item has different type", type(v1), type(v2)

    # iterable checking
    elif isinstance(v1, dict):
        r = r and isdicsimilar(v1, v2, verb=verb, dtol=dtol)

    elif isinstance(v1, list):
        r = r and islistsimilar(v1, v2, verb=verb, dtol=dtol)

    # numeric type
    elif isinstance(v1, (int, float)):
        if abs(v1 - v2) > dtol:
            r = False
            if verb:
                print "Key mismatch:", v1, v2

    # all other types: just check if equal
    else:
        if v1 != v2:
            r = False
            if verb:
                print "Key mismatch:", v1, v2

    return r


def isdicsimilar(dic1, dic2, verb=False, dtol=DTOL, ignore_pipe_display=False):
    """
    Compare two dictionaries for differences

    Float and int types compared within dtol. Lists and dictionaries are
    checked recursively all other checked by simple equivalence

    Parameters
    ----------
    dic1 : dict
        First dictionary of NMR parameters.
    dic2 : dict
        Second dictionary of NMR parameters
    verb : bool, optional
        Set True for verbose reporting.
    dtol : float, optional
        Maximum allowable difference between int and float elements if dic1
        and dic2.

    Returns
    -------
    r1 : bool
        True is dic1 and dic2 are similar, False if they differ.

    """
    # create copies of the two dictionaries
    dic1 = dict(dic1)
    dic2 = dict(dic2)

    # set return value to True
    r = True

    # create sets
    kset1 = set(dic1.keys())
    kset2 = set(dic2.keys())
    dset = set.difference(kset1, kset2)
    iset = set.intersection(kset1, kset2)

    if ignore_pipe_display is True:
        iset.discard('FDMIN')
        iset.discard('FDMAX')
        iset.discard('FDDISPMIN')
        iset.discard('FDDISPMAX')
        iset.discard('FDSCALEFLAG')

    # print out any keys not in both dictionaries
    if len(dset) != 0:
        r = False
        if verb:
            print "Keys not in both dictionaries:", dset

    # loop over keys in both sets
    for k in iset:
        v1, v2 = dic1[k], dic2[k]

        if not isitemsimilar(v1, v2, verb=verb, dtol=dtol):
            print "For key:", k
            r = False

    return r


def islistsimilar(l1, l2, verb=False, dtol=DTOL):
    """
    Compare two lists (or iterable) for differences

    See :py:func:`isdicsimilar` for Parameters.

    """
    # set return value to True
    r = True

    # print out any keys not in both dictionaries
    if len(l1) != len(l2):
        r = False
        if verb:
            print "Lists not of same length:", len(l1), len(l2)

    # loop over keys in both sets
    for v1, v2 in zip(l1, l2):
        if not isitemsimilar(v1, v2, verb=verb, dtol=dtol):
            r = False

    return r

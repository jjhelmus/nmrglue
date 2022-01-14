"""
analysisbase provides general purpose analysis functions and classes used by
several nmrglue.analysis modules
"""

import numpy as np

pi = np.pi


# helper functions
def neighbors(pt, shape, structure):
    """
    Generate a list of all neighbors to a point.

    Parameters
    ----------
    pt : tuple of ints
        Index of the point to find neighbors of.
    shape : tuple of ints
        Shape of the region.
    structure : ndarray of bools
        Structure element that defines connections.

    Returns
    -------
    pts : list of int tuples
        List of tuples which represent indices for all points neighboring pt.
        Edges are treated as stopping points.

    """
    # set middle of structure to False
    s = np.copy(structure)  # copy structure
    middle = [int(np.floor(i / 2.)) for i in s.shape]  # find middle
    s.flat[np.ravel_multi_index(middle, s.shape)] = False
    offsets = np.argwhere(s) - middle

    # loop over the offset adding all valid points
    pts = []
    for offset in offsets:
        npt = pt - offset
        if valid_pt(npt, shape):
            pts.append(tuple(npt))
    return pts


def valid_pt(pt, shape):
    """
    Determine if a point (indices) is valid for a given shaped
    """
    for i, j in zip(pt, shape):
        if i < 0:   # index is not negative
            return False
        if i >= j:    # index is less than j
            return False
    return True

dimension_names = ['A', 'Z', 'Y', 'X']


# utility functions
def find_limits(pts):
    """
    Find the limits which outline the provided list of points

    Parameters
    ----------
    pts : list of int tuples
        List of points [(z0, y0, x0), (z1, y1, x1), ...]

    Returns
    -------
    min : ndarray
        Array of minimum indices: array([zmin, ymin, xmin]
    max : ndarray
        Array of maximum indices: array([zmin, ymin, xmin]

    See Also
    --------
    limits2slice : Create a list of slices from min, max limits

    """
    arr_pts = np.array(pts)
    return np.min(arr_pts, 0), np.max(arr_pts, 0)


def limits2slice(limits):
    """
    Create a set of slice objects given an array of min, max limits.

    Parameters
    ----------
    limits: tuple, (ndarray, ndarray)
        Two tuple consisting of array of the minimum and maximum indices.

    Returns
    -------
    slices : list
        List of slice objects which return points between limits

    See Also
    --------
    find_limits : Find the minimum and maximum limits from a list of points.
    slice2limits : Find a minimum and maximum limits for a list of slices.

    """
    mins, maxs = limits
    return tuple([slice(i, j + 1) for i, j in zip(mins, maxs)])


def slice2limits(slices):
    """
    Create a tuple of minimum, maximum limits from a set of slices.

    Parameters
    ----------
    slices : list
        List of slice objects which return points between limits

    Returns
    -------
    limits: tuple, (ndarray, ndarray)
        Two tuple consisting of array of the minimum and maximum indices.

    See Also
    --------
    limits2slice : Find a list of slices given minimum and maximum limits.
    """
    mins = [s.start for s in slices]
    maxs = [s.stop - 1 for s in slices]
    return mins, maxs


def squish(r, axis):
    """
    Squish array along an axis.

    Determine the sum along all but one axis for an array.

    Parameters
    ----------
    r : ndarray
        Array to squish.
    axis : int
        Axis of r to squish along.

    Returns
    -------
    s : 1D ndarray
        Array r squished into a single dimension.

    """
    # put axis to be squished as the last axis
    N = int(r.ndim)
    r = r.swapaxes(axis, N - 1)

    # sum along leading axis N-1 times
    for i in range(N - 1):
        r = r.sum(0)
    return r


# Windowing classes
class ndwindow(object):
    """
    An N-dimensional iterator to slice arrays into windows.

    Given the shape of an array and a window size, an 'ndwindow' instance
    iterators over tuples of slices which slice an the array into wsize
    sub-arrays.  At each iteration, the index of the center of the sub-array
    is incremented by one along the last dimension.  Array borders are ignored
    so the resulting sub-array can be smaller than wsize.  If wsize contains
    even values the window is off center containing an additional point with
    lower index.

    Parameters
    ----------
    size : tuple of ints
        Size of array to generate tuples of slices from.
    wsize : tuple of ints
        Window/sub-array size. Size of the area to select from array.  This is
        the maximum size of the window.

    Examples
    --------

    >>> a = np.arange(12).reshape(3,4)
    >>> for s in ndwindow(a.shape,(3,3)):
    ...     print(a[s])
    [[0 1]
     [4 5]]
    [[0 1 2]
     [4 5 6]]
    [[1 2 3]
     [5 6 7]]
    [[2 3]
     [6 7]]
    [[0 1]
     [4 5]
     [8 9]]
    [[ 0  1  2]
     [ 4  5  6]
     [ 8  9 10]]
    [[ 1  2  3]
     [ 5  6  7]
     [ 9 10 11]]
    [[ 2  3]
     [ 6  7]
     [10 11]]
    [[4 5]
     [8 9]]
    [[ 4  5  6]
     [ 8  9 10]]
    [[ 5  6  7]
     [ 9 10 11]]
    [[ 6  7]
     [10 11]]

    See Also
    --------
    ndwindow_index : Iterator of a ndwindow and index of the window center
    ndwindow_inside : Iterator over equal sized windows in the array.

    """
    def __init__(self, shape, wsize):
        """ Set up the ndwindow object """
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(shape)
        wsize = np.array(wsize)
        self.sub = np.ceil((wsize - 1.) / 2.)
        self.add = wsize - 1. - self.sub

    def __next__(self):
        """ next iterator. """
        return self.next()

    def next(self):
        """ x.next() -> the next value, or raise StopIteration """
        center = self.ndindex.next()
        start = [max(0, i - j) for i, j in zip(center, self.sub)]
        stop = [i + j + 1 for i, j in zip(center, self.add)]
        return tuple([slice(x, y) for x, y in zip(start, stop)])

    def __iter__(self):
        """ x.__iter__() <==> iter(x) """
        return self


class ndwindow_index(object):
    """
    An N-dimensional iterator object which returns the index of the window
    center and a :py:class:`ndwindow` slice array.  See :py:class:`ndwindow`
    for additional documentation.

    This class is equivalent to:

    for slices, index in zip(np.ndindex(shape), ndwindow(shape,wshape)):
        return (index, slice)

    See Also
    --------
    ndwindow: Iterator over only the window slices.
    ndwindow_inside : Iterator over equal sized windows in the array.

    """
    def __init__(self, shape, wsize):
        """ Set up the object """
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(shape)
        wsize = np.array(wsize)
        self.sub = np.ceil((wsize - 1.) / 2.)
        self.add = wsize - 1. - self.sub

    def __next__(self):
        """ next iterator. """
        return self.next()

    def next(self):
        """ x.next() -> the next value, or raise StopIteration """
        center = self.ndindex.next()
        start = [max(0, i - j) for i, j in zip(center, self.sub)]
        stop = [i + j + 1 for i, j in zip(center, self.add)]
        return center, tuple([slice(x, y) for x, y in zip(start, stop)])

    def __iter__(self):
        """ x.__iter__() <==> iter(x) """
        return self


class ndwindow_inside(object):
    """
    An N-dimensional iterator to slice arrays into uniform size windows.

    Given the shape of an array and a window size, an 'ndwindow_inside'
    instance iterators over tuples of slices which slice an the array into
    uniform size wsize windows/sub-arrays.  At each iteration, the index of
    the top left of the sub-array is incremented by one along the last
    dimension until the resulting windows would extend past the array border.
    All sub-arrays are equal sized (wsize).

    Parameters
    ----------
    size : tuple of ints
        Size of array to generate tuples of slices from.
    wsize : tuple of ints
        Size of the area to select from array (widow size).

    Examples
    --------

    >>> a = np.arange(9).reshape(3,3)
    >>> for s in ndwindow_inside(a.shape,(2,2)):
    ...     print(a[s])
    [[0 1]
     [3 4]]
    [[1 2]
     [4 5]]
    [[3 4]
     [6 7]]
    [[4 5]
     [7 8]]

    See Also
    --------
    ndwindow : Iterator over non-uniform windows.
    ndwindow_inside_index : Iterator of a ndwindow_inside and the index of the
        window's top left point.

    """
    def __init__(self, shape, wsize):
        """ Set up the object """
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(
            tuple(np.array(shape) - np.array(wsize) + 1))
        self.wsize = wsize

    def __next__(self):
        """ next iterator. """
        return self.next()

    def next(self):
        """ x.next() -> the next value, or raise StopIteration """
        start = self.ndindex.next()
        stop = np.array(start) + np.array(self.wsize)
        return tuple([slice(x, y) for x, y in zip(start, stop)])

    def __iter__(self):
        """ x.__iter__() <==> iter(x) """
        return self


class ndwindow_inside_index(object):
    """
    An N-dimensional iterator object which returns the index of the window
    top-left and a :py:class:`ndwindow_inside` slice array.

    Similar to :py:class:`ndwindow_index` but reports top left index of
    window.

    See :py:class:`ndwindow_inside` and :py:class`ndwindow_index` for addition
    documentation.

    """
    def __init__(self, shape, wsize):
        " Set up the object """
        if len(shape) != len(wsize):
            raise ValueError("shape and wsize do match match")
        self.ndindex = np.ndindex(
            tuple(np.array(shape) - np.array(wsize) + 1))
        self.wsize = wsize

    def __next__(self):
        """ next iterator. """
        return self.next()

    def next(self):
        """ x.next() -> the next value, or raiseStopIteration """
        start = self.ndindex.next()
        stop = np.array(start) + np.array(self.wsize)
        return (start, tuple([slice(x, y) for x, y in zip(start, stop)]))

    def __iter__(self):
        """ x.__iter__() <==> iter(x) """
        return self

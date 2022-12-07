"""
nmrglue table functions.

nmrglue uses numpy records array as stores of various data (peak tables,
trajectories, etc).  This module provides functions to read and write records
arrays from disk.  Formatting of the numeric values is left to Python's str
function and only the data type need be specified.  In addition this module
contains functions to convert nmrglue's table format NMRPipe's table format.
"""

import numpy as np
from . import fileiobase


def pipe2glue(pcomments, pformat, rec):
    """
    Convert a NMRPipe table to a nmrglue table

    Parameters
    ----------
    pcomments : list
        List of NMRPipe comment lines.
    pformats : list
        List of NMRPipe table column formats strings.
    rec : recarray
        Records array with named fields.

    Returns
    -------
    comments : list
        List of comments
    rec : recarray
        Records array with named fields.

    """
    # add a "#" to the list of comments and we are done
    comments = ["# " + c for c in pcomments]
    return comments, rec


def glue2pipe(comments, rec):
    """
    Convert a nmrglue table to a NMRPipe table.

    Parameters
    ----------
    comments : list
        List of comments
    rec : recarray
        Records array with named fields.

    Returns
    -------
    pcomments : list
        List of NMRPipe comment lines.
    pformats : list
        List of NMRPipe table column formats strings.  This list is guessed
        from the data types and precision in the reconrds array.  This may not
        be the exact format desired, edit this to your liking.
    rec : recarray
        Records array with named fields.

    """
    # add REMARK to each comment
    pcomments = ["REMARK " + c for c in comments]

    # guess the pipe format strings
    pformat = [guess_pformat(rec[t]) for t in rec.dtype.names]

    return pcomments, pformat, rec


def guess_pformat(col):
    """
    Guess a NMRPipe table column format string given a column.

    Parameters
    ----------
    col : ndarray
        Array from a records array.

    Returns
    -------
    s : str
        String for formatting NMRPipe table.

    """
    kind = col.dtype.kind

    if kind in ('S', 'a'):  # string
        return '%s'

    if kind in ('i', 'u'):  # integer (signed or unsigned)
        # N is the number of digits in largest value, or 1
        N = max(np.ceil(np.log(np.abs(col).max()) / np.log(10)), 1)
        # +1 for sign
        return f'%{int(N) + 1}d'

    if kind == 'f':
        # will be either %+e or %N.3f, see if 'e' is in %g to determine
        if True in ['e' in '%g' % (v) for v in col]:
            return '%+e'
        else:
            N = max(np.ceil(np.log(np.abs(col).max()) / np.log(10)), 1)
            # +1 for sign, +1 for decimal points, +3 for precision
            return f'%{int(N) + 5}.3f'

    # remaining kinds: 'c' - complex, 'b' - boolean, 'U' - unicode, 'V' - void
    raise ValueError("unknown kind %s in column" % (kind))


def read(filename):
    """
    Read a nmrglue table file.

    Parameters
    ----------
    filename : str
        Filename of nmrglue table file to read.

    Returns
    -------
    comments : list
        List of comments (strings terminated with newline)
    rec : recarray
        Records array with named fields.

    """
    # pull out the comment lines from the file (start with #)
    f = open(filename, 'rb')
    comments = [l for l in f if l[0] == '#']
    f.close()

    # find the line beginning with # NAMES and parse out the column names
    nl = [i for i, l in enumerate(comments) if l[:7] == "# NAMES"]
    if len(nl) != 1:
        raise OSError("%s does not have a # NAMES line" % (filename))
    dtd = {'names': comments.pop(nl[0])[7:].split()}

    # find the line beginning with # DTYPE and parse out the column names
    dl = [i for i, l in enumerate(comments) if l[:9] == "# FORMATS"]
    if len(dl) != 1:
        raise OSError("%s does not have a # FORMATS line" % (filename))
    dtd['formats'] = comments.pop(dl[0])[9:].split()

    # return the data as a records array
    return comments, np.atleast_1d(np.recfromtxt(filename, dtype=dtd))


def write(filename, comments, rec, overwrite=False):
    """
    Write a nmrglue table to file.

    Parameters
    ----------
    filename : str
        Filename of file to write table to.
    comments : list
        List of comments (strings terminated with newline).
    rec : recarray
        Records array to write to file.
    overwrite : bool, optional
        True to overwrite file if it exists. False will raise an Warning if the
        file exists.

    """
    # open the file for writing
    f = fileiobase.open_towrite(filename, overwrite)

    # write out the comment lines at the top of the file
    for c in comments:
        f.write(c)

    # Determine the list of column names
    names = rec.dtype.names

    # Determine the list of column formats
    sizes = [rec.dtype[n].itemsize for n in names]
    kinds = [rec.dtype[n].kind for n in names]
    formats = [k + str(s) for s, k in zip(sizes, kinds)]

    # write out the NAMES and FORMATS lines
    f.write("# NAMES " + " ".join(names) + "\n")
    f.write("# FORMATS " + " ".join(formats) + "\n")

    # maximum string length for each column
    col_len = [max([len(str(i)) for i in rec[n]]) for n in names]

    # write out each line of the table
    for row in rec:
        s = " ".join([str(v).ljust(l) for v, l in zip(row, col_len)])
        f.write(s + "\n")
    f.close()


# Row functions
def insert_row(rec, N, row):
    """
    Insert a row into a records array before row number N.

    Parameters
    ----------
    rec : recarray
        Records array.
    N : int
        Row number to insert new row before.
    row : array_like
        Array or similar object which will be converted into a new row.

    Returns
    -------
    new_rec : recarray
        New records array with inserted row.

    """
    return np.insert(rec, N, row)


def append_row(rec, row):
    """
    Append a row to the end of a records array.

    Parameters
    ----------
    rec : recarray
        Records array.
    row : array_like
         Array or similar object which will be converted into a new row.

    Returns
    -------
    new_rec : recarray
        New records array with inserted row.

    """
    N = len(rec)
    return insert_row(rec, N, row)


def delete_row(rec, N):
    """
    Delete a row from a records array.

    Parameters
    ----------
    rec : recarray
        Records array.
    N : int
        Row number to delete.

    Returns
    -------
    new_rec : recarray
        New records array with row deleted.

    See Also
    --------
    reorder_rows : delete multiple rows in a single call.

    """
    return np.delete(rec, N)


def reorder_rows(rec, new_order):
    """
    Reorder or delete rows in a records array.

    This function can also be used to delete multiple rows from a records
    array, only the rows in the new_order list are retained in the new records
    array.

    Parameters
    ----------
    rec : recarray
        Records array.
    new_order : list
        List of row indices and order in new records array.  Only the rows in
        this list are retained in the new records array.  Therefore this
        function can also be used to delete multiple rows from a records
        array.

    Returns
    -------
    new_rec : recarray
        New records array with rows reordered.

    """
    return np.take(rec, new_order)


# Column functions
def append_column(rec, col, name=None, format=None):
    """
    Append a column to the end of a records array.

    Parameters
    ----------
    rec : recarray
        Records array.
    col : array_like
        Array or similar object which will be converted into the new column.
    name : str, optional
        Name of the column. If None col.dtypes.name will be used.
    format : dtype, optional
        Data type to convert the new column into before appending. Required if
        col is not an ndarray.

    Returns
    -------
    new_rec : recarray
        New records array with column appended.

    """
    N = len(rec.dtype.descr)
    return insert_column(rec, N, col, name, format)


def insert_column(rec, N, col, name=None, format=None):
    """
    Insert a column into a records array.

    Parameters
    ----------
    rec : recarray
        Records array.
    col : array_like
        Array or similar object which will be converted into the new column.
    N : int
        Column number to insert new column before.
    name : str, optional
        Name of the column. If None col.dtypes.name will be used.
    format : dtype, optional
        Data type to convert the new column into before appending. Required if
        col in not an ndarray.

    Returns
    -------
    new_rec : recarray
        New records array with column inserted.

    """
    col = np.array(col)

    # get name and format parameter from column if not provided
    if name is None:
        if col.dtype.names is not None:
            name = col.dtype.names
        else:
            raise ValueError("Must provide a name for the column")
    if format is None:
        format = col.dtype.str

    # insert the new column (name, format) to the table dtypes
    dtd = rec.dtype.descr
    dtd.insert(N, (name, format))

    # create the new table with an additional column
    new_rec = np.empty(rec.shape, dtd)

    # fill in the old columns
    for n in rec.dtype.names:
        new_rec[n] = rec[n]
    # and the new column
    new_rec[name] = col.astype(format)

    return np.atleast_1d(np.rec.array(new_rec))


def delete_column(rec, N):
    """
    Delete a column from a records array.

    Parameters
    ----------
    rec : recarray
        Records array.
    N : int
        Column number to delete.

    Returns
    -------
    new_rec : recarray
        New records array with column deleted.

    See Also
    --------
    reorder_columns : Delete multiple columns from a records array.

    """
    # remove the column from the list of columns.
    dtd = rec.dtype.descr
    dtd.pop(N)

    # create the new records array and fill it in
    new_rec = np.empty(rec.shape, dtd)
    for n in new_rec.dtype.names:
        new_rec[n] = rec[n]
    return np.atleast_1d(np.rec.array(new_rec))


def reorder_columns(rec, new_order):
    """
    Reorder or delete columns in a records array.

    Parameters
    ----------
    rec : recarray
        Records array.
    new_order : list
        List of column indices and order in new records array.  Only the
        columns in this list are retained in the new records array.
        Therefore this function can also be used to delete multiple columns
        from a records array.

    Returns
    -------
    new_rec : recarray
        New records array with columns reordered.

    """
    # reorder the dtype description list
    dtd = rec.dtype.descr
    new_dtd = [dtd[i] for i in new_order]

    # create the new array and fill it in.
    new_rec = np.empty(rec.shape, new_dtd)
    for n in new_rec.dtype.names:
        new_rec[n] = rec[n]
    return np.atleast_1d(np.rec.array(new_rec))

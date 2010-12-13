""" 
nmrglue table functions 

nmrglue uses numpy records array as stores of various data (peak tables, 
trajectories,etc).  This module provides functions to read and write records
arrays from disk.  Formatting of the numeric values is left to Python's str
function and only the data type need be specified.  In addition this module
contains functions to convert nmrglue's comments to and from NMRPipe's 
pcomments,pformats lists.

"""

import numpy as np
import fileiobase


def pipe2glue(pcomments,pformat,rec):
    """
    Convert a NMRPipe table to nmrglue table
    
    Parameters:

    * pcomments List of NMRPipe comment lines.
    * pformats  List of NMRPipe table column formats strings.
    * rec       Records array with named fields.
    
    Returns: comments,rec

    * comments  List of comments
    * rec       Records array with named fields.

    """
    # add a "#" to the list of comments and we are done
    comments = ["# "+c for c in pcomments]
    return comments,rec
    

def glue2pipe(comments,rec):
    """
    Convert a nmrglue table to a NMRPipe table

    The pformats list is a guess from data type and precision in the records
    array.  You may want to edit this to your liking.

    Parameters:

    * comments  List of comments
    * rec       Records array with named fields.

    Returns: pcomments,pformat,rec

    * pcomments List of NMRPipe comment lines.
    * pformats  List of NMRPipe table column formats strings.
    * rec       Records array with named fields.

    """
    # add REMARK to each comment
    pcomments = ["REMARK "+c for c in comments]

    # guess the pipe format strings
    pformat = [guess_pformat(rec[t]) for t in rec.dtype.names]

    return pcomments,pformat,rec


def guess_pformat(col):
    """
    Guess a NMRPipe table column format string given a column

    Parameters:

    * col   Array from records array

    Returns string for formatting NMRPipe table

    """
    kind = col.dtype.kind
    
    if kind=='S' or kind=='a':  # string
        return '%s'
    
    if kind=='i' or kind=='u':  # integer (signed or unsigned)
        # N is the number of digits in largest value, or 1
        N = max(np.ceil(np.log(np.abs(col).max())/np.log(10)),1)
        # +1 for sign
        return '%{0}d'.format(int(N+1)) 

    if kind=='f':
        # will be either %+e or %N.3f, see if 'e' is in %g to determine
        if True in ['e' in '%g'%(v) for v in col]:
            return '%+e'
        else:
            N = max(np.ceil(np.log(np.abs(col).max())/np.log(10)),1)
            # +1 for sign, +1 for decimal points, +3 for precision  
            return '%{0}.3f'.format(int(N+5))

    # remaining kinds: 'c' - complex, 'b' - boolean, 'U' - unicode, 'V' - void
    raise ValueError("unknown kind %s in column"%(kind))


def read(filename):
    """
    Read a table (.tbl) file.

    Parameters:

    * filename  Name of table file to read 

    Returns: (comments,rec)

    * comments  List of comments (strings terminated with newline)
    * rec       Records array with named fields.

    """
    # pull out the comment lines from the file (start with #)
    f = open(filename,'r')
    comments = [l for l in f if l[0]=='#']
    f.close()
    
    # find the line beginning with # NAMES and parse out the column names
    nl = [i for i,l in enumerate(comments) if l[:7]=="# NAMES"]
    if len(nl)!=1:
        raise IOError("%s does not have a # NAMES line"%(filename))
    dtd = {'names':comments.pop(nl[0])[7:].split()}
    
    # find the line beginning with # DTYPE and parse out the column names
    dl = [i for i,l in enumerate(comments) if l[:9]=="# FORMATS"]
    if len(dl)!=1:
        raise IOError("%s does not have a # FORMATS line"%(filename))
    dtd['formats'] = comments.pop(dl[0])[9:].split()
    
    # return the data as a records array
    return comments,np.atleast_1d(np.recfromtxt(filename,dtype=dtd))


def write(filename,comments,rec,overwrite=False):
    """
    Write a nmrglue table to file (.tbl).

    Parameters:
    
    * filename  Name of table file to write.
    * comments  List of comments (strings terminated with newline).
    * rec       Records array to write to file.
    * overwrite Set True to overwrite file if it exists.

    """
    # open the file for writing
    f = fileiobase.open_towrite(filename,overwrite)
    
    # write out the comment lines at the top of the file
    for c in comments:
        f.write(c)
    
    # Determind the list of column names
    names = rec.dtype.names
    
    # Determind the list of column formats
    sizes = [rec.dtype[n].itemsize for n in names]
    kinds = [rec.dtype[n].kind for n in names]
    formats = [k+str(s) for s,k in zip(sizes,kinds)]
    
    # write out the NAMES and FORMATS lines
    f.write("# NAMES "+" ".join(names)+"\n")
    f.write("# FORMATS "+" ".join(formats)+"\n")

    # maximum string length for each column
    col_len = [max([len(str(i)) for i in rec[n]]) for n in names]

    # write out each line of the table
    for row in rec:
        s = " ".join([str(v).ljust(l) for v,l in zip(row,col_len)])
        f.write(s+"\n")

    f.close()

# Row functions (these are easy)

def insert_row(rec,N,row):
    """ 
    Insert a row into a records array before row number N.
    
    Parameters:

    * rec   Records array.
    * N     Row number to insert new row before, integer.
    * row   Tuple, etc which will be converted into a new row.
    
    Returns: new_rec (new records array with inserted row)

    """
    return np.insert(rec,N,row)


def append_row(rec,row):
    """ 
    Append a row to the end of a records array
    
    Parameters:

    * rec   Records array.
    * row   Tuple, etc which will be converted into a new row.

    Returns: new_rec (new records array with appeneded row)

    """
    N = len(rec)
    return insert_row(rec,N,row)


def delete_row(rec,N):
    """ 
    Delete row N from records array.
    
    Parameters:

    * rec   Records array.
    * N     Row number to delete, integer.
    
    Returns: new_rec (new records array with row deleted)

    Use reorder_rows to delete multiple rows in a single call.

    """
    return np.delete(rec,N)


def reorder_rows(rec,new_order):
    """
    Reorder rows in a records array.

    This function can also be used to delete multiple rows from a records 
    array, only the rows in the new_order list are retained in the new records
    array.

    Parameters:
    
    * rec   Records array.
    * new_order List of row indices and order in new records array.

    Returns: new_rec (new records array with row reordered)

    """
    return np.take(rec,new_order)


# Column functions (harder)

def append_column(rec,col,name=None,format=None):
    """ 
    Append a column to the end of a records array.
    
    Parameters:

    * rec       Records array.
    * col       Array which will be converted into the new column.
    * name      Name of the column (if Name must be given in col.dtype.names)
    * format    Data type to convert the new column into before appending.

    Returns: new_rec (new records array with column appended)

    """
    N = len(rec.dtype.descr)
    return insert_column(rec,N,col,name,format)


def insert_column(rec,N,col,name=None,format=None):
    """ 
    Insert a column into a records array before column number N.

    Parameters:

    * rec       Records array.
    * col       Array which will be converted into the new column.
    * N         Number of the column to insert new column before.
    * name      Name of the column (if Name must be given in col.dtype.names)
    * format    Data type to convert the new column into before appending.

    Returns: new_rec (new records array with column inserted)

    """
    col = np.array(col)
    
    # get name and format parameter from column if not provided
    if name==None:
        if col.dtype.names!=None:
            name=col.dtype.names
        else:
            raise ValueError("Must provide a name for the column")
    if format==None:
        format = col.dtype.str
    
    # insert the new column (name, format) to the table dtypes
    dtd = rec.dtype.descr
    dtd.insert(N,(name,format))

    # create the new table with an additional column
    new_rec = np.empty( rec.shape,dtd)
   
    # fill in the old columns
    for n in rec.dtype.names: 
        new_rec[n]=rec[n]
    # and the new column
    new_rec[name] = col.astype(format)
    
    return np.atleast_1d(np.rec.array(new_rec))


def delete_column(rec,N):
    """ 
    Delete a column from a records array
    
    Parameters:

    * rec   Records array.
    * N     Number of the column to delete.

    Returns: new_rec (new records array with column inserted)
    
    """
    # remove the column from the list of columns.
    dtd = rec.dtype.descr
    dtd.pop(N)

    # create the new records array and fill it in
    new_rec = np.empty(rec.shape,dtd)
    for n in new_rec.dtype.names:
        new_rec[n] = rec[n]
    return np.atleast_1d(np.rec.array(new_rec))


def reorder_columns(rec,new_order):
    """ 
    Reorder columns in a records array
    
    Parameters:

    * rec           Records array.
    * new_order     List of column indices to order in new records array.

    Returns: new_rec (new records array with columns reordered)

    """
    # reorder the dtype description list
    dtd = rec.dtype.descr
    new_dtd = [dtd[i] for i in new_order]

    # create the new array and fill it in.
    new_rec = np.empty(rec.shape,new_dtd)
    for n in new_rec.dtype.names:
        new_rec[n] = rec[n]
    return np.atleast_1d(np.rec.array(new_rec))


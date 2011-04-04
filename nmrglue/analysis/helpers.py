"""
Helper functions
"""

from ..fileio import table

# bound functions

def delta(v,d):
    """
    For limits of v +/- d for array v
    """
    return zip( v-d,v+d)

def limit(vmin,vmax,npeaks):
    """
    Fixed limits
    """
    return [(vmin,vmax)]*npeaks

def scale_params(rec,prefix,first,last):
    """
    scale lineshape parameters in rec with name prefix+[0,1,2,...]
    """
    #return zip( *[rec[prefix+str(z)] for z in range(first,last+1)] )
    return [rec[prefix+str(z)] for z in range(first,last+1)]

def no_limits(nparams,npeaks):
    """
    """
    return [[(None,None)]*npeaks]*nparams

def no_limits_amp(npeaks):
    return [(None,None)]*npeaks

# misc

def super_zip(l):
    return zip(*[zip(*i) for i in l])

# misc

def scale_columns(prefix,first,last):
    """
    """
    return [prefix+str(i) for i in range(first,last+1)]


# table packing/unpacking

def add_to_table(rec,columns,column_names):
    for col,col_name in zip(columns,column_names):
        rec = table.append_column(rec,col,name=col_name)
    return rec

def pack_table(pbest,abest,iers,rec,param_columns,amp_column,ier_column=None):

    # pack the amplitudes
    rec[amp_column] = abest

    # pack the parameters
    for dbest,dcolumns in zip(zip(*pbest),param_columns):
        for p,c in zip(zip(*dbest),dcolumns):
            rec[c] = p

    # pack the iers
    if ier_column!=None:
        rec[ier_column] = iers

def unpack_table(rec,param_columns,amp_column):
    
    params = zip( *[zip(*[rec[c] for c in dc]) for dc in param_columns])
    amps = rec[amp_column]
    return params,amps

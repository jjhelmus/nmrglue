"""
Helper functions for lineshape fitting.
"""


# bound functions

def delta(v,d):
    """
    Limits of v +/- d for array v and d
    """
    return zip( v-d,v+d)

def limit(vmin,vmax,npeaks):
    """
    Fixed limits, vmin to vmax for all peaks.
    """
    return [(vmin,vmax)]*npeaks

def scale_params(rec,prefix,first,last):
    """
    scale lineshape parameters in rec with name prefixX with X from from 
    first to last (inclusive)
    """
    #return zip( *[rec[prefix+str(z)] for z in range(first,last+1)] )
    return [rec[prefix+str(z)] for z in range(first,last+1)]

def no_limits(nparams,npeaks):
    """
    No limits on nparameters for npeaks.
    """
    return [[(None,None)]*npeaks]*nparams

def no_limits_amp(npeaks):
    """
    No limits for amplitudes for npeaks
    """
    return [(None,None)]*npeaks

# misc

def super_zip(l):
    """
    zip a list after zipping each element.  Useful for lineshape parameter and
    bounds.
    """
    return zip(*[zip(*i) for i in l])

# misc

def scale_columns(prefix,first,last):
    """
    Create a list of scale parameter column names with name prefixX with X
    from first to last (inclusive)
    """
    return [prefix+str(i) for i in range(first,last+1)]



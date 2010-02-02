"""
Misc. functions
"""

import numpy as np

# default tolerences
rtol = 1.001e-01
atol = 1.001e-01
dtol = 1.001e-01

def pair_similar(dic1,data1,dic2,data2,verb=False,atol=atol,rtol=rtol,
    dtol=dtol):
    """
    Check a dic,data pair against a second dic,data pair for differences

    Returns tuple of Booleans indicating agreement between data and 
    dictionaries
    """
    r1 = isdatasimilar(data1,data2,verb,atol,rtol)
    r2 = isdicsimilar(dict(dic1),dict(dic2),verb,dtol)

    return r1,r2

def isdatasimilar(data1,data2,verb=False,atol=atol,rtol=rtol):
    """
    Check that data is equal within tolerance
    """

    r = True
    if data1.dtype != data2.dtype:
        r = False
        if verb:
            print "Dtypes do not match",data1.dtype,data2.dtype
    if data1.shape != data2.shape:
        r = False
        if verb:
            print "Shapes do not match",data1.shape,data2.shape
    if np.allclose(data1,data2,rtol=rtol,atol=atol) == False:
        r = False
        if verb:
            print "Data does not match"
    return r

def isdicsimilar(dic1,dic2,verb=False,dtol=dtol):
    """
    Compare two dictionaries for differences

    float and int types compared within dtol
    lists and dictionaries checked recursively
    all other checked by simple equivalence

    """

    # create copies of the two dictionaries
    dic1 = dict(dic1)
    dic2 = dict(dic2)
    
    # set return value to True
    r = True

    # create sets
    kset1 = set(dic1.keys())
    kset2 = set(dic2.keys())
    dset = set.difference(kset1,kset2)
    iset = set.intersection(kset1,kset2)

    # print out any keys not in both dictionaries
    if len(dset) !=0:
        r = False
        if verb:
            print "Keys not in both dictionaries:",dset

    # loop over keys in both sets
    for k in iset:
        v1,v2 = dic1[k],dic2[k]

        # type checking
        if type(v1) != type(v2):
            r=False
            if verb:
                print "Key has different type",k,type(dic1[k]),type(dic2[k])
            
        # iterable checking
        if isinstance(v1,dict):
            #if verb:
            #    print "Checking sub-dictionary:",k
            r = r and isdicsimilar(v1,v2,verb=verb,dtol=dtol)
        
        elif isinstance(v1,list):
            #if verb:
            #    print "Checking sub-list:",k
            r = r and islistsimilar(v1,v2,verb=verb,dtol=dtol)

        # numeric type
        elif isinstance(v1,int) or isinstance(v1,float):
            if abs(v1-v2) > dtol:
                r = False
                if verb:
                    print "Key mismatch:",k,v1,v2

        # all other type just check if equal
        else:
            if v1!=v2:
                r = False
                if verb:
                    print "Key mismatch:",k,v1,v2

    return r

def islistsimilar(l1,l2,verb=False,dtol=dtol):
    """
    Compare two lists (or iterable) for differences

    see isdicsimilar

    """

    # set return value to True
    r = True

    # print out any keys not in both dictionaries
    if len(l1) != len(l2):
        r = False
        if verb:
            print "Lists not of same length:",len(l1),len(l2)

    # loop over keys in both sets
    for v1,v2 in zip(l1,l2):
        
        # type checking
        if type(v1) != type(v2):
            r=False
            if verb:
                print "Item has different type",v1,v2
            
        # iterable checking
        if isinstance(v1,dict):
            #if verb:
            #    print "Checking sub-dictionary"
            r = r and isdicsimilar(v1,v2,verb=verb,dtol=dtol)
        
        if isinstance(v1,list):
            #if verb:
            #    print "Checking sub-list"
            r = r and islistsimilar(v1,v2,verb=verb,dtol=dtol)

        # numeric type
        elif isinstance(v1,int) or isinstance(v1,float):
            if abs(v1-v2) > dtol:
                r = False
                if verb:
                    print "Item mismatch:",v1,v2

        # all other type just check if equal
        else:
            if v1!=v2:
                r = False
                if verb:
                    print "Item mismatch:",v1,v2

    return r

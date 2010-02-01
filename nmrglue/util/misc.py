"""
Misc. functions
"""

import numpy as np


rtol = 1.0000000000000001e-01   # 1.0001E-6 is the goal, 1.e-3 is good
atol = 0.001
dtol = 0.001

def check_pair(dic1,data1,dic2,data2,verb=False):

    r1 = isdataclose(data1,data2,verb)
    r2 = isdicsame(dict(dic1),dict(dic2),verb)

    return r1,r2

def isdataclose(data1,data2,verb=False):

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

def isdicsame(dic1,dic2,verb=False):
    
    dic1 = cleandic(dict(dic1))
    dic2 = cleandic(dict(dic2))

    r = True
    kset1 = set(dic1.keys())
    kset2 = set(dic2.keys())

    dset = set.difference(kset1,kset2)
    iset = set.intersection(kset1,kset2)

    if len(dset) != 0:
        r = False
        if verb:
            print "Keys not in both sets:",dset

    for k in iset:
        #print k,dic1[k],dic2[k]
        if type(dic1[k]) == np.float32:
            if abs(dic1[k] - dic2[k]) > dtol:
                r = False
                if verb:
                    print "Key mismatch:",k,dic1[k],dic2[k]
        else:
            if dic1[k] != dic2[k]:
                r = False
                if verb:
                    print "Key mismatch:",k,dic1[k],dic2[k]

    #for k in iset:
    #    if dic1[k] != dic2[k]:
    #        r = False
    #        if verb:
    #            print "Key mismatch:",k,dic1[k],dic2[k]
    return r

def cleandic(dic):
  
    if dic.has_key("FDSLICECOUNT"):
        del(dic["FDSLICECOUNT"])
    if dic.has_key("RAWFDATA"):
        del(dic["RAWFDATA"])
    if dic.has_key("CURFILE"):
        del(dic["CURFILE"])
    return dic

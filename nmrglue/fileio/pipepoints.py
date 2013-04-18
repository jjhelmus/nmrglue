"""
Calculates the points in a nmr spectrum, given
a dictionary with data information, and structured loaded sparky list
"""
import numpy as np
import numpy.lib.recfunctions as rfn
def pts(stepX,highX,stepY,highY,HN,N):
    ptsX=(highX-HN)/stepX; # $process[2]=HN-proton
    ptsY=(highY-N)/stepY; # $process[1]=N-nitrogen
    return(ptsX,ptsY)

def calc(dic,slist):
    if not isinstance(slist, np.core.records.recarray): # 
        print "SPARKY list is not given in numpy.core.records.recarray"
        print "Use: nmrglue.fileio.sparkylist.read_HN_N(peakfile) "
        return ()
    frqX=dic['FDF2OBS']; frqY=dic['FDF1OBS']
    sizeX=dic['FDSIZE']; sizeY=dic['FDSPECNUM']
    origX=dic['FDF2ORIG']; origY=dic['FDF1ORIG']
    swX=dic['FDF2SW']; swY=dic['FDF1SW']
    #
    stepX=swX/frqX/sizeX;
    highX=(origX+swX)/frqX;
    stepY=swY/frqY/sizeY;
    highY=(origY+swY)/frqY;
    ptsX,ptsY = pts(stepX,highX,stepY,highY,slist['HN'],slist['N'])
    out = zip(ptsX,ptsY)
    out = np.array(out, dtype=[('ptsX',float),('ptsY',float)])
    out = rfn.merge_arrays((slist,out),flatten=True)
    out = out.view(np.recarray)
    return(out)

#! /usr/bin/env python

import nmrglue as ng

# create the iterator (xiter) which will loop over YX planes
xiter = ng.pipe.iter3D("../common_data/3d_pipe/data/test%03d.fid",'x','x')

print "Processing XY planes..."
for i,(dic,plane) in enumerate(xiter):
    
    # process the direct dimention (x)
    dic,plane = ng.pipe_proc.zf(dic,plane,auto=True)
    dic,plane = ng.pipe_proc.ft(dic,plane,auto=True)
    dic,plane = ng.pipe_proc.ps(dic,plane,p0=0.0,p1=0.0)
    dic,plane = ng.pipe_proc.di(dic,plane)

    # process the first indirect (y) dimention
    dic,plane = ng.pipe_proc.tp(dic,plane)
    dic,plane = ng.pipe_proc.zf(dic,plane,auto=True)
    dic,plane = ng.pipe_proc.ft(dic,plane,auto=True)
    dic,plane = ng.pipe_proc.ps(dic,plane,p0=0.0,p1=0.0)
    dic,plane = ng.pipe_proc.di(dic,plane)
    dic,plane = ng.pipe_proc.tp(dic,plane)
   
    # write out the plane
    #print "Writing out Z slice:",i,"/",xiter.i_max
    xiter.write("./ft/test%03d.ft2",plane,dic)


print "Processing ZX planes..."

# create the iterator (ziter) which will loop over XZ planes
ziter = ng.pipe.iter3D("./ft/test%03d.ft2",'z','z')

for i,(dic,plane) in enumerate(ziter):

    # process the z-dim
    dic,plane = ng.pipe_proc.zf(dic,plane,auto=True)
    dic,plane = ng.pipe_proc.ft(dic,plane,auto=True)
    dic,plane = ng.pipe_proc.ps(dic,plane,p0=-92.0,p1=65.0)
    dic,plane = ng.pipe_proc.di(dic,plane)

    # write out the plane
    #print "Writing out Y slice:",i,"/",ziter.i_max
    ziter.write("./ft/test%03d.ft3",plane,dic)

print "Done"
# check against a file processed with NMRPipe
dic1,data1 = ng.pipe.read_lowmem("../common_data/3d_pipe/ft/test%03d.ft3")
dic2,data2 = ng.pipe.read_lowmem("./ft/test%03d.ft3")
# we'll check two 2D slices as the whole 3D takes a long time
print ng.misc.pair_similar(dic1,data1[5],dic2,data2[5],verb=True)
print ng.misc.pair_similar(dic1,data1[52],dic2,data2[52],verb=True)


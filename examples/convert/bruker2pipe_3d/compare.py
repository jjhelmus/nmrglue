#! /usr/bin/env python

import nmrglue as ng

# check the Python conversion against the conversion using NMRPipe's var2pipe
pdic, pdata = ng.pipe.read("./data/3d_pipe%03d.fid")
pdic2, pdata2 = ng.pipe.read("./bruker_3d/fid/test%03d.fid")
r1, r2 = ng.misc.pair_similar(pdic, pdata, pdic2, pdata2, verb=True)

print "==================="
print "Summary:"

if r1:
    print "Data arrays are similar."
else:
    print "Data arrays differ as listed above."

if r2:
    print "Spectral parameters are similar."
else:
    print "Spectral parameters differ as listed above."

#! /usr/bin/env python

import nmrglue as ng

# check the Python conversion against the conversion using NMRPipe's var2pipe
sdic,sdata = ng.sparky.read("sparky_2d.ucsf")
sdic2,sdata2 = ng.sparky.read("nmrpipe_2d/sparky_2d.ucsf")
r1, r2 = ng.misc.pair_similar(sdic, sdata, sdic2, sdata2, verb=True)
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

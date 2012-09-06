#! /usr/bin/env python

import nmrglue as ng

# check the Python conversion against the conversion using NMRPipe's var2pipe
pdic, pdata = ng.pipe.read('1d_pipe.fid')
pdic2, pdata2 = ng.pipe.read('agilent_1d/test.fid')
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

#!/sbinlab2/software/python-enthought-dis/epd-7.3-2-rh5-x86_64/bin/python2.7
import sys, os
import nmrglue as ng

if len(sys.argv) < 2:
    sys.exit('Usage: %s sparky.list' % sys.argv[0])

peakfile = ng.fileio.sparkylist.read_HN_N(sys.argv[1])
max_resi = max(peakfile['resi'])
SEQ = ["-" for x in range(max_resi)]
for peak in peakfile:
    resn = peak['resn']
    oneAA = peak['oneAA']
    AA = peak['AA']
    resi = peak['resi']
    SEQ[resi-1]="%s"%oneAA
print ">%s|PDBID|CHAIN|SEQUENCE"%sys.argv[1]
print ''.join(SEQ)

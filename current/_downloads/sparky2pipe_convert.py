import nmrglue as ng

# read in the Sparky file
sdic, sdata = ng.sparky.read('data.ucsf')

# convert to NMRPipe format
C = ng.convert.converter()
C.from_sparky(sdic, sdata)
pdic, pdata = C.to_pipe()

# write results to NMRPipe file
ng.pipe.write('data.ft2', pdic, pdata)

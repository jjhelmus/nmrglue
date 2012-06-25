#!/bin/csh

# The first three have FDMIN/FDMAX differences
#nmrPipe -in ./1D_freq_real.dat  \
#| nmrPipe  -fn FSH  -ls 1 -sw      \
#-ov -out fsh1.dat

#nmrPipe -in ./1D_freq_real.dat  \
#| nmrPipe  -fn FSH  -rs 8 -sw      \
#-ov -out fsh2.dat

#nmrPipe -in ./1D_freq_real.dat  \
#| nmrPipe  -fn FSH  -rs 6.7 -sw    \
#-ov -out fsh3.dat

# NMRPipe performs a Hilbert transform here, pipe_proc does not
#nmrPipe -in ./1D_freq_complex.dat  \
#| nmrPipe  -fn FSH  -ls 9.5  -sw   \
#-ov -out fsh4.dat

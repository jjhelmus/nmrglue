#!/bin/csh

# process 2D complex data with Gaussian and cos^2 windows
nmrPipe -in time_complex.fid                   \
|   nmrPipe -fn GMB -gb 0.1 -lb -8 -c 0.5             \
|   nmrPipe -fn ZF -auto                \
|   nmrPipe -fn FT -auto                     \
|   nmrPipe -fn PS -p0 -36 -p1 0 -di       \
|   nmrPipe -fn TP                           \
|   nmrPipe -fn SP -off 0.5 -pow 2.0 -c 0.5  \
|   nmrPipe -fn ZF -auto            \
|   nmrPipe -fn FT -auto                     \
|   nmrPipe -fn PS -p0 -7 -p1 0 -di       \
-ov -out 2d_complex_processing1.dat

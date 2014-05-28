#!/bin/csh

# process 2D mixed mode data
nmrPipe -in time_real.fid                   \
|   nmrPipe -fn GMB -gb 0.1 -lb -8 -c 0.5             \
|   nmrPipe -fn ZF -auto                \
|   nmrPipe -fn FT -alt                 \
|   nmrPipe -fn PS -p0 0 -p1 0       \
|   nmrPipe -fn TP -hyper                    \
|   nmrPipe -fn SP -off 0.5 -pow 2.0 -c 0.5  \
|   nmrPipe -fn ZF -auto            \
|   nmrPipe -fn FT -auto                     \
|   nmrPipe -fn PS -p0 0 -p1 0 -di     \
-ov -out 2d_mixed_processing1.dat

nmrPipe -in time_real.fid                   \
|   nmrPipe -fn EM -lb 8           \
|   nmrPipe -fn ZF -auto                \
|   nmrPipe -fn FT -auto -di                    \
|   nmrPipe -fn TP                       \
|   nmrPipe -fn SP -off 0.5 -pow 1.0 -c 0.5  \
|   nmrPipe -fn ZF -auto            \
|   nmrPipe -fn FT -auto                     \
|   nmrPipe -fn MC    \
-ov -out 2d_mixed_processing2.dat

#!/bin/csh

nmrPipe -in ./1D_time_real.fid \
| nmrPipe -fn HT -ps90-180 \
-ov -out ht41.dat

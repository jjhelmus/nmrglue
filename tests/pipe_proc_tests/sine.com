#! /bin/sh

nmrPipe -in ./time_complex.fid                              \
| nmrPipe -fn SINE -off 0.35 -end 0.98 -pow 2 -c 1.0  \
-ov -out sine1.dat

nmrPipe -in ./time_complex.fid                              \
| nmrPipe -fn SINE -off 0.10 -end 0.75 -pow 1 -c 0.5  -size 200 -one \
-ov -out sine2.dat

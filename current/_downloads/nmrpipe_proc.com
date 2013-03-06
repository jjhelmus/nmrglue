#!/bin/csh

bruk2pipe -in ./expnmr_00001_1/fid \
  -bad 0.0 -noaswap -DMX -decim 32 -dspfvs 12 -grpdly 0  \
  -xN             32768  \
  -xT             16384  \
  -xMODE            DQD  \
  -xSW         4807.692  \
  -xOBS         400.132  \
  -xCAR           4.697  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./test.fid -verb -ov

nmrPipe -in test.fid \
| nmrPipe -fn ZF -auto                      \
| nmrPipe -fn FT                            \
| nmrPipe -fn PS -p0 -22.0 -p1 0.0 -di      \
   -out test.ft2 -verb -ov

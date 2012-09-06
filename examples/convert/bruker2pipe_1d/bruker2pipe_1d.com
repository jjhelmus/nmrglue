#!/bin/csh

bruk2pipe -in ./bruker_1d/fid \
  -bad 0.0 -noaswap -AMX -decim 16 -dspfvs 12 -grpdly 0  \
  -xN              4096  \
  -xT              2048  \
  -xMODE            DQD  \
  -xSW        10000.000  \
  -xOBS         600.133  \
  -xCAR           4.773  \
  -xLAB              1H  \
  -ndim               1  \
  -out ./bruker_1d/test.fid -verb -ov

#!/bin/csh

var2pipe -in ./agilent_1d/fid -noaswap \
  -xN              3000  \
  -xT              1500  \
  -xMODE        Complex  \
  -xSW        50000.000  \
  -xOBS         125.681  \
  -xCAR          99.000  \
  -xLAB             C13  \
  -ndim               1  \
  -out ./agilent_1d/test.fid -verb -ov

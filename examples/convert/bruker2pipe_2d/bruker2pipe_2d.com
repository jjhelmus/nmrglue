#!/bin/csh

bruk2pipe -in ./bruker_2d/ser \
  -bad 0.0 -noaswap -AMX -skip -1 -decim 16 -dspfvs 12 \
  -xN              1536  -yN               600  \
  -xT               650  -yT               300  \
  -xMODE        Complex  -yMODE        Complex  \
  -xSW        11061.947  -ySW         4000.000  \
  -xOBS         800.134  -yOBS         201.204  \
  -xCAR           4.773  -yCAR          58.742  \
  -xLAB              1H  -yLAB             13C  \
  -ndim               2  -aq2D          States  \
  -out ./bruker_2d/test.fid -verb -ov


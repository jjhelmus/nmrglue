#!/bin/csh

var2pipe -in ./agilent_2d/fid -noaswap  \
  -xN              3000  -yN               332  \
  -xT              1500  -yT               166  \
  -xMODE        Complex  -yMODE        Complex  \
  -xSW        50000.000  -ySW         5555.556  \
  -xOBS         125.691  -yOBS          50.648  \
  -xCAR          55.000  -yCAR         120.000  \
  -xLAB             13C  -yLAB             15N  \
  -ndim               2  -aq2D          States  \
  -out ./agilent_2d/test.fid -verb -ov

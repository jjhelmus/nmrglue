#!/bin/csh

var2pipe -in ./agilent_2d_tppi/fid \
 -noaswap  \
  -xN              2800  -yN               600  \
  -xT              1400  -yT               600  \
  -xMODE        Complex  -yMODE           TPPI  \
  -xSW        50000.000  -ySW        33333.333  \
  -xOBS         125.681  -yOBS         125.681  \
  -xCAR         101.274  -yCAR         101.274  \
  -xLAB            C13x  -yLAB            C13y  \
  -ndim               2  -aq2D            TPPI  \
  -out ./agilent_2d_tppi/test.fid -verb -ov

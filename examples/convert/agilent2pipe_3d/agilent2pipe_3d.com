#!/bin/csh

var2pipe -in ./agilent_3d/fid -noaswap -aqORD 1 \
  -xN              2500  -yN                88  -zN               128  \
  -xT              1250  -yT                44  -zT                64  \
  -xMODE        Complex  -yMODE         States  -zMODE         States  \
  -xSW        50000.000  -ySW         2777.778  -zSW         5555.556  \
  -xOBS         125.676  -yOBS          50.648  -zOBS         125.676  \
  -xCAR          56.000  -yCAR         120.000  -zCAR          56.000  \
  -xLAB              CX  -yLAB               N  -zLAB              CA  \
  -ndim               3  -aq2D          States                         \
  -out ./agilent_3d/data/test%03d.fid -verb -ov

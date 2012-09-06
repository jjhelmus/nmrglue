#!/bin/csh

bruk2pipe -in ./bruker_3d/ser  \
  -bad 0.0 -noaswap -AMX -decim 16 -dspfvs 12 -grpdly 0  \
  -xN              1536  -yN               128  -zN               116  \
  -xT               650  -yT                64  -zT                58  \
  -xMODE            DQD  -yMODE         States  -zMODE         States  \
  -xSW        11061.947  -ySW         2500.000  -zSW         5555.556  \
  -xOBS         800.134  -yOBS          81.086  -zOBS         201.204  \
  -xCAR           4.784  -yCAR         119.787  -zCAR          55.743  \
  -xLAB              1H  -yLAB             15N  -zLAB             13C  \
  -ndim               3  -aq2D          States                         \
  -out ./bruker_3d/fid/test%03d.fid -verb -ov

#!/bin/csh

xyz2pipe -in ./nmrpipe_3d/ft/test%03d.ft3 -x > ./nmrpipe_3d/full.ft3
pipe2ucsf ./nmrpipe_3d/full.ft3 ./nmrpipe_3d/sparky_3d.ucsf
rm ./nmrpipe_3d/full.ft3

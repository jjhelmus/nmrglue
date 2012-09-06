#!/bin/csh

pipe2ucsf ./nmrpipe_2d/test.ft2 ./nmrpipe_2d/sparky_2d.ucsf

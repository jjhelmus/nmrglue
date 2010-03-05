#! /usr/bin/env python
# Moves files listed in files_to_move.txt to _build directory

import shutil
import os

f = open("files_to_move.txt")
files = [x.strip() for x in f]
f.close()

for file in files:
    src = file
    dst = file.replace("./examples/el/","_build/html/examples/el/")

    # make directories as needed
    d,f = os.path.split(dst)
    if d !='' and os.path.exists(d) == False:
        os.makedirs(d)

    # copy the file
    shutil.copy(src,dst)

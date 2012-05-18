# /usr/bin/env python

import glob

fnames = glob.glob('generated/*.rst')

def fix_function_file(fname):
    
    # parse the filename for module and function names
    good_mod_name = ".".join(fname.split('/')[1].split('.')[0:3])
    bad_mod_name = ".".join(good_mod_name.split('.')[0:2])

    good_func_name = fname.split('.')[-2]
    bad_func_name = '.'.join(fname.split('.')[-3:-1])

    # read in the bad text
    f = open(fname, 'r')
    text = f.read()
    f.close()
   
    if 'autofunction::' not in text:
        # fix the text
        text = text.replace(bad_mod_name, good_mod_name)
        text = text.replace('automethod::', 'autofunction::')
        text = text.replace(bad_func_name, good_func_name)

        # write out the good text
        f = open(fname, 'w')
        f.write(text)
        f.close()

    return

classes = ['bruker_nd']

for fname in fnames:
    name = fname.split('.')[-2]
    if name in classes:
        print "Fixing class file:", fname
        pass
    else:
        print "Fixing function file:", fname
        fix_function_file(fname)


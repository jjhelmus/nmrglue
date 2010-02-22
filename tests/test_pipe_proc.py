#! /usr/bin/env python

# Test pipe_proc functions against nmrPipe.
# This script tests functions in the pipe_proc module against
# NMRPipe scripts to see if they produce similar files.

# The usage of the script is two-fold.  Calling it with no arguments
# runs all standard pipe_proc tests in the pipe_proc_tests folder
# If the script is called with one or more arguments, then the tests listed
# as arguments are executed.

import os
import sys
import nmrglue.fileio.pipe as pipe
import nmrglue.util.misc as util

verb1=True   # basic verbosity
verb2=True   # verbosity of util functions

test_dir = './pipe_proc_tests'
execfile(os.path.join(test_dir,"all_tests.py"))

# load suite with the desired test(s)
if len(sys.argv) > 1:
    suite = dict()
    for k in sys.argv[1:]:
        suite[k] = eval(k)

test_passed = 0
test_failed = 0
list_fail = []
print "Beginning pipe_proc test suite"

os.chdir(test_dir)

for key in suite:   # loop over requested tests

    k_test = True
    if verb1:
        print "---------------------------------"

    d = suite[key]
    # run the nmrPipe program
    os.system(d["sh"])
    # run the python program
    os.system(d["py"])
    print key,"...",
    
    if verb1:
        print ""

    # check each file for match
    for f1,f2 in zip(d["f"][::2],d["f"][1::2]):
        
        if verb1:
            print "checking",f1,f2
        dic1,data1 = pipe.read(f1)
        dic2,data2 = pipe.read(f2)
        r1,r2 = util.pair_similar(dic1,data1,dic2,data2,verb2)
        if verb1:
            if r1:
                print "Data: Pass",
            else:
                print "Data: Fail",

            if r2:
                print "Dictionary: Pass"
            else: 
                print "Dictionary: Fail"

        k_test = k_test and r1 and r2
        # remove the files
        os.remove(f1)
        os.remove(f2)
    if k_test:
        print "Passed."
        test_passed = test_passed + 1
    else:
        print "Failed!"
        test_failed = test_failed + 1
        list_fail.append(key)

print "----------------------------------"
print "Total tests:",test_passed+test_failed
print "Tests passed:",test_passed
if test_failed != 0:
    print "Tests failed:",test_failed
    print list_fail

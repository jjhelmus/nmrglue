"""
xcpy - interfaces with external programs from within Bruker-TopSpin

SYNOPSIS
    \t\txcpy
    \t\txcpy [OPTIONS] 
    \t\txcpy [OPTIONS] [SCRIPTNAME]

Description

"""
import sys
import os
from ConfigParser import SafeConfigParser
from subprocess import Popen, PIPE, STDOUT


def check_jython():
    """
    Checks whether Jython is being used to run this script
    from within Topspin

    """
    # some of the functions defined in the global namespace
    # by Topspin which are also used by xcpy
    topspin_inbuilts = ['MSG', 'INPUT_DIALOG', 'CURDATA']

    g = globals().keys()
    for function in topspin_inbuilts:
        if function in g:
            pass
        else:
            raise Exception('This file is meant to be executed \
                             using Jython from within Topspin')
    return True


def topspin_location():
    """
    Gets Topspin home directory. Also serves to check
    whether the script is being executed from within
    Topspin or externally.

    """
    try:
        # topspin >= 3.6
        toppath = sys.getBaseProperties()["XWINNMRHOME"]
    except:
        try:
            # topspin >= 3.1 and <= 3.5
            toppath = sys.registry["XWINNMRHOME"]
        except:
            # topspin < 3.1
            toppath = sys.getEnviron()["XWINNMRHOME"]
    
    # if all the above fail, there should be an exception raised and
    # the function should not return anything
    return toppath


def read_cfg(filename):
    """
    Reads in the configuration file

    """
    config = SafeConfigParser()
    try:
        config.read(filename)
        cpyname = config.get('xcpy', 'cpython')
        scripts_location = config.get('xcpy', 'scripts_location')
    except:
        cpyname, scripts_location = "", ""

    return cpyname, scripts_location


def write_cfg(outfile, infile):
    """
    Writes or overwrites a configuration file

    """
    if infile is not None:
        cpyname, scripts_location = read_cfg(infile)
    else:
        cpyname, scripts_location = "", ""

    cpyname, scripts_location = INPUT_DIALOG(
            "XCPy Configuration", 
            "Please Click on OK to write this configuration.", 
            ["CPython Executable", "CPython Scripts Location"], 
            [cpyname, scripts_location], 
            ["",""], 
            ["1", "1"])

    if not cpyname or not scripts_location: 
        MSG("Invalid configartion specified. Config file not written")
    else:
        config = SafeConfigParser()
        config.add_section('xcpy')
        config.set('xcpy', 'cpython', cpyname)
        config.set('xcpy', 'scripts_location', scripts_location)

        with open(outfile, "w") as f:
            config.write(f)
        MSG('Written Configuration file at: ' + outfile)



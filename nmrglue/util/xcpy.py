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


def exists(filename, raise_error=False):
    """
    Checks whether a file exists either returns False or
    raises an exception

    """
    if os.path.exists(filename):
        return True
    elif raise_error:
        raise Exception('{} not found'.format(filename))
    else:
        return False


def current_data():
    """
    Returns the current EXPNO and PROCNO open in Topspin,
    if executed when a data folder is open

    """
    cd = CURDATA()
    current_dir = os.path.join(cd[3], cd[0])
    current_expno = cd[1]
    current_procno = cd[2]

    return [current_dir, current_expno, current_procno]


def run(cpython, script, pass_current_folder=True, use_shell=None, dry=None):
    """
    Runs a cpython script

    """
    if pass_current_folder == True:
        cd = current_data()
    else:
        cd = []

    if use_shell is None:
        if os.name == 'nt':
            use_shell = True
        else:
            use_shell = False

    args = [cpython, script] + cd

    if dry:
        MSG('The following command will be executed: \n' + ' '.join(args))
        process = None

    else:
        process = Popen(args, stdin=PIPE, stderr=STDOUT, shell=use_shell)
        process.stdin.close()
    
    return process

    
def verify_completion(process):
    """
    Verify that the output is correct
    """

    if process is not None:
        errmsg = [line for line in iter(process.stdout.readline, '')]

        if not errmsg:
            MSG('Script executed to completion')
        else:
            MSG(''.join(errmsg))

    else:
        return
        
    
if __name__ == "__main__":
    
    # check whether running in Topspin/Jython and get the location if yes
    if check_jython():
        toppath = topspin_location()

    # this is where the config file must be store
    config_file = os.path.join(toppath, 'exp', 'stan', 'nmr', 'py', 'user', 'xcpy.cfg')

    # get arguments passed to xcpy 
    argv = sys.argv
    if len(argv) == 1:
        argv.append('--info')

    # return docstring
    if argv[1] in ['-i', '--info']: 
        if len(argv) > 2:
            MSG('Opening Documentation. All other options ignored. Press OK')
        MSG(__doc__)

    # check if configuration exists
    elif not os.path.exists(config_file):
        MSG("Configuration file does not exists. Will open an input dialog to write it")
        write_cfg(config_file)

    # if configuration settings are to be changed
    elif argv[1] in ['-s', '--settings']:
        if len(argv) > 2:
            MSG('Opening Configuration Settings. All other options ignored. Press OK')
        write_cfg(config_file, config_file)

    else: 
        if '--use-shell' in argv:
            use_shell = True
        else:
            use_shell = False

        if '--no-args' in argv:
            pass_current_folder = False
        else:
            pass_current_folder = True

        if '--dry-run' in argv:
            dry = True
        else:
            dry = False

        # read configuration
        cpyname, folder = read_cfg(config_file)

        # see if script is there and then run
        if exists(cpyname) and exists(folder):
            scriptname = os.path.join(folder, argv[1])

            if exists(scriptname):
                process = run(cpyname, scriptname, pass_current_folder, use_shell, dry)
                verify_completion(process)
            else:
                scriptname = scriptname + '.py'
                if exists(scriptname, raise_error=True):
                    process = run(cpyname, scriptname, pass_current_folder, use_shell, dry)
                    verify_completion(process)

#! /usr/bin/env python
# pipe_proc module tests, calls scripts in pipe_proc_tests directory.

from __future__ import print_function

import os
import nmrglue as ng

# Framework functions


def _perform_test(glue_script, pipe_script, glue_files, pipe_files,
                  ignore_pipe_display=False):
    """
    """
    cwd_backup = os.getcwd()    # save the current working directory

    # decent into the pipe_proc_tests directory
    script_dir, script_fname = os.path.split(os.path.realpath(__file__))
    os.chdir(os.path.join(script_dir, 'pipe_proc_tests'))

    # execute the scripts and compare the files
    exec(open(glue_script).read())
    os.system(pipe_script)

    for glue_file, pipe_file in zip(glue_files, pipe_files):

        glue_dic, glue_data = ng.pipe.read(glue_file)
        pipe_dic, pipe_data = ng.pipe.read(pipe_file)
        os.remove(glue_file)
        os.remove(pipe_file)
        r1, r2 = ng.misc.pair_similar(glue_dic, glue_data, pipe_dic,
                                      pipe_data, True,
                                      ignore_pipe_display=ignore_pipe_display)
        print(glue_file, pipe_file, r1, r2)
        assert r1 is True
        assert r2 is True

    # return to the backup-ed directory.
    os.chdir(cwd_backup)
    return


def _standard_args(func_name, num_files):
    """ Generate a set of standard  """
    glue_files = [func_name+str(i + 1)+'.glue' for i in range(num_files)]
    pipe_files = [func_name+str(i + 1)+'.dat' for i in range(num_files)]
    pipe_script = './' + func_name + '.com'
    glue_script = './' + func_name + '.py'
    return glue_script, pipe_script, glue_files, pipe_files


def _standard_test(func_name, num_files, ignore_pipe_display=False):
    return _perform_test(*_standard_args(func_name, num_files),
                         ignore_pipe_display=ignore_pipe_display)


#########################
# Apodization functions #
#########################


def test_apod():
    """ APOD function """
    return _standard_test('apod', 7)


def test_em():
    """ EM function """
    return _standard_test('em', 2)


def test_gm():
    """ GM function """
    return _standard_test('gm', 3)


def test_gmb():
    """ GMB function """
    return _standard_test('gmb', 2)


def test_jmod():
    """ JMOD function """
    return _standard_test('jmod', 2)


def test_sp():
    """ SP function """
    return _standard_test('sp', 2)


def test_sine():
    """ SINE function """
    return _standard_test('sine', 2)


def test_tm():
    """ TM function """
    return _standard_test('tm', 2)


def test_tri():
    """ TRI function """
    return _standard_test('tri', 2)

###################
# Shift functions #
###################


def test_rs():
    """ RS function """
    return _standard_test('rs', 5)


def test_ls():
    """ LS function """
    return _standard_test('ls', 5)


def test_cs():
    """ CS function """
    return _standard_test('cs', 8)


# XXX fsh test 1-4 fail
def test_fsh():
    """ FSH function """
    return _standard_test('fsh', 0)


##############
# Transforms #
##############


def test_ft():
    """ FT function """
    return _standard_test('ft', 8)


def test_rft():
    """ RFT function """
    #return _standard_test('rft', 14)    # XXX tests 9-11 fail
    glue_files = ['rft1.glue', 'rft2.glue', 'rft3.glue', 'rft4.glue',
                  'rft5.glue', 'rft6.glue', 'rft7.glue', 'rft8.glue',
                  'rft12.glue', 'rft13.glue', 'rft14.glue']
    pipe_files = ['rft1.dat', 'rft2.dat', 'rft3.dat', 'rft4.dat',
                  'rft5.dat', 'rft6.dat', 'rft7.dat', 'rft8.dat',
                  'rft12.dat', 'rft13.dat', 'rft14.dat']
    pipe_script = './rft.com'
    glue_script = './rft.py'
    return _perform_test(glue_script, pipe_script, glue_files, pipe_files)


def test_ha():
    """ HA function """
    return _standard_test('ha', 2)


def test_ht():
    """ HT function """
    #return _standard_test('ht', 8)  # XXX test 4 fails
    glue_files = ['ht1.glue', 'ht2.glue', 'ht3.glue', 'ht5.glue', 'ht6.glue',
                  'ht7.glue', 'ht8.glue']
    pipe_files = ['ht1.dat', 'ht2.dat', 'ht3.dat', 'ht5.dat', 'ht6.dat',
                  'ht7.dat', 'ht8.dat']
    pipe_script = './ht.com'
    glue_script = './ht.py'
    return _perform_test(glue_script, pipe_script, glue_files, pipe_files)


##########################
# Standard NMR Functions #
##########################


def test_ps():
    """ PS function """
    return _standard_test('ps', 6)


def test_tp():
    """ TP function """
    return _standard_test('tp', 9, ignore_pipe_display=True)


def test_ytp():
    """ YTP function """
    return _standard_test('ytp', 3)


def test_xy2yx():
    """ XY2YX function """
    return _standard_test('xy2yx', 3)


def test_zf():
    """ ZF function """
    return _standard_test('zf', 5)


###################
# Basic Utilities #
###################


def test_add():
    """ ADD function """
    return _standard_test('add', 4)     # Note that test 5 fails intensionally


def test_dx():
    """ DX function """
    return _standard_test('dx', 1)


def test_ext():
    """ EXT function """
    return _standard_test('ext', 11, ignore_pipe_display=True)


def test_integ():
    """ INTEG function """
    return _standard_test('integ', 1)


def test_mc():
    """ MC function """
    return _standard_test('mc', 2)


def test_mir():
    """ MIR function """
    return _standard_test('mir', 14)


def test_mult():
    """ MULT function """
    return _standard_test('mult', 3)    # Note that test 4 fails intensionally


def test_rev():
    """ REV function """
    return _standard_test('rev', 3)


def test_set():
    """ SET function """
    return _standard_test('set', 4)


def test_shuf():
    """ SHUF function """
    return _standard_test('shuf', 7)


def test_sign():
    """ SIGN function """
    return _standard_test('sign', 8)


########
# Misc #
########


def test_coadd():
    """ COADD function """
    return _standard_test('coadd', 2)


def test_coad():
    """ COAD function """
    return _standard_test('coad', 2)


def test_dev():
    """ DEV function """
    return _standard_test('dev', 0)


def test_null():
    """ NULL function """
    return _standard_test('null', 2)


def test_qart():
    """ QART function """
    return _standard_test('qart', 2)


def test_qmix():
    """ QMIX function """
    return _standard_test('qmix', 2)


def test_smo():
    """ SMO function """
    return _standard_test('smo', 3)


def test_zd():
    """ ZD function """
    return _standard_test('zd', 4)


def test_save():
    """ SAVE function """
    return _standard_test('save', 2)


######################
# Baseline functions #
######################


def test_base():
    """ BASE function """
    return _standard_test('base', 7)


def test_cbf():
    """ CBF function """
    return _standard_test('cbf', 4)

#########
# Units #
#########


def test_units():
    """ Units """
    return _standard_test('units', 17)

#####################
# Inetgration Tests #
#####################


def test_2D_complex_processing():
    """ 2D complex mode processing pipeline """
    return _standard_test('2d_complex_processing', 1)


####################################################
# Test which are known to fail for various reasons.#
####################################################

# Tests which fail because of differences between NMRPipe and nmrglue
"""
def test_known_fail():
    pipe_files = ['add5.dat',
                  'mult4.dat',
                  'shuf8.dat', 'shuf9.dat', 'shuf10.dat']
    glue_files = ['add5.glue',
                  'mult4.glue',
                  'shuf8.glue', 'shuf9.glue', 'shuf10.glue']
    pipe_script = 'known_fail.com'
    glue_script = 'known_fail.py'
    # also 'shuf11', 'shuf12' and dev1' test all fail intensionally.
    return _perform_test(glue_script, pipe_script, glue_files, pipe_files)
"""

# Test which fail and should be fixed
"""
def test_to_fix():
    pipe_files = ["fsh1.dat", "fsh2.dat", "fsh3.dat", "fsh4.dat"
                  "rft9.dat", "rft10.dat", "rft11.dat",
                  "ht4.dat"]
    glue_files = ["fsh1.glue", "fsh2.glue", "fsh3.glue", "fsh4.glue"
                  "rft9.glue", "rft10.glue", "rft11.glue",
                  "ht4.glue"]
    pipe_script = 'to_fix.com'
    glue_script = 'to_fix.py'
    return _perform_test(glue_script, pipe_script, glue_files, pipe_files)
"""

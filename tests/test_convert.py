""" Tests for the fileio.convert submodule """

import tempfile
import os
import shutil
import glob

import nmrglue as ng
from numpy.testing import assert_array_equal

from setup import DATA_DIR

# subroutines
def check_sdic(dic1, dic2, exclude=None, v=False):
    """ Check two Sparky parameter dictionaries against each other """ 
    if exclude is None:
        exclude = []
    # check axis dictionaries
    if 'w1' in dic1:
        check_pdic(dic1["w1"], dic2["w1"], exclude=exclude, v=v)
        exclude.append('w1')

    if 'w2' in dic1:
        check_pdic(dic1["w2"], dic2["w2"], exclude=exclude, v=v)
        exclude.append('w2')

    if 'w3' in dic1:
        check_pdic(dic1["w3"], dic2["w3"], exclude=exclude, v=v)
        exclude.append('w3')

    # check remainder
    check_pdic(dic1, dic2, exclude=exclude, v=v)

def check_dic(dic1, dic2, exclude=None, v=False):
    """ Compare two parameter dictionaries """
    if exclude is None:
        exclude = []
    e1 = [k for k in dic1.keys() if k not in dic2.keys() and k not in exclude]
    e2 = [k for k in dic2.keys() if k not in dic2.keys() and k not in exclude]

    if v:
        print e1
    assert len(e1) == 0
    if v:
        print e2 
    assert len(e2) == 0

    if v:
        for k in dic1.keys():
            if k in exclude:
                continue
            if dic1[k] != dic2[k]:
                print k, dic1[k], dic2[k]

    for k in dic1.keys():
        if k in exclude:
            continue
        assert dic1[k] == dic2[k]

def check_pdic(dic1, dic2, exclude=None, v=False):
    """ Compare two NMRPipe parameter dictionaries """
    if exclude is None:
        exclude = []
    e1 = [k for k in dic1.keys() if k not in dic2.keys() and k not in exclude]
    e2 = [k for k in dic2.keys() if k not in dic2.keys() and k not in exclude]

    if v:
        print e1
    assert len(e1) == 0
    if v:
        print e2 
    assert len(e2) == 0

    if v:
        for k in dic1.keys():
            if k in exclude:
                continue
            if type(dic1[k]) == str or type(dic1[k]) == list: 
                if dic1[k] != dic2[k]:
                    print k, dic1[k], dic2[k]
            elif abs(dic1[k] - dic2[k]) >= 0.002:
                print k, dic1[k], dic2[k], dic1[k] - dic2[k]

    for k in dic1.keys():
        if k in exclude:
            continue
        if type(dic1[k]) == str or type(dic1[k]) == list:
            assert dic1[k] == dic2[k]
        else:
            assert abs(dic1[k] - dic2[k]) <= 0.002



bad_varian_keys = ["procpar"]
bad_pipe_keys = ["FDYEAR", "FDMONTH", "FDDAY", "FDHOURS", "FDMINS", "FDSECS"]
bad_bruker_keys = ["pprog", "acqus", "acqu2s"]
bad_sparky_keys = ['bsize', 'extended', 'date', 'owner']

# tests

def test_agilent_1d():
    """ 1D time agilent, pipe <-> agilent, pipe """
    # prepare agilent converter
    vdic, vdata = ng.varian.read(DATA_DIR + "agilent_1d")
    uvdic = ng.varian.guess_udic(vdic, vdata)
    vC = ng.convert.converter()
    vC.from_varian(vdic, vdata, uvdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "agilent_1d/test.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # agilent -> agilent
    cdic, cdata = vC.to_varian()
    assert_array_equal(vdata, cdata)
    check_dic(vdic, cdic, bad_varian_keys)
    
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write(td, cdic, cdata)
    rdic, rdata = ng.varian.read(td)
    assert_array_equal(vdata, rdata)
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

    # agilent -> pipe
    cdic, cdata = vC.to_pipe()
    assert_array_equal(pdata, cdata)
    #check_pdic(pdic, cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    #check_pdic(pdic, cdic)   # XXX don't check dictionary
    os.remove(tf)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata)
    check_pdic(pdic, cdic, bad_pipe_keys)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    check_pdic(pdic, cdic, bad_pipe_keys)
    os.remove(tf)

    # pipe -> agilent
    cdic, cdata = pC.to_varian()
    assert_array_equal(vdata, cdata)
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write(td, cdic, cdata)
    rdic, rdata = ng.varian.read(td)
    assert_array_equal(vdata, rdata)
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

def test_agilent_2d():
    """ 2D time agilent, pipe <-> agilent, pipe """
    # prepare Varian converter
    vdic, vdata = ng.varian.read(DATA_DIR + "agilent_2d/")
    uvdic = ng.varian.guess_udic(vdic, vdata)
    vC = ng.convert.converter()
    vC.from_varian(vdic, vdata, uvdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "agilent_2d/test.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # varian -> varian
    cdic, cdata = vC.to_varian()
    assert_array_equal(vdata, cdata)
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write(td, cdic, cdata)
    rdic, rdata = ng.varian.read(td)
    assert_array_equal(vdata, rdata)
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

    # varian -> pipe
    cdic, cdata = vC.to_pipe()
    assert_array_equal(pdata, cdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    os.remove(tf)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata)
    check_pdic(pdic, cdic, bad_pipe_keys)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    check_pdic(pdic, cdic, bad_pipe_keys)
    os.remove(tf)

    # pipe -> varian
    cdic, cdata = pC.to_varian()
    assert_array_equal(vdata, cdata)
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write(td, cdic, cdata)
    rdic, rdata = ng.varian.read(td)
    assert_array_equal(vdata, rdata)
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

def test_agilent_3d():
    """ 3D time agilent, pipe <-> agilent, pipe """
    # prepare Agilent converter
    vdic, vdata = ng.varian.read(DATA_DIR + "agilent_3d")
    uvdic = ng.varian.guess_udic(vdic, vdata)
    vC = ng.convert.converter()
    vC.from_varian(vdic, vdata, uvdic)

    # prepare NMRPipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "agilent_3d/data/test%03d.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # agilent -> agilent
    cdic, cdata = vC.to_varian()
    assert_array_equal(vdata, cdata)
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write(td, cdic, cdata)
    rdic, rdata = ng.varian.read(td)
    assert_array_equal(vdata, rdata)
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

    # agilent -> pipe
    cdic, cdata = vC.to_pipe()
    assert_array_equal(pdata, cdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata)
    bpk = list(bad_pipe_keys)
    bpk.append("FDDISPMAX")     # nmrglue doesn't update the MIN/MAX values
    bpk.append("FDMIN")
    bpk.append("FDDISPMIN")
    bpk.append("FDSCALEFLAG")
    bpk.append("FDMAX")
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    check_pdic(pdic, cdic, bpk)
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> agilent
    cdic, cdata = pC.to_varian()
    assert_array_equal(vdata, cdata)
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write(td, cdic, cdata)
    rdic, rdata = ng.varian.read(td)
    assert_array_equal(vdata, rdata)
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

def test_bruker_1d():
    """ 1D time bruker, pipe <-> bruker, pipe """
    # prepare Bruker converter
    bdic, bdata = ng.bruker.read(DATA_DIR + "bruker_1d")
    ubdic = ng.bruker.guess_udic(bdic, bdata)
    bC = ng.convert.converter()
    bC.from_bruker(bdic, bdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "bruker_1d/test.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # bruker -> bruker
    cdic, cdata = bC.to_bruker()
    assert_array_equal(bdata, cdata)
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write(td, cdic, cdata)
    rdic, rdata = ng.bruker.read(td)
    assert_array_equal(bdata, rdata)
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    shutil.rmtree(td) 

    # bruker -> pipe
    cdic, cdata = bC.to_pipe()
    assert_array_equal(pdata, cdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    os.remove(tf)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata)
    #bpdk = list(bad_pipe_keys)
    #bpdk.append("FDF2ORIG")         # these are roundoff errors
    #bpdk.append("FDFLTORDER")       # not a problem when written
    check_pdic(pdic, cdic, bad_pipe_keys)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    check_pdic(pdic, cdic, bad_pipe_keys)
    os.remove(tf)

    # pipe -> bruker
    cdic, cdata = pC.to_bruker()
    assert_array_equal(bdata, cdata)
    check_dic(bdic, cdic, bad_bruker_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write(td, cdic, cdata)
    rdic, rdata = ng.bruker.read(td)
    assert_array_equal(bdata, rdata)
    check_dic(bdic, cdic, bad_bruker_keys)
    shutil.rmtree(td)

def test_bruker_2d():
    """ 2D time bruker, pipe <-> bruker, pipe """
    # prepare Bruker converter
    bdic, bdata = ng.bruker.read(DATA_DIR + "bruker_2d")
    ubdic = ng.bruker.guess_udic(bdic, bdata)
    bC = ng.convert.converter()
    bC.from_bruker(bdic, bdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "bruker_2d/test.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # bruker -> bruker
    cdic, cdata = bC.to_bruker()
    assert_array_equal(bdata, cdata)
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write(td, cdic, cdata)
    rdic, rdata = ng.bruker.read(td)
    assert_array_equal(bdata, rdata)
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    shutil.rmtree(td) 

    # bruker -> pipe
    cdic, cdata = bC.to_pipe()
    assert_array_equal(pdata, cdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    os.remove(tf)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata)
    bpk = list(bad_pipe_keys)
    bpk.append("FDREALSIZE")    # NMRPipe corrects sizes for oversampling
    bpk.append("FDF2APOD")      
    bpk.append("FDF2TDSIZE")    
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    check_pdic(pdic, cdic, bpk, v=True)
    os.remove(tf)

    # pipe -> bruker
    cdic, cdata = pC.to_bruker()
    assert_array_equal(bdata, cdata)
    check_dic(bdic, cdic, bad_bruker_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write(td, cdic, cdata)
    rdic, rdata = ng.bruker.read(td)
    assert_array_equal(bdata, rdata)
    check_dic(bdic, cdic, bad_bruker_keys)
    shutil.rmtree(td)

def test_bruker_3d():
    """ 3D time bruker, pipe <-> bruker, pipe """
    # prepare Bruker converter
    bdic, bdata = ng.bruker.read(DATA_DIR + "bruker_3d")
    ubdic = ng.bruker.guess_udic(bdic, bdata)
    bC = ng.convert.converter()
    bC.from_bruker(bdic, bdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "bruker_3d/fid/test%03d.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # bruker -> bruker
    cdic, cdata = bC.to_bruker()
    assert_array_equal(bdata, cdata)
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write(td, cdic, cdata)
    rdic, rdata = ng.bruker.read(td)
    assert_array_equal(bdata, rdata)
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    shutil.rmtree(td) 

    # bruker -> pipe
    cdic, cdata = bC.to_pipe()
    assert_array_equal(pdata, cdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata)
    bpk = list(bad_pipe_keys)
    bpk.append("FDDISPMAX")     # nmrglue doesn't update the MIN/MAX values
    bpk.append("FDMIN")
    bpk.append("FDDISPMIN")
    bpk.append("FDSCALEFLAG")
    bpk.append("FDMAX")
    bpk.append("FDREALSIZE")    # NMRPipe corrects sizes for oversampling
    bpk.append("FDF2APOD")      
    bpk.append("FDF2TDSIZE")    
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    check_pdic(pdic, cdic, bpk, v=True)
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> bruker
    cdic, cdata = pC.to_bruker()
    assert_array_equal(bdata, cdata)
    check_dic(bdic, cdic, bad_bruker_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write(td, cdic, cdata)
    rdic, rdata = ng.bruker.read(td)
    assert_array_equal(bdata, rdata)
    check_dic(bdic, cdic, bad_bruker_keys)
    shutil.rmtree(td)

def test_sparky_2d():
    """ 2D freq sparky, pipe <-> sparky, pipe """
    # prepare Sparky converter
    sdic, sdata = ng.sparky.read(DATA_DIR + "sparky_2d/data.ucsf")
    ubdic = ng.sparky.guess_udic(sdic, sdata)
    sC = ng.convert.converter()
    sC.from_sparky(sdic, sdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "nmrpipe_2d/test.ft2")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # sparky -> sparky
    cdic, cdata = sC.to_sparky()
    assert_array_equal(sdata, cdata)
    check_sdic(sdic, cdic, bad_sparky_keys, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write(tf, cdic, cdata)
    rdic, rdata = ng.sparky.read(tf)
    assert_array_equal(sdata, rdata)
    check_sdic(sdic, rdic, bad_sparky_keys, v=True)
    os.remove(tf)

    # sparky -> pipe
    cdic, cdata = sC.to_pipe()
    assert_array_equal(pdata, cdata[:])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata[:])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    os.remove(tf)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata[:])
    bpk = list(bad_pipe_keys)
    bpk.append("FDF1TDSIZE")    # we lose all processing information
    bpk.append("FDF1APOD")
    bpk.append("FDF1ZF")
    bpk.append("FDF1C1")
    bpk.append("FDF1APODCODE")
    bpk.append("FDF1APODQ1")
    bpk.append("FDF1APODQ2")
    bpk.append("FDF1APODQ3")
    bpk.append("FDF2P0")
    bpk.append("FDF2APOD")
    bpk.append("FDF2ZF")
    bpk.append("FDF2APODCODE")
    bpk.append("FDF2APODQ1")
    bpk.append("FDF2APODQ2")
    bpk.append("FDF2APODQ3")
    bpk.append("FDSLICECOUNT") 
    bpk.append("FDREALSIZE")
    bpk.append("FDF2TDSIZE")
    
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata[:])
    check_pdic(pdic, cdic, bpk)
    os.remove(tf)

    # pipe -> sparky
    cdic, cdata = pC.to_sparky()
    assert_array_equal(sdata, cdata[:])
    check_sdic(sdic, cdic, bad_sparky_keys, True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write(tf, cdic, cdata)
    rdic, rdata = ng.sparky.read(tf)
    assert_array_equal(sdata, rdata[:])
    check_sdic(sdic, rdic, bad_sparky_keys, True)
    os.remove(tf)


def test_sparky_3d():
    """ 3D freq sparky, pipe <-> sparky, pipe """
    # prepare Sparky converter
    sdic, sdata = ng.sparky.read(DATA_DIR + "sparky_3d/data.ucsf")
    ubdic = ng.sparky.guess_udic(sdic, sdata)
    sC = ng.convert.converter()
    sC.from_sparky(sdic, sdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "nmrpipe_3d/ft/test%03d.ft3")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # sparky -> sparky
    cdic, cdata = sC.to_sparky()
    assert_array_equal(sdata, cdata)
    check_sdic(sdic, cdic, bad_sparky_keys, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write(tf, cdic, cdata)
    rdic, rdata = ng.sparky.read(tf)
    assert_array_equal(sdata, rdata)
    check_sdic(sdic, rdic, bad_sparky_keys, v=True)
    os.remove(tf)

    # sparky -> pipe
    cdic, cdata = sC.to_pipe()
    assert_array_equal(pdata, cdata[:])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata[:])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    for f in glob.glob(tf[:-4]+"*"):
        os.remove(f)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata[:])
    bpk = list(bad_pipe_keys)
    bpk.append("FDDISPMAX")     # nmrglue doesn't update the MIN/MAX values
    bpk.append("FDMIN")
    bpk.append("FDDISPMIN")
    bpk.append("FDSCALEFLAG")
    bpk.append("FDMAX")
    
    bpk.append("FDF1TDSIZE")    # we lose all processing information
    bpk.append("FDF1APOD")
    bpk.append("FDF1ZF")
    
    bpk.append("FDF2TDSIZE")
    bpk.append("FDF2APOD")
    bpk.append("FDF2ZF")
 
    bpk.append("FDF3TDSIZE")
    bpk.append("FDF3APOD")
    bpk.append("FDF3ZF")
    bpk.append("FDF3P0")
    bpk.append("FDF3P1")
    
    bpk.append("FDREALSIZE")
    
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata[:])
    check_pdic(pdic, cdic, bpk)
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> sparky
    cdic, cdata = pC.to_sparky()
    assert_array_equal(sdata, cdata[:])
    bsk = list(bad_sparky_keys)
    bsk.append("nucleus")           # the NMRPipe nucleus labels are passed
                                    # to sparky not 13C, 15N, etc
    check_sdic(sdic, cdic, bsk, True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write(tf, cdic, cdata)
    rdic, rdata = ng.sparky.read(tf)
    assert_array_equal(sdata, rdata[:])
    check_sdic(sdic, rdic, bsk, True)
    os.remove(tf)

def test_pipe_1d():
    """ 1D freq pipe <-> pipe """
    # prepare Pipe converter
    pdic, pdata = ng.pipe.read(DATA_DIR + "nmrpipe_1d/test.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata, cdata)
    check_pdic(pdic, cdic, bad_pipe_keys)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read(tf)
    assert_array_equal(pdata, rdata)
    check_pdic(pdic, cdic, bad_pipe_keys)
    os.remove(tf)

#### lowmemory tests ###
def test_sparky_2d_lowmem():
    """ 2D freq sparky, pipe <-> sparky, pipe low memory"""
    # prepare Sparky converter
    sdic, sdata = ng.sparky.read_lowmem(DATA_DIR + "sparky_2d/data.ucsf")
    ubdic = ng.sparky.guess_udic(sdic, sdata)
    sC = ng.convert.converter()
    sC.from_sparky(sdic, sdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read_lowmem(DATA_DIR + "nmrpipe_2d/test.ft2")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # sparky -> sparky
    cdic, cdata = sC.to_sparky()
    assert_array_equal(sdata[0:3, 0:4], cdata[0:3, 0:4])
    check_sdic(sdic, cdic, bad_sparky_keys, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.sparky.read_lowmem(tf)
    assert_array_equal(sdata[0:3, 0:4], rdata[0:3, 0:4])
    check_sdic(sdic, rdic, bad_sparky_keys, v=True)
    os.remove(tf)
    
    # sparky -> pipe
    cdic, cdata = sC.to_pipe()
    assert_array_equal(pdata[0:3, 0:4], cdata[0:3, 0:4])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".") 
    # XXX trace by trace writing from sparky is very slow
    ng.pipe.write(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:3, 0:4], rdata[0:3, 0:4])
    #check_pdic(pdic,cdic,v=True)   # XXX don't check dictionary
    os.remove(tf)
    
    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata[0:3, 0:4], cdata[0:3, 0:4])
    bpk = list(bad_pipe_keys)
    bpk.append("FDF1TDSIZE")    # we lose all processing information
    bpk.append("FDF1APOD")
    bpk.append("FDF1ZF")
    bpk.append("FDF1C1")
    bpk.append("FDF1APODCODE")
    bpk.append("FDF1APODQ1")
    bpk.append("FDF1APODQ2")
    bpk.append("FDF1APODQ3")
    bpk.append("FDF2P0")
    bpk.append("FDF2APOD")
    bpk.append("FDF2ZF")
    bpk.append("FDF2APODCODE")
    bpk.append("FDF2APODQ1")
    bpk.append("FDF2APODQ2")
    bpk.append("FDF2APODQ3")
    bpk.append("FDSLICECOUNT") 
    bpk.append("FDREALSIZE")
    bpk.append("FDF2TDSIZE")
    
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:3, 0:4], rdata[0:3, 0:4])
    check_pdic(pdic, cdic, bpk)
    os.remove(tf)

    # pipe -> sparky
    cdic, cdata = pC.to_sparky()
    assert_array_equal(sdata[0:3, 0:4], cdata[0:3, 0:4])
    check_sdic(sdic, cdic, bad_sparky_keys, True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.sparky.read_lowmem(tf)
    assert_array_equal(sdata[0:3, 0:4], rdata[0:3, 0:4])
    check_sdic(sdic, rdic, bad_sparky_keys, True)
    os.remove(tf)

def test_sparky_3d_lowmem():
    """ 3D freq sparky, pipe <-> sparky, pipe low memory"""
    # prepare Sparky converter
    sdic, sdata = ng.sparky.read_lowmem(DATA_DIR + "sparky_3d/data.ucsf")
    ubdic = ng.sparky.guess_udic(sdic, sdata)
    sC = ng.convert.converter()
    sC.from_sparky(sdic, sdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read_lowmem(DATA_DIR + "nmrpipe_3d/ft/test%03d.ft3")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # sparky -> sparky
    cdic, cdata = sC.to_sparky()
    assert_array_equal(sdata[0:3, 0:4, 0:5], cdata[0:3, 0:4, 0:5])
    check_sdic(sdic, cdic, bad_sparky_keys, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.sparky.read_lowmem(tf)
    assert_array_equal(sdata[0:3, 0:4, 0:5], rdata[0:3, 0:4, 0:5])
    check_sdic(sdic, rdic, bad_sparky_keys, v=True)
    os.remove(tf)

    # sparky -> pipe
    cdic, cdata = sC.to_pipe()
    assert_array_equal(pdata[0:3, 0:4, 0:5], cdata[0:3, 0:4, 0:5])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    # XXX trace by trace writing from sparky is very slow.
    #tf = tempfile.mktemp(dir=".") + "%03d"
    #ng.pipe.write_lowmem(tf, cdic, cdata)
    #rdic, rdata = ng.pipe.read_lowmem(tf)
    #assert_array_equal(pdata, rdata[:])
    #check_pdic(pdic, cdic)   # XXX don't check dictionary
    #for f in glob.glob(tf[:-4] + "*"):
    #    os.remove(f)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata[0:3, 0:4, 0:5], cdata[0:3, 0:4, 0:5])
    bpk = list(bad_pipe_keys)
    bpk.append("FDDISPMAX")     # nmrglue doesn't update the MIN/MAX values
    bpk.append("FDMIN")
    bpk.append("FDDISPMIN")
    bpk.append("FDSCALEFLAG")
    bpk.append("FDMAX")
    
    bpk.append("FDF1TDSIZE")    # we lose all processing information
    bpk.append("FDF1APOD")
    bpk.append("FDF1ZF")
    
    bpk.append("FDF2TDSIZE")
    bpk.append("FDF2APOD")
    bpk.append("FDF2ZF")
 
    bpk.append("FDF3TDSIZE")
    bpk.append("FDF3APOD")
    bpk.append("FDF3ZF")
    bpk.append("FDF3P0")
    bpk.append("FDF3P1")
    
    bpk.append("FDREALSIZE")
    
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:3, 0:4, 0:5], rdata[0:3, 0:4, 0:5])
    check_pdic(pdic, cdic, bpk)
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> sparky
    cdic, cdata = pC.to_sparky()
    assert_array_equal(sdata[0:3, 0:4, 0:5], cdata[0:3, 0:4, 0:5])
    bsk = list(bad_sparky_keys)
    bsk.append("nucleus")           # the NMRPipe nucleus labels are passed
                                    # to sparky not 13C, 15N, etc
    check_sdic(sdic, cdic, bsk, True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.sparky.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.sparky.read_lowmem(tf)
    assert_array_equal(sdata[0:3, 0:4, 0:5], rdata[0:3, 0:4, 0:5])
    check_sdic(sdic, rdic, bsk, True)
    os.remove(tf)

def test_bruker_2d_lowmem():
    """ 2D time bruker, pipe <-> bruker, pipe low memory"""
    # prepare Bruker converter
    bdic, bdata = ng.bruker.read_lowmem(DATA_DIR + "bruker_2d")
    ubdic = ng.bruker.guess_udic(bdic, bdata)
    bC = ng.convert.converter()
    bC.from_bruker(bdic, bdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read_lowmem(DATA_DIR + "bruker_2d/test.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # bruker -> bruker
    cdic, cdata = bC.to_bruker()
    assert_array_equal(bdata[0:3, 100:200], cdata[0:3, 100:200])
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write_lowmem(td, cdic, cdata)
    rdic, rdata = ng.bruker.read_lowmem(td)
    assert_array_equal(bdata[0:3, 100:200], rdata[0:3, 100:200])
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    shutil.rmtree(td) 

    # bruker -> pipe
    cdic, cdata = bC.to_pipe()
    assert_array_equal(pdata[0:3, 100:200], cdata[0:3, 100:200])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:3, 100:200], rdata[0:3, 100:200])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    os.remove(tf)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata[0:3, 100:200], cdata[0:3, 100:200])
    bpk = list(bad_pipe_keys)
    bpk.append("FDREALSIZE")    # NMRPipe corrects sizes for oversampling
    bpk.append("FDF2APOD")      
    bpk.append("FDF2TDSIZE")    
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:3, 100:200], rdata[0:3, 100:200])
    check_pdic(pdic, cdic, bpk, v=True)
    os.remove(tf)

    # pipe -> bruker
    cdic, cdata = pC.to_bruker()
    assert_array_equal(bdata[0:3, 100:200], cdata[0:3, 100:200])
    check_dic(bdic, cdic, bad_bruker_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write_lowmem(td, cdic, cdata)
    rdic, rdata = ng.bruker.read_lowmem(td)
    assert_array_equal(bdata[0:3, 100:200], rdata[0:3, 100:200])
    check_dic(bdic, cdic, bad_bruker_keys)
    shutil.rmtree(td)

def test_bruker_3d_lowmem():
    """ 3D time bruker, pipe <-> bruker, pipe low memory"""
    # prepare Bruker converter
    bdic, bdata = ng.bruker.read_lowmem(DATA_DIR + "bruker_3d")
    ubdic = ng.bruker.guess_udic(bdic, bdata)
    bC = ng.convert.converter()
    bC.from_bruker(bdic, bdata, ubdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read_lowmem(DATA_DIR + "bruker_3d/fid/test%03d.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # bruker -> bruker
    cdic, cdata = bC.to_bruker()
    assert_array_equal(bdata[0:2, 0:3, 100:200], cdata[0:2, 0:3, 100:200])
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write_lowmem(td, cdic, cdata)
    rdic, rdata = ng.bruker.read_lowmem(td)
    assert_array_equal(bdata[0:2, 0:3, 100:200], rdata[0:2, 0:3, 100:200])
    check_dic(bdic, cdic, bad_bruker_keys, v=True)
    shutil.rmtree(td) 

    # bruker -> pipe
    cdic, cdata = bC.to_pipe()
    assert_array_equal(pdata[0:2, 0:3, 100:200], cdata[0:2, 0:3, 100:200])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:2, 0:3, 100:200], rdata[0:2, 0:3, 100:200])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata[0:2, 0:3, 100:200], cdata[0:2, 0:3, 100:200])
    bpk = list(bad_pipe_keys)
    bpk.append("FDDISPMAX")     # nmrglue doesn't update the MIN/MAX values
    bpk.append("FDMIN")
    bpk.append("FDDISPMIN")
    bpk.append("FDSCALEFLAG")
    bpk.append("FDMAX")
    bpk.append("FDREALSIZE")    # NMRPipe corrects sizes for oversampling
    bpk.append("FDF2APOD")      
    bpk.append("FDF2TDSIZE")    
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:2, 0:3, 100:200], rdata[0:2, 0:3, 100:200])
    check_pdic(pdic, cdic, bpk, v=True)
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> bruker
    cdic, cdata = pC.to_bruker()
    assert_array_equal(bdata[0:2, 0:3, 100:200], cdata[0:2, 0:3, 100:200])
    check_dic(bdic, cdic, bad_bruker_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.bruker.write_lowmem(td, cdic, cdata)
    rdic, rdata = ng.bruker.read_lowmem(td)
    assert_array_equal(bdata[0:2, 0:3, 100:200], rdata[0:2, 0:3, 100:200])
    check_dic(bdic, cdic, bad_bruker_keys)
    shutil.rmtree(td)

def test_agilent_2d_lowmem():
    """ 2D time agilent, pipe <-> agilent, pipe low memory"""
    # prepare Varian converter
    vdic, vdata = ng.varian.read_lowmem(DATA_DIR + "agilent_2d")
    uvdic = ng.varian.guess_udic(vdic, vdata)
    vC = ng.convert.converter()
    vC.from_varian(vdic, vdata, uvdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read_lowmem(DATA_DIR + "agilent_2d/test.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # varian -> varian
    cdic, cdata = vC.to_varian()
    assert_array_equal(vdata[0:4, 0:20], cdata[0:4, 0:20])
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write_lowmem(td, cdic, cdata)
    rdic, rdata = ng.varian.read_lowmem(td)
    assert_array_equal(vdata[0:4, 0:20], rdata[0:4, 0:20])
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

    # varian -> pipe
    cdic, cdata = vC.to_pipe()
    assert_array_equal(pdata[0:4, 0:20], cdata[0:4, 0:20])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:4, 0:20], rdata[0:4, 0:20])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    os.remove(tf)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata[0:4, 0:20], cdata[0:4, 0:20])
    check_pdic(pdic, cdic, bad_pipe_keys)
    # write and readback
    tf = tempfile.mktemp(dir=".")
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:4, 0:20], rdata[0:4, 0:20])
    check_pdic(pdic, cdic, bad_pipe_keys)
    os.remove(tf)

    # pipe -> varian
    cdic, cdata = pC.to_varian()
    assert_array_equal(vdata[0:4, 0:20], cdata[0:4, 0:20])
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write_lowmem(td, cdic, cdata)
    rdic, rdata = ng.varian.read_lowmem(td)
    assert_array_equal(vdata[0:4, 0:20], rdata[0:4, 0:20])
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

def test_agilent_3d_lowmem():
    """ 3D time agilent, pipe <-> agilent, pipe low memory"""
    # prepare Agilent converter
    vdic, vdata = ng.varian.read_lowmem(DATA_DIR + "agilent_3d")
    uvdic = ng.varian.guess_udic(vdic, vdata)
    vC = ng.convert.converter()
    vC.from_varian(vdic, vdata, uvdic)

    # prepare Pipe converter
    pdic, pdata = ng.pipe.read_lowmem(DATA_DIR + "agilent_3d/data/test%03d.fid")
    updic = ng.pipe.guess_udic(pdic, pdata)
    pC = ng.convert.converter()
    pC.from_pipe(pdic, pdata, updic)

    # agilent -> agilent
    cdic, cdata = vC.to_varian()
    assert_array_equal(vdata[0:3, 0:4, 0:20], cdata[0:3, 0:4, 0:20])
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write_lowmem(td, cdic, cdata)
    rdic, rdata = ng.varian.read_lowmem(td)
    assert_array_equal(vdata[0:3, 0:4, 0:20], rdata[0:3, 0:4, 0:20])
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

    # agilent -> pipe
    cdic, cdata = vC.to_pipe()
    assert_array_equal(pdata[0:3, 0:4, 0:20], cdata[0:3, 0:4, 0:20])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:3, 0:4, 0:20], rdata[0:3, 0:4, 0:20])
    #check_pdic(pdic,cdic)   # XXX don't check dictionary
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> pipe
    cdic, cdata = pC.to_pipe()
    assert_array_equal(pdata[0:3, 0:4, 0:20], cdata[0:3, 0:4, 0:20])
    bpk = list(bad_pipe_keys)
    bpk.append("FDDISPMAX")     # nmrglue doesn't update the MIN/MAX values
    bpk.append("FDMIN")
    bpk.append("FDDISPMIN")
    bpk.append("FDSCALEFLAG")
    bpk.append("FDMAX")
    check_pdic(pdic, cdic, bpk, v=True)
    # write and readback
    tf = tempfile.mktemp(dir=".") + "%03d"
    ng.pipe.write_lowmem(tf, cdic, cdata)
    rdic, rdata = ng.pipe.read_lowmem(tf)
    assert_array_equal(pdata[0:3, 0:4, 0:20], rdata[0:3, 0:4, 0:20])
    check_pdic(pdic, cdic, bpk)
    for f in glob.glob(tf[:-4] + "*"):
        os.remove(f)

    # pipe -> varian
    cdic, cdata = pC.to_varian()
    assert_array_equal(vdata[0:3, 0:4, 0:20], cdata[0:3, 0:4, 0:20])
    check_dic(vdic, cdic, bad_varian_keys)
    # write and readback
    td = tempfile.mkdtemp(dir=".")
    ng.varian.write_lowmem(td, cdic, cdata)
    rdic, rdata = ng.varian.read_lowmem(td)
    assert_array_equal(vdata[0:3, 0:4, 0:20], rdata[0:3, 0:4, 0:20])
    check_dic(vdic, cdic, bad_varian_keys)
    shutil.rmtree(td)

# To skip a particular test, uncomment the test stub below

from nose.exc import SkipTest

#def test_agilent_1d():
#    """ 1D time agilent, pipe <-> agilent, pipe """
#    raise SkipTest

#def test_agilent_2d():
#    """ 2D time agilent, pipe <-> agilent, pipe """
#    raise SkipTest

#def test_agilent_3d():
#    """ 3D time agilent, pipe <-> agilent, pipe """
#    raise SkipTest

#def test_bruker_1d():
#    """ 1D time bruker, pipe <-> bruker, pipe """
#    raise SkipTest

#def test_bruker_2d():
#    """ 2D time bruker, pipe <-> bruker, pipe """
#    raise SkipTest

#def test_bruker_3d():
#    """ 3D time bruker, pipe <-> bruker, pipe """
#    raise SkipTest

#def test_pipe_1d():
#    """ 1D freq pipe <-> pipe """
#    raise SkipTest

#def test_sparky_2d():
#    """ 2D freq sparky, pipe <-> sparky, pipe """
#    raise SkipTest

#def test_sparky_3d():
#    """ 3D freq sparky, pipe <-> sparky, pipe """
#    raise SkipTest

#def test_sparky_2d_lowmem():
#    """ 2D freq sparky, pipe <-> sparky, pipe low memory"""
#    raise SkipTest

#def test_sparky_3d_lowmem():
#    """ 3D freq sparky, pipe <-> sparky, pipe low memory"""
#    raise SkipTest

#def test_bruker_2d_lowmem():
#    """ 2D time bruker, pipe <-> bruker, pipe low memory"""
#    raise SkipTest

#def test_bruker_3d_lowmem():
#    """ 3D time bruker, pipe <-> bruker, pipe low memory"""
#    raise SkipTest

#def test_agilent_2d_lowmem():
#    """ 2D time agilent, pipe <-> agilent, pipe low memory"""
#    raise SkipTest

#def test_agilent_3d_lowmem():
#    """ 3D time agilent, pipe <-> agilent, pipe low memory"""
#    raise SkipTest

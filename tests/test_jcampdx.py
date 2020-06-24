""" Tests for the fileio.jcampdx submodule """

import os
import numpy as np
import nmrglue as ng
from setup import DATA_DIR


def test_jcampdx1():
    '''JCAMP-DX read: format testset'''
    cases = []

    # these 4 files have exactly the same datapoints
    cases.append("BRUKAFFN.DX")  # 0
    cases.append("BRUKPAC.DX")  # 1
    cases.append("BRUKSQZ.DX")  # 2
    cases.append("TEST32.DX")  # 3

    # the rest have the same data than the first ones but slightly
    # different values probably due to some processing steps
    cases.append("BRUKNTUP.DX")  # 4
    cases.append("BRUKDIF.DX")  # 5
    cases.append("TESTSPEC.DX")  # 6
    cases.append("TESTNTUP.DX")  # 7

    npoints_target = 16384

    # target data values are from the BRUKAFFN.DX file
    # which is most human-readable

    # some values from the beginning:
    target_first20 = [
        2259260, -5242968, -7176216, -1616072,
        10650432, 4373926, -3660824, 2136488,
        8055988, 1757248, 3559312, 1108422,
        -5575546, -233168, -1099612, -4657542,
        -4545530, -1996712, -5429568, -7119772]

    # some values from the middle
    target_7144_7164 = [
        4613558, 9603556, 11823620, 3851634,
        14787192, 17047672, 34585306, 69387092,
        307952794, 542345870, 143662704, 52472738,
        29157730, 12017988, 10142310, -6518692,
        11292386, 4692342, 2839598, 4948336]

    # some values from the end:
    target_last20 = [
        -2731004, 2823836, 1542934, 8096410,
        1143092, -5356388, 4028632, 121858,
        3829486, 5562002, -3851528, 919686,
        1060812, -4446420, -716388, 2080534,
        7145886, 11400102, 5012684, 1505988]

    for i, case in enumerate(cases):

        print(case)

        # read
        casepath = os.path.join(DATA_DIR, "jcampdx", case)
        dic, data = ng.jcampdx.read(casepath)

        # check some dic entries:
        assert "DATATYPE" in dic
        assert "DATA TYPE" not in dic
        assert "END" not in dic
        if "BRUK" in case:
            assert dic["DATATYPE"][0] == "NMR Spectrum"
            assert dic["$SOLVENT"][0] == "<MeOH>"  # single $ is not comment
            assert "Bruker" not in dic
            assert "Bruker NMR JCAMP-DX V1.0" in dic["_comments"]  # comment

        # check data point counts:
        if "NTUP" not in case:  # no ##NPOINTS in NTUPLES format
            npoints_dic = int(dic["NPOINTS"][0])
            assert npoints_dic == npoints_target
            npoints_read = len(data)
        else:  # NTUPLES has both real & imag arrays
            npoints_read = len(data[0])  # R
            assert len(data[1]) == npoints_target  # I
            data = data[0]
        assert npoints_read == npoints_target

        # check data:
        epsilon_e = 1e-9
        epsilon_r = 15000
        for idx in range(20):
            print(target_first20[idx], data[idx])
            if i < 4:  # exactly same data
                assert np.abs(target_first20[idx]-data[idx]) < epsilon_e
            else:  # roughly same data
                assert np.abs(target_first20[idx]-data[idx]) < epsilon_r

        for idx in range(20):
            dslice = data[7144:7164]
            print(target_7144_7164[idx], dslice[idx])
            if i < 4:  # exactly same data
                assert np.abs(target_7144_7164[idx]-dslice[idx]) < epsilon_e
            else:  # roughly same data
                assert np.abs(target_7144_7164[idx]-dslice[idx]) < epsilon_r

        for idx in range(20):
            dslice = data[-20:]
            print(target_last20[idx], dslice[idx])
            if i < 4:  # exactly same data
                assert np.abs(target_last20[idx]-dslice[idx]) < epsilon_e
            else:  # roughly same data
                assert np.abs(target_last20[idx]-dslice[idx]) < epsilon_r

        # check udic:
        udic = ng.jcampdx.guess_udic(dic, data)
        assert np.abs(udic[0]["obs"]-100.4) < epsilon_e
        assert np.abs(udic[0]["sw"]-24038.5) < epsilon_e
        assert udic[0]["size"] == npoints_target
        assert udic[0]["label"] == "13C"


def test_jcampdx2():
    '''JCAMP-DX read: miscellaneous files '''
    cases = []
    # check the following values from each read:
    # npoints, first, last, freq, sweep
    # note: first and last are raw values from datalines for convenience,
    # i.e. not scaled with YFACTORS
    cases.append(("TESTFID.DX", 16384, 573, -11584, 100.4, 0.6815317))
    cases.append(("bruker1.dx", 16384, -5, -51, 200.13, 4098.3606557377))
    cases.append(("bruker2.dx", 16384, 42, 422, 300.13336767, 6024.096385479))
    cases.append(("bruker3.dx", 16384, 22, -313, 300.13336729, 6024.096385479))
    cases.append(("aug07.dx",
                  4096, -22288, -25148, 400.13200065, 4006.41025641027))
    cases.append(("aug07b.dx",
                  4096, -324909, -205968, 400.13200065, 4006.41025641027))
    cases.append(("aug07c.dx",
                  4096, -322709, -64216, 400.13200065, 4006.41025641027))
    cases.append(("aug07d.dx",
                  4096, -21208, 3029, 400.13200065, 4006.41025641027))
    cases.append(("aug07e.dx",
                  4096, -501497, 79397, 400.13200065, 4006.41025641027))

    epsilon = 1e-9

    for case in cases:
        print(case[0])
        # read
        casepath = os.path.join(DATA_DIR, "jcampdx", case[0])
        dic, data = ng.jcampdx.read(casepath)
        if isinstance(data, list):
            data = data[0]  # for data with both R&I, check only R

        # since first and last are raw values, do yfactor
        # back-scaling here
        is_ntuples = ng.jcampdx.get_is_ntuples(dic)
        if is_ntuples:
            yfactor_r, _yfactor_i = ng.jcampdx.find_yfactors(dic)
            data = data / yfactor_r
        else:
            yfactor = float(dic["YFACTOR"][0])
            data = data / yfactor

        # check data
        assert len(data) == case[1]
        assert np.abs(data[0]-case[2]) < epsilon
        assert np.abs(data[-1]-case[3]) < epsilon

        # check udic
        udic = ng.jcampdx.guess_udic(dic, data)
        assert np.abs(udic[0]["obs"]-case[4]) < epsilon
        assert np.abs(udic[0]["sw"]-case[5]) < epsilon


def test_jcampdx_dicstructure():
    '''JCAMP-DX read: ensure correct dic structure '''

    casepath = os.path.join(DATA_DIR, "jcampdx", "dicstructure.jdx")
    dic, data = ng.jcampdx.read(casepath)

    # check data
    assert len(data) == 8
    assert data[-1] == 8.0

    # check dic:
    assert dic["DATATYPE"][0] == "NMR SPECTRUM"
    assert dic[".OBSERVENUCLEUS"][0] == "^1H"
    assert "_datatype_LINK" in dic
    assert "_datatype_NMRPEAKTABLE" in dic
    assert "_datatype_NMRSPECTRUM" not in dic
    assert "_datatype_NA" in dic
    assert len(dic["_datatype_LINK"]) == 2
    assert len(dic["_datatype_NMRPEAKTABLE"]) == 1
    assert len(dic["_datatype_NA"]) == 1
    assert dic["_datatype_LINK"][0]["DUMMYENTRY"][0] == "1.0"
    assert dic["_datatype_LINK"][1]["DUMMYENTRY"][0] == "2.0"
    assert dic["_datatype_NMRPEAKTABLE"][0]["DUMMYENTRY"][0] == "3.0"
    assert dic["_datatype_NA"][0]["DUMMYENTRY"][0] == "4.0"
    assert dic["_datatype_NA"][0]["_comments"][0] == "comment line"
    assert dic["_datatype_NA"][0]["_comments"][1] == "another comment"

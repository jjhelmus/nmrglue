""" Tests for the fileio.rs2d submodule """

import os
import numpy as np
import nmrglue as ng
from setup import DATA_DIR


def test_rs2d():
    '''RS2D read: format testset. Expects to find the provided test cases under 
    data/rs2d/Installer_data'''

    # case tuple is:
    #  0       1      2       3,          4            5        6
    # (folder, shape, is_fid, is_complex, basefreqMHz, sweepHz, offsetHz,
    #  7        8        9
    #  nucleus, solvent, temperature)

    # tested values are collected manually from header.xml files

    cases = []

    # -------------------------------------------------------------------------

    # 1D 1H (FID)
    bf_hz = 4.001293456641604E8  # BASE_FREQ
    bf_mhz = bf_hz/1e6
    # OBSERVED_FREQUENCY (for sw conversion!)
    tf_hz = 4.0013094618154305E8
    tf_mhz = tf_hz/1e6
    car = tf_hz - bf_hz
    cases.append(("1", (16384,), True, True, (bf_mhz,),
                 (10.000439579964194 * tf_mhz,), (car,), ("1H",), "CDCl3", 298))

    # 1D 1H (PROCESSED)
    cases.append(("1/Proc/0", (32768,), False, True, (bf_mhz,),
                 (10.000439579964194 * tf_mhz,), (car,), ("1H",), "CDCl3", 298))

    # -------------------------------------------------------------------------

    # 1D 13C (FID)
    bf_hz = 1.006126140544557E8
    bf_mhz = bf_hz/1e6
    tf_hz = 1.006216691897206E8  # OBSERVED_FREQUENCY (for sw conversion!)
    tf_mhz = tf_hz/1e6
    car = tf_hz - bf_hz
    cases.append(("2", (32768,), True, True, (bf_mhz,),
                 (200.00598044404157 * tf_mhz,), (car,), ("13C",), "Benzene", 298))

    # 1D 13C (PROCESSED)
    cases.append(("2/Proc/0", (65536,), False, True, (bf_mhz,),
                  (200.00598044404157 * tf_mhz,), (car,), ("13C",), "Benzene", 298))

    # -------------------------------------------------------------------------

    # 2D HSQC (FID)
    bf_hz_1 = 4.001293856773854E8
    bf_mhz_1 = bf_hz_1/1e6
    bf_hz_2 = 1.006126140544557E8
    bf_mhz_2 = bf_hz_2/1e6

    tf_hz_1 = 4.0013138632431376E8  # OBSERVED_FREQUENCY (for sw conversion!)
    tf_mhz_1 = tf_hz_1/1e6
    tf_hz_2 = 1.0062066306358005E8  # TRANSMIT_FREQ_2 (for sw conversion!)
    tf_mhz_2 = tf_hz_2/1e6
    car_1 = tf_hz_1 - bf_hz_1
    car_2 = tf_hz_2 - bf_hz_2
    cases.append(("17", (256, 1024), True, True, (bf_mhz_2, bf_mhz_1),
                  (199.9544840071294 * tf_mhz_2, 11.993142972152109 * tf_mhz_1),
                 (car_2, car_1), ("13C", "1H"), "Benzene", 298))

    # 2D HSQC (PROCESSED)
    cases.append(("17/Proc/0", (256, 2048), False, True, (bf_mhz_2, bf_mhz_1),
                  (199.9544840071294 * tf_mhz_2, 11.993142972152109 * tf_mhz_1),
                 (car_2, car_1), ("13C", "1H"), "Benzene", 298))

    # -------------------------------------------------------------------------

    epsilon = 1e-9

    for case in cases:

        folder = case[0]
        print("RS2D case:", folder)

        # read
        casepath = os.path.join(DATA_DIR, "rs2d", "Installer_Data", folder)
        dic, data = ng.rs2d.read(casepath)

        # check data:
        assert data.shape == case[1]

        # check udic:
        udic = ng.rs2d.guess_udic(dic, data)
        for dim in range(len(data.shape)):
            assert udic[dim]["time"] is case[2]
            assert udic[dim]["freq"] is not case[2]
            assert np.abs(udic[dim]["obs"] - case[4][dim]) < epsilon
            assert np.abs(udic[dim]["sw"] - case[5][dim]) < epsilon
            assert np.abs(udic[dim]["car"] - case[6][dim]) < epsilon
            assert udic[dim]["size"] == case[1][dim]
            assert udic[dim]["label"] == case[7][dim]

        # data type (complex/real) tested only for the first dimension
        assert udic[len(data.shape)-1]["complex"] is case[3]

        # check some rs2d dict values:
        assert dic["SOLVENT"]["value"] == case[8]
        temperature = float(dic["SAMPLE_TEMPERATURE"]["value"])
        assert np.abs(temperature - case[9]) < epsilon

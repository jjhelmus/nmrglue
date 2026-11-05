"""Pipe digital-filter metadata through format conversion."""

import numpy as np

import nmrglue as ng


def _one_d_pipe():
    udic = ng.fileiobase.create_blank_udic(1)
    dic = ng.pipe.create_dic(udic)
    data = np.zeros(8, dtype="complex64")
    return dic, data


def test_pipe_round_trip_keeps_fddmxval():
    dic, data = _one_d_pipe()
    dic["FDDMXVAL"] = 71.625

    converter = ng.convert.converter()
    converter.from_pipe(dic, data)
    out_dic, out_data = converter.to_pipe()

    assert out_dic["FDDMXVAL"] == 71.625
    assert np.array_equal(data, out_data)


def test_pipe_round_trip_keeps_negative_fddmxval():
    dic, data = _one_d_pipe()
    dic["FDDMXVAL"] = -12.5

    converter = ng.convert.converter()
    converter.from_pipe(dic, data)
    out_dic, _out_data = converter.to_pipe()

    assert out_dic["FDDMXVAL"] == -12.5


def test_missing_fddmxval_stays_at_the_pipe_default():
    dic, data = _one_d_pipe()
    del dic["FDDMXVAL"]

    converter = ng.convert.converter()
    converter.from_pipe(dic, data)
    out_dic, _out_data = converter.to_pipe()

    assert out_dic["FDDMXVAL"] == 0.0


def test_later_non_pipe_load_does_not_keep_fddmxval():
    dic, data = _one_d_pipe()
    dic["FDDMXVAL"] = 71.625
    udic = ng.pipe.guess_udic(dic, data)

    converter = ng.convert.converter()
    converter.from_pipe(dic, data)
    converter.from_universal(udic, data)
    out_dic, _out_data = converter.to_pipe()

    assert out_dic["FDDMXVAL"] == 0.0

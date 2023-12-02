"""
Functions for reading and writing Jeol JDF files

"""

import locale
import io

__developer_info__ = """
Jeol data format information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:

"""

import os
from warnings import warn

import numpy as np

from . import fileiobase
from ..process import proc_base

import struct


def read(fname):
    """
    _summary_

    Parameters
    ----------
    fname : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """

    with open(fname, "rb") as f:
        data = f.read()

    header = parse_jeol(data)

    return header


class IOBuffer:

    def __init__(self, buffer, conversion_table, endian=None):
        """
        _summary_

        Parameters
        ----------
        buffer : _type_
            _description_
        conversion_table : _type_
            _description_
        endian : _type_, optional
            _description_, by default None

        """
        self.buffer = buffer
        self.conversion_table = conversion_table
        self.position = 0
        self.set_big_endian() # default for the header

    def set_big_endian(self):
        self.endian = 'big'
        self.e = ">" 

    def set_little_endian(self):
        self.endian = "little"
        self.e = "<"

    def skip(self, count):
        self.position += count

    def read_chars(self, count):
        chars = struct.unpack_from(f"{self.e}{count}s", self.buffer, self.position)
        self.position += count
        return chars[0].decode("utf-8")

    def read_int8(self):
        value = struct.unpack_from(f"{self.e}b", self.buffer, self.position)
        self.position += 1
        return value[0]

    def read_uint8(self):
        value = struct.unpack_from(f"{self.e}B", self.buffer, self.position)
        self.position += 1
        return value[0]

    def read_uint16(self):
        value = struct.unpack_from(f"{self.e}H", self.buffer, self.position)
        self.position += 2
        return value[0]

    def read_int16(self):
        value = struct.unpack_from(f"{self.e}h", self.buffer, self.position)
        self.position += 2
        return value[0]

    def read_uint32(self):
        value = struct.unpack_from(f"{self.e}I", self.buffer, self.position)
        self.position += 4
        return value[0]

    def read_int32(self):
        value = struct.unpack_from(f"{self.e}i", self.buffer, self.position)
        self.position += 4
        return value[0]

    def read_float64(self):
        value = struct.unpack_from(f"{self.e}d", self.buffer, self.position)
        self.position += 8
        return value[0]

    def read_byte(self):
        """ reads a single byte """
        value = struct.unpack_from(f"B", self.buffer, self.position)
        self.position += 1
        return value[0]

    def read_boolean(self):
        """ reads a boolean """
        value = struct.unpack_from(f"B", self.buffer, self.position)
        self.position += 1
        return bool(value[0])

    def get_array(self, read_func, count, *args, **kwargs):
        return [getattr(self, read_func)(*args, **kwargs) for _ in range(count)]

    def get_unit(self, size):
        """
        parses the byte that gives the unit of the value of a 
        parameter and gives the power

        Parameters
        ----------
        size : _type_
            _description_

        Returns
        -------
        _type_
            _description_
        """
        unit = []
        for i in range(size):
            byte = self.read_byte()
            prefix = self.conversion_table.prefix[byte >> 4]
            power = byte & 0b00001111
            base = self.conversion_table.base[self.read_int8()]
            unit.append({"prefix": prefix, "power": power, "base": base})
        return unit

    def get_string(self, count):
        return self.read_chars(count).replace("\x00", "")


def parse_jeol(buffer):
    buffer = IOBuffer(buffer, conversion_table=ConversionTable)
    buffer, header = read_header(buffer, conversion_table=ConversionTable)
    buffer, params = read_parameters(
        buffer, header["param_start"], header["endian"], ConversionTable
    )

    return {"header": header, "parameters": params}


def read_header(buffer, conversion_table):
    """
    _summary_

    Parameters
    ----------
    buffer : _type_
        _description_
    conversion_table : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    t = conversion_table  # short for convinience
    header = {}

    header["file_identifier"] = buffer.read_chars(8)
    header["endian"] = t.endianness[buffer.read_int8()]
    header["major_version"] = buffer.read_uint8()
    header["minor_version"] = buffer.read_uint16()
    header["ndim"] = buffer.read_uint8()

    info = buffer.read_byte()
    header["dim_exist"] = [bool(int(x)) for x in format(info, "08b")]

    info = buffer.read_byte()
    header["type"] = t.data_type[info >> 6]
    header["format"] = t.data_type[info & 0b00111111]

    header["instrument"] = t.instruments[buffer.read_int8()]
    header["translate"] = [buffer.read_int8() for i in range(8)]
    header["axis_type"] = [t.axis_type[i] for i in buffer.get_array("read_int8", 8)]
    header["units"] = buffer.get_unit(8)
    header["title"] = buffer.get_string(124)

    info = []
    for i in buffer.get_array("read_uint8", 4):
        # print(bytes(i))
        info.append(t.data_axis_ranged[i >> 4])
        info.append(t.data_axis_ranged[i & 0b00001111])
    header["data_axis_ranged"] = info

    header["npoints"] = buffer.get_array("read_uint32", 8)
    header["offset_start"] = buffer.get_array("read_uint32", 8)
    header["offset_stop"] = buffer.get_array("read_uint32", 8)
    header["axis_start"] = buffer.get_array("read_float64", 8)
    header["axis_stop"] = buffer.get_array("read_float64", 8)

    info = buffer.get_array("read_byte", 4)
    header["creation_date"] = {
        "year": 1990 + (info[0] >> 1),
        "month": ((info[0] << 3) & 0b00001000) + (info[1] >> 5),
        "day": info[2] & 0b00011111,
    }

    info = buffer.get_array("read_byte", 4)
    header["revision_date"] = {
        "year": 1990 + (info[0] >> 1),
        "month": ((info[0] << 3) & 0b00001000) + (info[1] >> 5),
        "day": info[2] & 0b00011111,
    }

    header["node_name"] = buffer.get_string(16)
    header["site"] = buffer.get_string(128)
    header["author"] = buffer.get_string(128)
    header["comment"] = buffer.get_string(128)
    header["axis_title"] = buffer.get_array("get_string", 8, 32)
    header["base_frequency"] = buffer.get_array("read_float64", 8)
    header["zero_point"] = buffer.get_array("read_float64", 8)
    header["reversed"] = buffer.get_array("read_boolean", 8)

    buffer.skip(3)
    header["annotation_flag"] = bool(buffer.read_byte() >> 7)
    header["history_used"] = buffer.read_uint32()
    header["history_length"] = buffer.read_uint32()
    header["param_start"] = buffer.read_uint32()
    header["param_length"] = buffer.read_uint32()
    header["list_start"] = buffer.get_array("read_uint32", 8)
    header["list_length"] = buffer.get_array("read_uint32", 8)
    header["data_start"] = buffer.read_uint32()
    header["data_length"] = (buffer.read_uint32() << 32) | (buffer.read_uint32())
    header["context_start"] = (buffer.read_uint32() << 32) | (buffer.read_uint32())
    header["context_length"] = buffer.read_uint32()
    header["annotation_start"] = (buffer.read_uint32() << 32) | (buffer.read_uint32())
    header["annotation_length"] = buffer.read_uint32()
    header["size"] = (buffer.read_uint32() << 32) | (buffer.read_uint32())
    header["unit_location"] = buffer.get_array(("read_uint8"), 8)

    info = []
    for _ in range(2):
        unit = []
        scaler = buffer.read_int16()
        for _ in range(5):
            byte = buffer.read_int16()
            unit.append(byte)
        info.append({"scaler": scaler, "unit": unit})
    header["compound_unit"] = info

    return buffer, header


def read_parameters(buffer, param_start, endianness, conversion_table):
    """
    _summary_

    Parameters
    ----------
    buffer : _type_
        _description_
    param_start : _type_
        _description_
    endianness : _type_
        _description_
    conversion_table : _type_
        _description_

    Returns
    -------
    _type_
        _description_
    """
    t = conversion_table

    print(buffer.position)
    if endianness == "little_endian":
        buffer.set_little_endian()

    buffer.position = param_start
    print(buffer.position, param_start)

    params = {
        "parameter_size": buffer.read_uint32(),
        "low_index": buffer.read_uint32(),
        "high_index": buffer.read_uint32(),
        "total_size": buffer.read_uint32(),
    }
    print(params)

    param_array = []
    for p in range(params["high_index"]):
        buffer.skip(4)
        scaler = buffer.read_int16()
        unit = buffer.get_unit(5)
        # print(unit)
        buffer.skip(16)
        value_type = t.value_type[buffer.read_int32()]
        # print(p)
        buffer.position -= 20
        # print(buffer.position)

        value = None
        if value_type == "string":
            value = buffer.get_string(16).replace(" ", "")

        elif value_type == "integer":
            value = buffer.read_int32()
            buffer.skip(12)

        elif value_type == "float":
            value = buffer.read_float64()
            buffer.skip(8)

        elif value_type == "complex":
            value = buffer.read_float64() + 1j * buffer.read_float64()

        elif value_type == "infinity":
            value = buffer.read_int32()
            buffer.skip(12)

        else:
            buffer.skip(16)

        buffer.skip(4)

        name = buffer.get_string(28).replace(" ", "")

        param_array.append(
            {
                "name": name.lower(),
                "scaler": scaler,
                "unit": unit,
                "value": value,
                "value_type": value_type,
            }
        )

    params["array"] = param_array

    return buffer, params


# Conversion Table for Jeol Data Format
class ConversionTable:
    endianness = {
        0: "big_endian",
        1: "little_endian",
    }

    instruments = {
        0: None,
        1: "gsx",
        2: "alpha",
        3: "eclipse",
        4: "mass_spec",
        5: "compiler",
        6: "other_nmr",
        7: "unknown",
        8: "gemini",
        9: "unity",
        10: "aspect",
        11: "ux",
        12: "felix",
        13: "lambda",
        14: "ge_1280",
        15: "ge_omega",
        16: "chemagnetics",
        17: "cdff",
        18: "galactic",
        19: "triad",
        20: "generic_nmr",
        21: "gamma",
        22: "jcamp_dx",
        23: "amx",
        24: "dmx",
        25: "eca",
        26: "alice",
        27: "nmrpipe",
        28: "simpson",
    }

    data_type = {
        0: "float64",
        1: "float32",
        2: "reserved",
        3: "reserved",
    }

    data_format = {
        1: "one_d",
        2: "two_d",
        3: "three_d",
        4: "four_d",
        5: "five_d",
        6: "six_d",
        7: "seven_d",
        8: "eight_d",
        9: "not for NMR data formats",
        10: "not for NMR data formats",
        11: "not for NMR data formats",
        12: "small_two_d",
        13: "small_three_d",
        14: "small_four_d",
    }

    axis_type = {
        0: None,
        1: "real",
        2: "tppi",
        3: "complex",
        4: "real_complex",
        5: "envelope",
    }

    prefix = {
        -8: "yotta",
        -6: "exa",
        -7: "zetta",
        -5: "pecta",
        -4: "tera",
        -3: "giga",
        -2: "mega",
        -1: "kilo",
        0: "none",
        1: "milli",
        2: "micro",
        3: "nano",
        4: "pico",
        5: "femto",
        6: "atto",
        7: "zepto",
        15: "None"
    }

    unit_prefix_table = {
        "Yotta": 24,
        "Exa": 21,
        "Zetta": 18,
        "Pecta": 15,
        "Tera": 12,
        "Giga": 9,
        "Mega": 6,
        "Kilo": 3,
        "None": 0,
        "Milli": -3,
        "Micro": -6,
        "Nano": -9,
        "Pico": -12,
        "Femto": -15,
        "Atto": -18,
        "Zepto": -21,
    }

    base = {
        0: None,
        1: "abundance",
        2: "ampere",
        3: "candela",
        4: "celsius",
        5: "coulomb",
        6: "degree",
        7: "electronvolt",
        8: "farad",
        9: "sievert",
        10: "gram",
        11: "gray",
        12: "henry",
        13: "hertz",
        14: "kelvin",
        15: "joule",
        16: "liter",
        17: "lumen",
        18: "lux",
        19: "meter",
        20: "mole",
        21: "newton",
        22: "ohm",
        23: "pascal",
        24: "percent",
        25: "point",
        26: "ppm",
        27: "radian",
        28: "second",
        29: "siemens",
        30: "steradian",
        31: "tesla",
        32: "volt",
        33: "watt",
        34: "weber",
        35: "decibel",
        36: "dalton",
        37: "thompson",
        38: "ugeneric",
        39: "lpercent",
        40: "ppt",
        41: "ppb",
        42: "index",
    }

    data_axis_ranged = {
        0: "ranged",
        1: "listed",
        2: "sparse",
        3: "listed",
    }

    value_type = {
        0: "string",
        1: "integer",
        2: "float",
        3: "complex",
        4: "infinity",
    }

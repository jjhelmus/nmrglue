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
    with open(fname, "rb") as f:
        data = f.read()

    header = parse_jeol(data)

    return header


class IOBuffer:
    def __init__(self, buffer):
        self.buffer = buffer
        self.position = 0

    def set_big_endian(self):
        pass

    def read_chars(self, count):
        chars = struct.unpack_from(f"{count}s", self.buffer, self.position)
        self.position += count
        return chars[0].decode("utf-8")

    def read_int8(self):
        value = struct.unpack_from(">b", self.buffer, self.position)
        self.position += 1
        return value[0]

    def read_uint8(self):
        value = struct.unpack_from(">B", self.buffer, self.position)
        self.position += 1
        return value[0]

    def read_uint16(self):
        value = struct.unpack_from(">H", self.buffer, self.position)
        self.position += 2
        return value[0]

    def read_byte(self):
        value = struct.unpack_from("B", self.buffer, self.position)
        self.position += 1
        return value[0]

    def get_array(self, read_func, count):
        return [getattr(self, read_func)() for _ in range(count)]

    def get_unit(self, size):
        unit = []
        for i in range(size):
            byte = self.read_byte()
            prefix = Table.prefix_table[byte >> 4]
            power = byte & 0b00001111
            base = Table.base_table[self.read_int8()]
            unit.append({'prefix': prefix, 'power': power, 'base': base})
        return unit

    def get_string(self, count):
        return self.read_chars(count).replace("\u0000", "")


def parse_jeol(buffer):
    io_buffer = IOBuffer(buffer)

    # read header section
    header = {}

    header["file_identifier"] = io_buffer.read_chars(8)
    header["endian"] = Table.endianness[io_buffer.read_int8()]
    header["major_version"] = io_buffer.read_uint8()
    header["minor_version"] = io_buffer.read_uint16()
    header["ndim"] = io_buffer.read_uint8()

    data_dimension_exist_byte = io_buffer.read_byte()
    header["dim_exist"] = [
        bool(int(x)) for x in format(data_dimension_exist_byte, "08b")
    ]

    info = io_buffer.read_byte()
    header["data_type"] = Table.data_type_table[info >> 6]
    header["data_format"] = Table.data_type_table[info & 0b00111111]

    header["data_instrument"] = Table.instruments[io_buffer.read_int8()]
    header["translate"] = [io_buffer.read_int8() for i in range(8)]
    header["axis_type"] = [Table.axis_type[i] for i in io_buffer.get_array("read_int8", 8)]
    header["data_units"] = io_buffer.get_unit(8)
    header["title"] = io_buffer.get_string(124)

    return header


# Conversion Table for Jeol Data Format
class Table:
    endianness = {
        0: "big_endian",
        1: "little_endian",
    }

    instruments = {
        0: "NONE",
        1: "GSX",
        2: "ALPHA",
        3: "ECLIPSE",
        4: "MASS_SPEC",
        5: "COMPILER",
        6: "OTHER_NMR",
        7: "UNKNOWN",
        8: "GEMINI",
        9: "UNITY",
        10: "ASPECT",
        11: "UX",
        12: "FELIX",
        13: "LAMBDA",
        14: "GE_1280",
        15: "GE_OMEGA",
        16: "CHEMAGNETICS",
        17: "CDFF",
        18: "GALACTIC",
        19: "TRIAD",
        20: "GENERIC_NMR",
        21: "GAMMA",
        22: "JCAMP_DX",
        23: "AMX",
        24: "DMX",
        25: "ECA",
        26: "ALICE",
        27: "NMR_PIPE",
        28: "SIMPSON",
    }

    data_type_table = {
        0: "float64",
        1: "float32",
        2: "Reserved",
        3: "Reserved",
    }

    data_format_table = {
        1: "One_D",
        2: "Two_D",
        3: "Three_D",
        4: "Four_D",
        5: "Five_D",
        6: "Six_D",
        7: "Seven_D",
        8: "Eight_D",
        9: "not for NMR data formats",
        10: "not for NMR data formats",
        11: "not for NMR data formats",
        12: "Small_Two_D",
        13: "Small_Three_D",
        14: "Small_Four_D",
    }

    axis_type = {
        0: "None",
        1: "Real",
        2: "TPPI",
        3: "Complex",
        4: "Real_Complex",
        5: "Envelope",
    }

    prefix_table = {
        "-8": "Yotta",
        "-6": "Exa",
        "-7": "Zetta",
        "-5": "Pecta",
        "-4": "Tera",
        "-3": "Giga",
        "-2": "Mega",
        "-1": "Kilo",
        0: "None",
        1: "Milli",
        2: "Micro",
        3: "Nano",
        4: "Pico",
        5: "Femto",
        6: "Atto",
        7: "Zepto",
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

    base_table = {
        0: "None",
        1: "Abundance",
        2: "Ampere",
        3: "Candela",
        4: "Celsius",
        5: "Coulomb",
        6: "Degree",
        7: "Electronvolt",
        8: "Farad",
        9: "Sievert",
        10: "Gram",
        11: "Gray",
        12: "Henry",
        13: "Hertz",
        14: "Kelvin",
        15: "Joule",
        16: "Liter",
        17: "Lumen",
        18: "Lux",
        19: "Meter",
        20: "Mole",
        21: "Newton",
        22: "Ohm",
        23: "Pascal",
        24: "Percent",
        25: "Point",
        26: "Ppm",
        27: "Radian",
        28: "Second",
        29: "Siemens",
        30: "Steradian",
        31: "Tesla",
        32: "Volt",
        33: "Watt",
        34: "Weber",
        35: "Decibel",
        36: "Dalton",
        37: "Thompson",
        38: "Ugeneric",
        39: "LPercent",
        40: "PPT",
        41: "PPB",
        42: "Index",
    }

    data_axis_ranged_table = {
        0: "Ranged",
        1: "Listed",
        2: "Sparse",
        3: "Listed",
    }

    value_type_table = {
        0: "String",
        1: "Integer",
        2: "Float",
        3: "Complex",
        4: "Infinity",
    }

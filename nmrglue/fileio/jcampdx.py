"""
Functions for reading 1D JCAMP-DX files.
"""

import os
import re
from warnings import warn

import numpy as np

from . import fileiobase

__developer_info__ = """
JCAMP-DX file format information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The format reference publications are available at:
<http://www.jcamp-dx.org/protocols.html>

Notes:
- Writing NMR data in JCAMP-DX format is not currently supported
- Multi-dimensional JCAMP-files are not currently supported
  (see <http://www.jcamp-dx.org/ndnmr-index.html#2dnmr%20testfiles>)

"""


def _getkey(keystr):
    '''
    Format key strings.
    From JCAMP-DX specs:
    "When LABELS are parsed, alphabetic
    characters are converted to upper case, and all spaces,
    dashes, slashes, and underlines are discarded. (XUNITS,
    xunits, xUNITS, and X-UNITS are equivalent.)"
    '''
    return (keystr.strip().upper().replace(" ", "")
            .replace("-", "").replace("_", "").replace("/", ""))


def _readrawdic(filename, read_err=None):
    '''
    Reads JCAMP-DX file to key-value dictionary, from which
    actual data is separated later.
    '''

    dic = {"_comments": []}  # create empty dictionary
    filein = open(filename, 'r', errors=read_err)

    currentkey = None
    currentvaluestrings = []

    for line in filein:

        # split comments
        commentsplit = line.split("$$", 1)
        actual = commentsplit[0].lstrip()
        if len(commentsplit) > 1:
            dic["_comments"].append(commentsplit[1])

        # continue with rest:
        if not actual:
            continue  # line had with nothing but comments

        # for multi-line data, linebreak must be restored if it has been
        # cut out with comments:
        if actual[-1] != "\n":
            actual = actual + "\n"

        # encountered new key:
        if actual[:2] == "##":

            # push previous key/value pair to dic
            # single value is continuous string including newlines
            # but there might be multiple values if the same key exist
            # multiple times thus values are collected to list
            if currentkey is not None and currentvaluestrings:
                key = _getkey(currentkey)
                value = "".join(currentvaluestrings)  # collapse
                if not value.strip():
                    warn("JCAMP-DX key without value:" + key)
                else:
                    try:
                        dic[key].append(value)
                    except KeyError:
                        dic[key] = [value]
                currentkey = None
                currentvaluestrings = []

            if actual[:5] == "##END":
                continue

            # try to split to key and value and check sanity
            keysplit = actual.split("=", 1)
            if len(keysplit) < 2:
                warn("Bad JCAMP-DX line, cant split key and value correctly:" +
                     line)
                continue
            keystr = keysplit[0][2:]  # remove "##" already here
            valuestr = keysplit[1]
            if not keystr:
                warn("Empty key in JCAMP-DX line:" + line)
                currentkey = None
                currentvaluestrings = []
                continue

            # split ok, init new key
            currentkey = keystr
            currentvaluestrings.append(valuestr)

        # line continues data of previous key, append to currentvaluestrings:
        else:
            if currentkey is None:
                warn("JCAMP-DX data line without associated key:" + line)
                continue

            currentvaluestrings.append(commentsplit[0])

    filein.close()

    return dic


###############################################################################
# digit dictionaries for pseudodigit parsing
_DIGITS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]
_SQZ_DIGITS = {"@": "0",
               "A": "1", "B": "2", "C": "3", "D": "4", "E": "5",
               "F": "6", "G": "7", "H": "8", "I": "9",
               "a": "-1", "b": "-2", "c": "-3", "d": "-4", "e": "-5",
               "f": "-6", "g": "-7", "h": "-8", "i": "-9"}
_DIF_DIGITS = {"%": "0",
               "J": "1", "K": "2", "L": "3", "M": "4", "N": "5",
               "O": "6", "P": "7", "Q": "8", "R": "9",
               "j": "-1", "k": "-2", "l": "-3", "m": "-4", "n": "-5",
               "o": "-6", "p": "-7", "q": "-8", "r": "-9"}
_DUP_DIGITS = {"S": "1", "T": "2", "U": "3", "V": "4", "W": "5",
               "X": "6", "Y": "7", "Z": "8", "s": "9"}
###############################################################################


def _detect_format(dataline):
    '''
    Detects and returns digit format:
    0  Normal
    1  Pseudodigits
    -1 Error
    '''

    # regexp to find & skip the first value of line, that never begins
    # with a pseudodigit in any format
    firstvalue_re = re.compile(
        "(\s)*([+-]?\d+\.?\d*|[+-]?\.\d+)([eE][+-]?\d+)?(\s)*")

    xy_re = re.compile('^[0-9\.]+,[ ]?[0-9\.]+')

    index = firstvalue_re.match(dataline).end()
    if index is None:
        return -1
    try:
        firstchar = dataline[index:index+1]
    except IndexError:
        return -1
    # detect the format from the first character of the second value in line
    if firstchar in _SQZ_DIGITS:
        return 1
    if firstchar in _DIF_DIGITS:
        return 1
    if firstchar in _SQZ_DIGITS:
        return 1

    if re.search(xy_re, dataline):
        return 2

    return 0


def _parse_affn_pac(datalines):
    ''' Parses datalines that do NOT contain any pseudodigits  '''

    # regexp explained:
    # -values may be delimited with whitespace, comma, or sign (+/-)
    # -may contain leading + or -
    # -base number may have decimal separator (.) or not
    # -if decimal separator is present, number can be given without leading
    #  zero (.1234) or decimals (123.)
    # -exponent (E/e) may be present
    value_re = re.compile("(\s|,)*([+-]?\d+\.?\d*|[+-]?\.\d+)([eE][+-]?\d+)?")

    data = []
    for dataline in datalines:
        linedata = []
        for match in value_re.finditer(dataline):
            base = match.group(2)
            exp = match.group(3)
            try:
                value = float(base + (exp if exp is not None else ""))
            except ValueError:
                warn("Data parsing failed at line:" + dataline)
                return None
            linedata.append(value)
        if len(linedata) > 1:
            data.extend(linedata[1:])  # ignore first column (X value)
    return data


def _append_value(data, value_to_append, isdif):
    '''
    Helper function for _finish_value: actual data push happens
    here based on isdif flag (direct value or difference from prev)
    '''
    if isdif:
        data.append(data[-1] + value_to_append)
    else:
        data.append(value_to_append)


def _finish_value(valuestr, currentmode, prev_value_to_append, data):
    '''
    Helper for _parse_pseudo:
    -Processes value in prev_value_to_append based on currentmode
    -Parses and returns next value_to_append for next round
    '''

    try:
        # squeeze format: number is added to data array as such
        if currentmode == 1:
            if prev_value_to_append is not None:
                _append_value(data, *prev_value_to_append)
            new_value_to_append = (float(valuestr), False)  # isdif=False
            return new_value_to_append, True
        # diff format: number is diff from previous entry
        elif currentmode == 2:
            if prev_value_to_append is not None:
                _append_value(data, *prev_value_to_append)
            new_value_to_append = (float(valuestr), True)  # isdif=True
            return new_value_to_append, True
        # duplicate format: push previous number (or diff) n times
        elif currentmode == 3:
            if prev_value_to_append is None:
                warn("Parse error: DUP entry without preceding value")
                return None, False
            dupcount = int(valuestr)
            for _i in range(dupcount):
                _append_value(data, *prev_value_to_append)
            new_value_to_append = None
            return new_value_to_append, True
        else:  # first value
            return None, True
    except ValueError:
        return None, False


def _parse_pseudo(datalines):
    ''' Parses datalines packed with pseudodigits  '''

    # regexp to find the first value of line, that never begins
    # with a pseudodigit (exponents are not allowed here)
    firstvalue_re = re.compile("(\s)*([+-]?\d+\.?\d*|[+-]?\.\d+)")

    data = []
    currentmode = 0
    valuestr = []
    skip_checkpoint = False

    # since the DUP mode entries may duplicate previous number n times,
    # we can't directly append newly read values to data. Instead, store
    # them here until next value is read:
    value_to_append = None

    for dataline in datalines:
        if not dataline:
            continue

        # ignore first value of line (X value)
        firstmatch = firstvalue_re.match(dataline)
        y_valuestring = dataline[firstmatch.end():]

        first_of_line = True

        # parse rest one char at a time
        for char in y_valuestring.strip():

            # char is digit = continues the current number
            if char in _DIGITS:
                valuestr.append(char)
                continue

            # char is pseudodigit = begin new number
            try:
                valuechar = _SQZ_DIGITS[char]
                newmode = 1
            except KeyError:
                try:
                    valuechar = _DIF_DIGITS[char]
                    newmode = 2
                except KeyError:
                    try:
                        valuechar = _DUP_DIGITS[char]
                        newmode = 3
                    except KeyError:
                        warn("Unknown pseudo-digit: " +
                             char + " at line: " + dataline)
                        return None

            # finish previous number
            valuestr = "".join(valuestr)

            # before updating value_to_append, store the mode of last value
            # actually appended to data.
            # this is needed for the DIF checkpoint removal below
            previous_is_dif = False
            if currentmode == 2:
                previous_is_dif = True
            elif currentmode == 3 and value_to_append[1]:
                previous_is_dif = True

            if not skip_checkpoint:
                # append number in value_to_append to data array if exists,
                # and update value_to_append
                value_to_append, success = _finish_value(valuestr,
                                                         currentmode,
                                                         value_to_append,
                                                         data)
                if not success:
                    warn("Data parsing failed at line:" + dataline)
                    return None

            # in DIF mode last of line is same than the first of next line
            # (= checkpoint). In such case raise a flag that this number
            # is to be skipped in _finish_value
            skip_checkpoint = first_of_line and previous_is_dif

            # init new number
            currentmode = newmode
            valuestr = [valuechar]
            first_of_line = False

    # read ended. finish last number:
    if not skip_checkpoint:
        valuestr = "".join(valuestr)
        value_to_append, success = _finish_value(valuestr,
                                                 currentmode,
                                                 value_to_append,
                                                 data)
        if not success:
            warn("Data parsing failed at last dataline")
            return None

    # append last number now in value_to_append:
    if value_to_append is not None:
        _append_value(data, *value_to_append)

    return data


def _parse_xy_xy(datalines):
    pts = []
    len_group_data = 0
    for dataline in datalines:
        if not dataline:
            continue
        xy_re = re.compile('[^ ][0-9\.]+, [0-9\.]+')
        group_data = re.findall(xy_re, dataline)
        len_group_data = len(group_data)
        if len_group_data == 0:
            xy_re = re.compile('[^ ][0-9\.]+,[0-9\.]+;')
            group_data = re.findall(xy_re, dataline)

        for data in group_data:
            clean_data = data.replace(', ', ',')
            clean_data = clean_data.replace(';', '')
            x, y = clean_data.split(',')
            pts.append([float(x), float(y)])
    return [pts]


def _parse_data(datastring):
    '''
    Creates numpy array from datalines
    '''
    probe_data = datastring[80:320]
    if ',' in probe_data and not('.' in probe_data): # fix comma as decimal points
        datastring = datastring.replace(',', '.')

    datalines = datastring.split("\n")
    headerline = datalines[0]

    datatype = "R"
    if "I..I" in headerline:
        datatype = "I"

    datalines = datalines[1:]  # get rid of the header line (e.g. (X++(Y..Y)))
    mode = _detect_format(datalines[0])
    if mode == 1:
        data = _parse_pseudo(datalines)
    elif mode == 0:
        data = _parse_affn_pac(datalines)
    elif mode == 2:
        if headerline == '(X++(Y..Y))':
            data = _parse_affn_pac(datalines)
        else:
            data = _parse_xy_xy(datalines)
    else:
        return None
    if data is None:
        return None
    return np.asarray(data, dtype="float64"), datatype


def get_is_ntuples(dic):
    '''
    Determine data class from dic: XYDATA or NTUPLES
    '''
    is_ntuples = False  # default is XYDATA
    try:
        dataclass = dic["DATACLASS"][0]
        if dataclass.strip() == "NTUPLES":
            is_ntuples = True
    except KeyError:
        pass
    return is_ntuples


def find_yfactors(dic):
    '''
    Helper to find yfactors from NTUPLES format.
    Returns YFactors in tuple with order (R,I)
    '''

    # determine data class:
    is_ntuples = get_is_ntuples(dic)

    if is_ntuples:
        # first check which column is R and I:
        index_r = None
        index_i = None
        try:
            symbols = dic["SYMBOL"][0].split(",")
            symbols = [s.strip() for s in symbols]
            index_r = symbols.index("R")
            index_i = symbols.index("I")
        except (KeyError, IndexError, ValueError):
            return (None, None)

        try:
            factors = dic["FACTOR"][0].split(",")
            factors = [s.strip() for s in factors]
            factor_r = float(factors[index_r])
            factor_i = float(factors[index_i])
        except (KeyError, IndexError, ValueError):
            return (None, None)

    return (factor_r, factor_i)


def _getdataarray(dic, show_all_data=False):
    '''
    Main function for data array parsing, input is the
    raw dictionary from _readrawdic
    '''

    data = None

    is_ntuples = get_is_ntuples(dic)

    if is_ntuples:  # NTUPLES
        valuelist = None
        try:
            valuelist = dic["DATATABLE"]
        except KeyError:
            is_ntuples = False
            warn("NTUPLES without DATA TABLEs. Trying XYDATA instead...")

        if valuelist:
            rdatalist = []
            idatalist = []
            for value in valuelist:
                parseret = _parse_data(value)
                if parseret is None:
                    return None
                data, datatype = parseret
                if datatype == "I":
                    idatalist.append(data)
                else:
                    rdatalist.append(data)

            if show_all_data:
                data = { 'real': rdatalist, 'imaginary': idatalist }
            else:
                if len(rdatalist) > 1:
                    warn("NTUPLES: multiple real arrays, returning first one only")
                if len(idatalist) > 1:
                    warn("NTUPLES: multiple imaginary arrays, \
                         returning first one only")
                if rdatalist:
                    if idatalist:
                        data = [rdatalist[0], idatalist[0]]
                    else:
                        data = rdatalist[0]
                else:
                    if idatalist:
                        data = [None, idatalist[0]]

    if data is None:  # XYDATA
        try:
            valuelist = dic["XYDATA"]
            if len(valuelist) == 1:
                data, datatype = _parse_data(valuelist[0])
            else:
                warn("Multiple XYDATA arrays in JCAMP-DX file, \
                     returning first one only")
        except KeyError:
            warn("XYDATA not found ")

    if data is None:  # PEAK TABLE
        try:
            valuelist = dic["PEAKTABLE"]
            if len(valuelist) == 1:
                data, datatype = _parse_data(valuelist[0])
            else:
                warn("Multiple PEAKTABLE arrays in JCAMP-DX file, \
                     returning first one only")
        except KeyError:
            warn("PEAKTABLE not found ")

    # apply YFACTOR to data if available
    if is_ntuples:
        yfactor_r, yfactor_i = find_yfactors(dic)
        if yfactor_r is None or yfactor_r is None:
            warn("NTUPLES: YFACTORs not applied, parsing failed")
        elif show_all_data:
            for i, _ in enumerate(data['real']):
                data['real'][i] = data['real'][i] * yfactor_r
            for i, _ in enumerate(data['imaginary']):
                data['imaginary'][i] = data['imaginary'][i] * yfactor_i
        else:
            data[0] = data[0] * yfactor_r
            data[1] = data[1] * yfactor_i
    else:
        try:
            yfactor = float(dic["YFACTOR"][0])
            data = data * yfactor
        except (ValueError, IndexError):
            warn("YFACTOR not applied, parsing failed")
        except KeyError:
            pass

    return data


def read(filename, show_all_data=False, read_err=None):
    """
    Read JCAMP-DX file

    Parameters
    ----------
    filename : str
        File to read from.

    Returns
    -------
    dic : dict
        Dictionary of parameters.
    data : ndarray
        Array of NMR data, or a list NMR data arrays in order [real, imaginary]
    """

    if os.path.isfile(filename) is not True:
        raise IOError("file %s does not exist" % (filename))

    # first read everything (including data array) to "raw" dictionary,
    # in which data values are read as raw strings including whitespace
    # and newlines
    dic = _readrawdic(filename, read_err)

    # find and parse NMR data array from raw dic
    data = _getdataarray(dic, show_all_data)

    # remove data tables from dic
    try:
        dic['XYDATA_OLD'] = dic["XYDATA"]
        del dic["XYDATA"]
    except KeyError:
        pass
    try:
        del dic["DATATABLE"]
    except KeyError:
        pass

    # clean dic values from leading and trailing whitespace
    for key, valuelist in dic.items():
        dic[key] = [value.strip() for value in valuelist]

    return dic, data


def _find_firstx_lastx(dic):
    '''
    Helper for guess_udic: seeks firstx and lastx for
    sweep calculation
    '''

    firstx = None
    lastx = None

    # determine data class:
    is_ntuples = get_is_ntuples(dic)

    if is_ntuples:
        # first check which column is X:
        index_x = None
        try:
            symbols = dic["SYMBOL"][0].split(",")
            symbols = [s.strip() for s in symbols]
            index_x = symbols.index("X")
        except (KeyError, IndexError, ValueError):
            warn("Cannot found X column on NTUPLES")
        if index_x is not None:
            try:
                firsts = dic["FIRST"][0].split(",")
                firsts = [s.strip() for s in firsts]
                firstx = float(firsts[index_x])
            except (KeyError, IndexError, ValueError):
                warn("Cannot parse FIRST (X) on NTUPLES")
            try:
                lasts = dic["LAST"][0].split(",")
                lasts = [s.strip() for s in lasts]
                lastx = float(lasts[index_x])
            except (KeyError, IndexError, ValueError):
                warn("Cannot parse LAST (X) on NTUPLES")

    # XYDATA (try always if not yet found)
    if firstx is None and lastx is None:
        try:
            firstx = float(dic["FIRSTX"][0])
        except ValueError:
            warn('Cannot parse "FIRSTX"')
        except KeyError:
            warn('No "FIRSTX" in file')
        try:
            lastx = float(dic["LASTX"][0])
        except ValueError:
            warn('Cannot parse "LASTX"')
        except KeyError:
            warn('No "LASTX" in file')

    return firstx, lastx


def guess_udic(dic, data):
    """
    Guess parameters of universal dictionary from dic, data pair.

    Parameters
    ----------
    dic : dict
        Dictionary of JCAMP-DX parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    udic : dict
        Universal dictionary of spectral parameters.
    """

    # create an empty universal dictionary
    udic = fileiobase.create_blank_udic(1)

    # update default values (currently only 1D possible)
    # "label"
    try:
        label_value = dic[".OBSERVENUCLEUS"][0].replace("^", "")
        udic[0]["label"] = label_value
    except KeyError:
        # sometimes INSTRUMENTAL PARAMETERS is used:
        try:
            label_value = dic["INSTRUMENTALPARAMETERS"][0].replace("^", "")
            udic[0]["label"] = label_value
        except KeyError:
            pass

    # "obs"
    try:
        obs_value = float(dic[".OBSERVEFREQUENCY"][0])
        udic[0]["obs"] = obs_value
    except ValueError:
        warn('Cannot parse ".OBSERVE FREQUENCY"')
    except KeyError:
        pass

    # "size"
    if isinstance(data, list):
        data = data[0]  # if list [R,I]
    if data is not None:
        udic[0]["size"] = len(data)
    else:
        warn('No data, cannot set udic size')

    # "sw"
    # by format specs, FIRSTX and LASTX should be always in Hz
    firstx, lastx = _find_firstx_lastx(dic)

    if firstx is not None and lastx is not None:
        udic[0]["sw"] = abs(lastx - firstx)
    else:
        warn('Cannot set udic sweep')

    # keys not found in standard&required JCAMP-DX keys and thus left default:
    # car, complex, encoding

    return udic

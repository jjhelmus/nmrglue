"""
Functions for reading JCAMP-DX files.

The interface is oriented towards 1D data. Multidimensional NTUPLES spectra
can be read, as the set of 1D pages they are stored as, but guess_udic
describes their direct dimension only.
"""

import os
import re
from warnings import warn

import numpy as np

from . import fileiobase

__developer_info__ = """
JCAMP-DX file format information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The format reference publications are available at:
http://www.jcamp-dx.org/protocols.html

Notes:
    * Writing NMR data in JCAMP-DX format is not currently supported.
    * Multi-dimensional JCAMP-files are not currently supported.
      See http://www.jcamp-dx.org/ndnmr-index.html#2dnmr%20testfiles

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


def _parsejcampdx(filename):
    '''
    Actual JCAMP-DX reading. Returns a list of data sections i.e. "blocks",
    each of them being a dictionary of JCAMP-DX tags
    '''

    # "blocks" in JCAMP-DX may be nested, and they begin with ##TITLE and
    # end with ##END. Keep active blocks in stack:
    blockstack = []
    # keep reference to active block, i.e. the one to which to read
    activeblock = None
    # when encountering ##END, push the ready dict to another list
    readyblocklist = []

    filein = open(filename, 'r', encoding="utf-8-sig", errors="replace")

    currentkey = None
    currentvaluestrings = []

    for line in filein:

        # split comments
        commentsplit = line.split("$$", 1)
        actual = commentsplit[0].lstrip()
        if len(commentsplit) > 1:
            if activeblock:
                activeblock["_comments"].append(commentsplit[1])

        # continue with rest:
        if not actual:
            continue  # line had with nothing but comments

        # for multi-line data, linebreak must be restored if it has been
        # cut out with comments:
        if actual[-1] != "\n":
            actual += "\n"

        # encountered new key:
        if actual[:2] == "##":

            # push previous key/value pair to active dictionary
            # single value is continuous string including newlines
            # but there might be multiple values if the same key exist
            # multiple times thus values are collected to list
            if currentkey is not None and currentvaluestrings:
                key = _getkey(currentkey)
                value = "".join(currentvaluestrings)  # collapse
                if not value.strip():
                    warn(f"JCAMP-DX key without value: {key}")
                else:
                    try:
                        activeblock[key].append(value)
                    except KeyError:
                        activeblock[key] = [value]
                currentkey = None
                currentvaluestrings = []

            if actual[:7] == "##TITLE":
                # begin new block dictionary
                activeblock = {"_comments": []}
                blockstack.append(activeblock)

            if actual[:5] == "##END":
                # finalize current block
                if activeblock:  # ensure that we had active block instead of too many ##ENDs
                    readyblocklist.append(blockstack.pop())
                if blockstack:
                    # continue reading to previous block in stack, if available
                    activeblock = blockstack[-1]
                else:
                    # otherwise mark it None to prevent reading
                    activeblock = None
                continue

            if activeblock:
                # try to split to key and value and check sanity
                keysplit = actual.split("=", 1)
                if len(keysplit) < 2:
                    warn("Bad JCAMP-DX line, can't split key and value correctly:" +
                         line)
                    continue
                keystr = keysplit[0][2:]  # remove "##" already here
                valuestr = keysplit[1]
                if not keystr:
                    warn(f"Empty key in JCAMP-DX line: {line}")
                    currentkey = None
                    currentvaluestrings = []
                    continue

                # split ok, init new key
                currentkey = keystr
                currentvaluestrings.append(valuestr)

        # line continues data of previous key, append to currentvaluestrings:
        else:
            if activeblock:
                if currentkey is None:
                    warn(f"JCAMP-DX data line without associated key: {line}")
                    continue

                currentvaluestrings.append(commentsplit[0])

    # push possible non-closed blocks
    while blockstack:
        readyblocklist.append(blockstack.pop())

    filein.close()

    return readyblocklist


def _readrawdic(filename):
    '''
    Reads entire JCAMP-DX file to dictionary, from which actual
    data is parsed later. Return value is a dictionary of different
    DATATYPEs, containing lists of actual block dictionaries(as it
    is possible have multiple entries of same DATATYPE in one file).
    '''

    # parse file to list of "blocks" i.e. separate data sections
    blocklist = _parsejcampdx(filename)

    # clean whitespace from entries, and remove empty entries
    cleandiclist = []
    for dic in blocklist:
        for key, valuelist in dic.items():
            dic[key] = [value.strip() for value in valuelist]
        for key, valuelist in dic.items():
            dic[key] = [value for value in valuelist if value]
        dic = {key: valuelist for key, valuelist in dic.items() if valuelist}
        if dic:
            cleandiclist.append(dic)

    returndic = {}
    # check DATATYPE entry of each block,
    # and build a dict of lists of dicts
    for dic in cleandiclist:
        try:
            datatypelist = dic["DATATYPE"]
            currdatatype = datatypelist[0].strip().upper().replace(" ", "")
            if len(datatypelist) > 1:
                # basically sections with multiple DATATYPES
                # are invalid, but we may still give it a try:
                for datatype in datatypelist:
                    cleandatatype = datatype.strip().upper().replace(" ", "")
                    if cleandatatype == "NMRSPECTRUM":
                        currdatatype = cleandatatype
                        break
                    if cleandatatype == "NMRFID":
                        currdatatype = cleandatatype
                        break
                # no SPECTRUM / FID found, just go with the first entry
        except KeyError:
            # no datatype in this section, use dummy
            currdatatype = "NA"

        # push to result dict
        key = "_datatype_"+currdatatype
        try:
            returndic[key].append(dic)
        except KeyError:
            returndic[key] = [dic]

    return returndic


###############################################################################
# digit dictionaries for pseudodigit parsing
_DIGITS = ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "."]
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
        r"(\s)*([+-]?\d+\.?\d*|[+-]?\.\d+)([eE][+-]?\d+)?(\s)*")

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
    if firstchar in _DUP_DIGITS:
        return 1
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
    value_re = re.compile(r"(\s|,)*([+-]?\d+\.?\d*|[+-]?\.\d+)([eE][+-]?\d+)?")

    data = []
    for dataline in datalines:
        linedata = []
        for match in value_re.finditer(dataline):
            base = match.group(2)
            exp = match.group(3)
            try:
                value = float(base + (exp if exp is not None else ""))
            except ValueError:
                warn(f"Data parsing failed at line: {dataline}")
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
    firstvalue_re = re.compile(r"(\s)*([+-]?\d+\.?\d*|[+-]?\.\d+)")

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
                        warn(f"Unknown pseudo-digit: {char} at line: {dataline}")
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
                    warn(f"Data parsing failed at line: {dataline}")
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


# variable list of a data table header, e.g. (X++(R..R)) or (F2++(Y..Y)):
# the first symbol is the abscissa, the second the dependent variable
_VARIABLE_LIST_RE = re.compile(
    r"\(\s*([A-Za-z][A-Za-z0-9]*)\s*\+\+\s*\(\s*([A-Za-z][A-Za-z0-9]*)\s*\.\.")


def _parse_variable_list(headerline):
    '''
    Finds the variable symbols declared by a data table header, e.g.
    "(X++(R..R))" gives ("X", "R") and "(F2++(Y..Y))" gives ("F2", "Y").
    Returns (None, None) when the header declares no variable list.
    '''
    match = _VARIABLE_LIST_RE.search(headerline)
    if match is None:
        return (None, None)
    return match.group(1), match.group(2)


def _parse_data(datastring):
    '''
    Creates numpy array from datalines
    '''
    datalines = datastring.split("\n")
    headerline = datalines[0]

    # the dependent variable is named by the header's variable list; headers
    # that declare none are assumed to hold the real component
    datatype = _parse_variable_list(headerline)[1]
    if datatype is None:
        datatype = "R"

    datalines = datalines[1:]  # get rid of the header line (e.g. (X++(Y..Y)))
    mode = _detect_format(datalines[0])
    if mode == 1:
        data = _parse_pseudo(datalines)
    elif mode == 0:
        data = _parse_affn_pac(datalines)
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


def find_factor(dic, symbol):
    '''
    Helper to find the scaling factor of the given symbol (e.g. "R", "I"
    or "Y") from the NTUPLES ##SYMBOL= and ##FACTOR= records, whose
    comma-separated entries correspond by position.
    Returns None if the factor cannot be determined.
    '''
    try:
        symbols = [s.strip() for s in dic["SYMBOL"][0].split(",")]
        factors = [s.strip() for s in dic["FACTOR"][0].split(",")]
        return float(factors[symbols.index(symbol)])
    except (KeyError, IndexError, ValueError):
        return None


def find_yfactors(dic):
    '''
    Helper to find yfactors from NTUPLES format.
    Returns YFactors in tuple with order (R,I)
    '''
    return (find_factor(dic, "R"), find_factor(dic, "I"))


def getdataarray(dic, show_all_data=False):
    '''
    Main function for data array parsing, input is the
    raw dictionary from _readrawdic

    Parameters
    ----------
    dic : dict
        Raw dictionary from _readrawdic.
    show_all_data : bool
        If True and data is NTUPLES, return all data arrays as a dict
        with keys 'real' and 'imaginary', each containing a list of arrays.
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
                # scale by the factor of this table's own dependent
                # symbol: an (X++(I..I)) table gets I's factor, an
                # (X++(R..R)) table R's, so factors can never be
                # applied to the wrong component
                factor = find_factor(dic, datatype)
                if factor is None:
                    warn("NTUPLES: no FACTOR found for symbol %s, "
                         "data not scaled" % datatype)
                else:
                    data = data * factor
                if datatype == "I":
                    idatalist.append(data)
                else:
                    rdatalist.append(data)
            if show_all_data:
                data = {'real': rdatalist, 'imaginary': idatalist}
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
            if len(valuelist) > 1:
                warn("Multiple XYDATA arrays in JCAMP-DX file, \
                     returning first one only")
            parseret = _parse_data(valuelist[0])
            if parseret is None:
                return None
            data, datatype = parseret
        except KeyError:
            warn("XYDATA not found ")

    if data is None:
        return None

    # apply YFACTOR to data if available
    # (NTUPLES data was already scaled per-table above)
    if not is_ntuples:
        try:
            yfactor = float(dic["YFACTOR"][0])
            data = data * yfactor
        except (ValueError, IndexError):
            warn("YFACTOR not applied, parsing failed")
        except KeyError:
            pass

    return data


def read(filename, show_all_data=False):
    """
    Read JCAMP-DX file

    Parameters
    ----------
    filename : str
        File to read from.
    show_all_data : bool
        If True and data is NTUPLES, return all data arrays as a dict
        with keys 'real' and 'imaginary', each containing a list of
        numpy arrays. If False (default), return only the first real
        and imaginary arrays.

    Returns
    -------
    dic : dict
        Dictionary of parameters. In the case of multiple data sections in
        file, parameters of first NMR SPECTRUM or NMR FID are read to base
        level and others are stored under _datatype_<DATATYPE> keys in the
        dictionary.
    data : ndarray or dict
        Array of NMR data, or a list of NMR data arrays in order
        [real, imaginary]. When show_all_data=True and data is NTUPLES,
        a dict with keys 'real' and 'imaginary' is returned.
    """

    if os.path.isfile(filename) is not True:
        raise OSError("file %s does not exist" % (filename))

    # first read everything (including data array) to "raw" dictionary,
    # in which data values are read as raw strings including whitespace
    # and newlines
    dic = _readrawdic(filename)

    # select the relevant data section, taking the first section of the most
    # preferred DATATYPE that yields data. Non-typed sections are tried
    # because the DATATYPE label is sometimes missing, and multidimensional
    # spectra come last so that a file offering both a 1D and an nD spectrum
    # keeps returning the 1D one.
    data = None
    correctdic = None
    for datatype in ("NMRSPECTRUM", "NMRFID", "NA", "NDNMRSPECTRUM"):
        for subdic in dic.get("_datatype_" + datatype, []):
            data = getdataarray(subdic, show_all_data)
            if data is not None:
                correctdic = subdic
                break
        if data is not None:
            break

    if data is None:
        warn("no data found either in XYDATA or NTUPLES format")

    if correctdic is not None:
        # remove correct dic from subdic lists:
        for key, subdiclist in dic.items():
            for subdic in subdiclist:
                if subdic is correctdic:
                    subdiclist = [d for d in subdiclist if d is not correctdic]
                    dic[key] = subdiclist

        # clean correct dic:
        # remove data tables
        try:
            del correctdic["XYDATA"]
        except KeyError:
            pass
        try:
            del correctdic["DATATABLE"]
        except KeyError:
            pass

        # push correct dic entries to base level of main dic
        for key, valuelist in correctdic.items():
            dic[key] = valuelist

    # clean main dic from possible empty entries
    dic = {key: value for key, value in dic.items() if value}

    return dic, data


def _find_abscissa_symbol(dic, symbols):
    '''
    Finds which of an NTUPLES variable list's symbols is the abscissa.
    One dimensional data calls it X, while multidimensional data names its
    dimensions instead. In that case the innermost, i.e. last, independent
    variable is the abscissa of every page: a file declaring (F1, F2, Y)
    stores (F2++(Y..Y)) data tables, one per value of F1.
    Returns None if it cannot be identified.
    '''
    if "X" in symbols:
        return "X"

    try:
        vartypes = [s.strip().upper() for s in dic["VARTYPE"][0].split(",")]
        independent = [symbol for symbol, vartype in zip(symbols, vartypes)
                       if vartype == "INDEPENDENT"]
        if independent:
            return independent[-1]
    except (KeyError, IndexError):
        pass

    # no VARTYPE to go on: fall back to the abscissa a data table declares,
    # available while reading but not after read() strips the tables
    try:
        return _parse_variable_list(dic["DATATABLE"][0].split("\n", 1)[0])[0]
    except (KeyError, IndexError):
        return None


def _find_firstx_lastx(dic):
    '''
    Helper for guess_udic: seeks firstx and lastx for
    sweep calculation. Also returns True/False if the
    data was in ppm.
    '''

    firstx = None
    lastx = None
    unitx = None
    isppm = False  # default to Hz

    # determine data class:
    is_ntuples = get_is_ntuples(dic)

    if is_ntuples:
        # first check which column is the abscissa:
        index_x = None
        try:
            symbols = dic["SYMBOL"][0].split(",")
            symbols = [s.strip() for s in symbols]
            index_x = symbols.index(_find_abscissa_symbol(dic, symbols))
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
            try:
                units = dic["UNITS"][0].split(",")
                units = [s.strip() for s in units]
                unitx = units[index_x]
            except (KeyError, IndexError, ValueError):
                warn("Cannot parse UNITS (X) on NTUPLES")

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
    if unitx is None:
        try:
            unitx = dic["XUNITS"][0]
        except KeyError:
            warn('No "XUNITS" in file')

    # flag ppm data
    if unitx is not None:
        isppm = unitx.upper() == "PPM"

    return firstx, lastx, isppm


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

    # the universal dictionary is one dimensional, so for multidimensional
    # data it can only describe the direct dimension
    try:
        if int(dic["NUMDIM"][0]) > 1:
            warn("Multidimensional data: udic describes the direct "
                 "dimension only")
    except (KeyError, IndexError, ValueError):
        pass

    # update default values
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
    obs_freq = None
    try:
        obs_freq = float(dic[".OBSERVEFREQUENCY"][0])
        udic[0]["obs"] = obs_freq
    except ValueError:
        warn('Cannot parse ".OBSERVE FREQUENCY"')
    except KeyError:
        pass

    # "size"
    if isinstance(data, dict):
        # show_all_data form: every page has the same length, so measure the
        # first one that is present
        pages = data.get("real") or data.get("imaginary") or [None]
        data = pages[0]
    if isinstance(data, list):
        data = data[0]  # if list [R,I]
    if data is not None:
        udic[0]["size"] = len(data)
    else:
        warn('No data, cannot set udic size')

    # "sw"
    # get firstx, lastx and unit
    firstx, lastx, isppm = _find_firstx_lastx(dic)

    # ppm data: convert to Hz
    if isppm:
        if obs_freq:
            firstx = firstx * obs_freq
            lastx = lastx * obs_freq
        else:
            firstx, lastx = (None, None)
            warn('Data is in ppm but have no frequency, cannot set udic sweep')

    if firstx is not None and lastx is not None:
        udic[0]["sw"] = abs(lastx - firstx)
    else:
        warn('Cannot set udic sweep')

    # keys not found in standard&required JCAMP-DX keys and thus left default:
    # car, complex, encoding

    return udic

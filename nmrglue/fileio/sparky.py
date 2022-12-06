"""
Functions for reading and writing Sparky (.ucsf) files.
"""


__developer_info__ = """
Sparky file format information
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Information on the Sparky file format can be found online at:
http://www.cgl.ucsf.edu/home/sparky/manual/files.html
and in the Sparky source file ucsffile.cc.

"""

import os
import struct
import datetime
from warnings import warn
try:
    from html.parser import HTMLParser
except ImportError:
    from HTMLParser import HTMLParser

import numpy as np

from . import fileiobase


# unit conversion function
def make_uc(dic, data, dim=-1):
    """
    Create a unit conversion object.

    Parameters
    ----------
    dic : dict
        Dictionary of Sparky parameters.
    data : ndarray
        Array of NMR data.
    dim : int, optional
        Dimension number to create unit conversion object for.  Default is for
        last dimension.

    Returns
    -------
    uc : unit conversion object.
        Unit conversion object for given dimension.

    """
    if dim == -1:
        dim = data.ndim - 1  # last dimension

    wdic = dic["w" + str(int(1 + dim))]

    size = float(wdic["npoints"])
    cplx = False
    sw = wdic["spectral_width"]
    obs = wdic["spectrometer_freq"]
    car = wdic["xmtr_freq"] * obs

    return fileiobase.unit_conversion(size, cplx, sw, obs, car)


# dictionary/data creation
def create_data(data):
    """
    Create a Sparky data array (recast into float32 array)
    """
    return np.array(data, dtype="float32")


def guess_udic(dic, data):
    """
    Guess parameter of universal dictionary from dic,data pair.

    Parameters
    ----------
    dic : dict
        Dictionary of Sparky parameters.
    data : ndarray
        Array of NMR data.

    Returns
    -------
    udic : dict
        Universal dictionary of spectral parameter.

    """
    # create an empty universal dictionary
    udic = fileiobase.create_blank_udic(data.ndim)

    # update default values
    for i in range(data.ndim):
        adic = dic["w" + str(i + 1)]
        udic[i]["size"] = data.shape[i]
        udic[i]["sw"] = adic['spectral_width']
        udic[i]["obs"] = adic['spectrometer_freq']
        udic[i]["car"] = adic['xmtr_freq'] * adic['spectrometer_freq']
        udic[i]["label"] = adic['nucleus']
        udic[i]["complex"] = False
        udic[i]["time"] = False
        udic[i]["freq"] = True

    return udic


def create_dic(udic, datetimeobj=datetime.datetime.now(), user='user'):
    """
    Create a Sparky parameter dictionary from universal dictionary.

    Parameters
    ----------
    udic : dict
        Universal dictionary of spectral parameters.
    datatimeobj : datetime object, optional
        Datetime to record in Sparky dictionary
    user : str, optional
        Username to record in Sparky dictionary. Default is 'user'

    Returns
    -------
    dic : dict
        Dictionary of Sparky parameters.

    """
    dic = dict()

    # determine shape of array
    shape = [udic[k]["size"] for k in range(udic["ndim"])]

    # populate the dictionary
    dic["ident"] = 'UCSF NMR'
    dic["naxis"] = udic["ndim"]
    dic["ncomponents"] = 1
    dic["encoding"] = 0
    dic["version"] = 2
    dic["owner"] = user
    dic["date"] = datetimeobj.ctime()
    dic["comment"] = ''
    dic["scratch"] = ''

    # calc a good tile shape
    tshape = calc_tshape(shape)

    # total number of tiles
    ntiles = 1
    for tlen, slen in zip(tshape, shape):
        ntiles *= np.ceil(float(slen) / tlen)

    # points in tile
    tpoints = np.array(tshape).prod()

    # data bytes
    dbytes = tpoints * ntiles * 4

    # total file size if data size plus leaders
    dic["seek_pos"] = int(dbytes + 180 + 128 * len(shape))

    # populate the dictionary with axis dictionaries
    for i, (tlen, dlen) in enumerate(zip(tshape, shape)):
        dic["w" + str(i + 1)] = create_axisdic(udic[i], tlen, dlen)
    return dic


def create_axisdic(adic, tlen, dlen):
    """
    Make an Sparky axis dictionary from a universal axis dictionary.

    Parameters
    ----------
    adic : dict
        Axis dictionary from a universal dictionary.
    tlen : int
        Tile length of axis.
    dlen : int
        Data length of axis.

    Returns
    -------
    sdic : dict
        Sparky axis dictionary

    """
    dic = dict()
    dic["nucleus"] = adic["label"]
    dic["spectral_shift"] = 0
    dic["npoints"] = int(dlen)
    dic["size"] = int(dlen)
    dic["bsize"] = int(tlen)
    dic["spectrometer_freq"] = float(adic["obs"])
    dic["spectral_width"] = float(adic["sw"])
    dic["xmtr_freq"] = float(adic["car"]) / dic["spectrometer_freq"]
    dic["zero_order"] = 0.0
    dic["first_order"] = 0.0
    dic["first_pt_scale"] = 0.0
    dic["extended"] = b'\x80'  # transform bit set
    return dic


def datetime2dic(datetimeobj, dic):
    """
    Add datetime object to dictionary
    """
    dic["date"] = datetimeobj.ctime()
    return dic


def dic2datetime(dic):
    """
    Create a datetime object from a Sparky dictionary
    """
    return datetime.datetime.strptime(dic["date"], "%a %b %d %H:%M:%S %Y")


def calc_tshape(shape, kbyte_max=128):
    """
    Calculate a tile shape from data shape.


    Parameters
    ----------
    shape : tuple
        Shape of NMR data (data.shape).
    kbyte_max : float or int
        Maximum tile size in Kilobytes.

    Returns
    -------
    tshape : tuple
        Shape of tile.

    """
    # Algorithm divides each dimension by 2 until under kbyte_max tile size.
    s = np.array(shape, dtype="int")
    i = 0
    while (s.prod() * 4. / 1024. > kbyte_max):
        s[i] = np.floor(s[i] / 2.)
        i = i + 1
        if i == len(s):
            i = 0
    return tuple(s)


# global read/write functions
def read(filename):
    """
    Read a Sparky file.

    Parameters
    ----------
    filename : str
        Filename of Sparky file to read.

    Returns
    -------
    dic : dict
        Dictionary of Sparky parameters.
    data : ndarray
        Array of NMR data.

    See Also
    --------
    read_lowmem : Sparky file reading with minimal memory usage.
    write : Write a Sparky file.

    """
    # open the file
    f = open(filename, 'rb')

    # determine the dimensionality
    n = fileheader2dic(get_fileheader(f))["naxis"]
    f.close()

    if n == 2:
        return read_2D(filename)
    if n == 3:
        return read_3D(filename)
    if n == 4:
        return read_4D(filename)

    raise ValueError("unknown dimensionality: %s" % n)


def read_lowmem(filename):
    """
    Read a Sparky file using minimal memory.

    Parameters
    ----------
    filename : str
        Filename of Sparky file to read.

    Returns
    -------
    dic : dict
        Dictionary of Sparky parameters.
    data : array_like
        Low memory object which can access NMR data on demand.

    See Also
    --------
    read : Read a Sparky file.
    write_lowmem : Write a Sparky file using minimal memory.

    """
    # open the file
    f = open(filename, 'rb')

    # determine the dimensionality
    n = fileheader2dic(get_fileheader(f))["naxis"]
    f.close()

    if n == 2:
        return read_lowmem_2D(filename)
    if n == 3:
        return read_lowmem_3D(filename)

    raise ValueError("unknown dimensionality: %s" % n)


def write(filename, dic, data, overwrite=False):
    """
    Write a Sparky file.

    Parameters
    ----------
    filename : str
        Filename of Sparky file to write to.
    dic : dict
        Dictionary of Sparky parameters.
    data : array_like
        Array of NMR data.
    overwrite : bool, optional
        Set True to overwrite files, False will raise a Warning if the file
        exists.

    See Also
    --------
    write_lowmem : Write a Sparky file using minimal amounts of memory.
    read : Read a Sparky file.

    """
    n = dic["naxis"]

    if n == 2:
        return write_2D(filename, dic, data, overwrite=overwrite)
    if n == 3:
        return write_3D(filename, dic, data, overwrite=overwrite)

    raise ValueError("unknown dimensionality: %s" % n)


def write_lowmem(filename, dic, data, overwrite=False):
    """
    Write a Sparky using minimum amounts of memory (tile by tile)

    Parameters
    ----------
    filename : str
        Filename of Sparky file to write to.
    dic : dict
        Dictionary of Sparky parameters.
    data : array_like.
        Array of NMR data.
    overwrite : bool, optional
        Set True to overwrite files, False will raise a Warning if the file
        exists.

    See Also
    --------
    write : Write a Sparky file.
    read_lowmem : Read a Sparky file using minimal amounts of memory.

    """
    # write also writes tile by tile...
    return write(filename, dic, data, overwrite)


# dimension specific reading/writing functions
def read_2D(filename):
    """
    Read a 2D sparky file. See :py:func:`read` for documentation.
    """
    seek_pos = os.stat(filename).st_size
    with open(filename, 'rb') as f:

        # read the file header
        dic = fileheader2dic(get_fileheader(f))

        # check for file size mismatch
        expected_seek_pos = dic["seek_pos"]
        if seek_pos != expected_seek_pos:
            warn(f'Bad file size in header {seek_pos} vs {expected_seek_pos}')

        # read the axis headers...
        for i in range(dic['naxis']):
            dic["w" + str(i + 1)] = axisheader2dic(get_axisheader(f))

        # read the data and untile
        lenY = dic["w1"]["npoints"]
        lenX = dic["w2"]["npoints"]
        lentY = dic["w1"]["bsize"]
        lentX = dic["w2"]["bsize"]
        data = get_data(f)
        data = untile_data2D(data, (lentY, lentX), (lenY, lenX))

        return dic, data


def write_2D(filename, dic, data, overwrite=False):
    """
    Write a 2D Sparky file. See :py:func:`write` for documentation.
    """
    # open the file for writing
    f = fileiobase.open_towrite(filename, overwrite)

    # write the file header
    put_fileheader(f, dic2fileheader(dic))

    # write the axis headers
    put_axisheader(f, dic2axisheader(dic["w1"]))
    put_axisheader(f, dic2axisheader(dic["w2"]))

    lentX = dic["w2"]["bsize"]
    lentY = dic["w1"]["bsize"]
    t_tup = (lentY, lentX)

    ttX = int(np.ceil(data.shape[1] / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(data.shape[0] / float(lentY)))  # total tiles in Y dim
    tt = ttX * ttY

    for i in range(int(tt)):
        put_data(f, find_tilen_2d(data, i, (t_tup)))

    f.close()
    return


def read_3D(filename):
    """
    Read a 3D Sparky file. See :py:func:`read` for documentation.
    """
    seek_pos = os.stat(filename).st_size
    with open(filename, 'rb') as f:

        # read the file header
        dic = fileheader2dic(get_fileheader(f))

        # check for file size mismatch
        expected_seek_pos = dic["seek_pos"]
        if seek_pos != expected_seek_pos:
            warn(f'Bad file size in header {seek_pos} vs {expected_seek_pos}')

        # read the axis headers...
        for i in range(dic['naxis']):
            dic["w" + str(i + 1)] = axisheader2dic(get_axisheader(f))

        # read the data and untile
        lenZ = dic["w1"]["npoints"]
        lenY = dic["w2"]["npoints"]
        lenX = dic["w3"]["npoints"]
        lentZ = dic["w1"]["bsize"]
        lentY = dic["w2"]["bsize"]
        lentX = dic["w3"]["bsize"]
        data = get_data(f)
        data = untile_data3D(data, (lentZ, lentY, lentX), (lenZ, lenY, lenX))

        return dic, data

def write_3D(filename, dic, data, overwrite=False):
    """
    Write a 3D Sparky file. See :py:func:`write` for documentation.
    """
    # open the file for writing
    f = fileiobase.open_towrite(filename, overwrite)

    # write the file header
    put_fileheader(f, dic2fileheader(dic))

    # write the axis headers
    put_axisheader(f, dic2axisheader(dic["w1"]))
    put_axisheader(f, dic2axisheader(dic["w2"]))
    put_axisheader(f, dic2axisheader(dic["w3"]))

    lentX = dic["w3"]["bsize"]
    lentY = dic["w2"]["bsize"]
    lentZ = dic["w1"]["bsize"]

    t_tup = (lentZ, lentY, lentX)

    ttX = int(np.ceil(data.shape[2] / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(data.shape[1] / float(lentY)))  # total tiles in Y dim
    ttZ = int(np.ceil(data.shape[0] / float(lentZ)))  # total tiles in Z dim

    tt = ttX * ttY * ttZ

    for i in range(int(tt)):
        put_data(f, find_tilen_3d(data, i, (t_tup)))
    f.close()
    return

def read_4D(filename):
    """
    Read a 4D Sparky file. See :py:func:`read` for documentation.
    """
    seek_pos = os.stat(filename).st_size

    with open(filename, 'rb') as f:

        # read the file header
        dic = fileheader2dic(get_fileheader(f))

        # check for file size mismatch
        expected_seek_pos = dic['seek_pos']
        if seek_pos != expected_seek_pos:
            warn(f'Bad file size in header {seek_pos} vs {expected_seek_pos}')

        # read the axis headers...
        for i in range(dic['naxis']):
            dic["w" + str(i + 1)] = axisheader2dic(get_axisheader(f))

        # read the data and untile
        lenA = dic["w1"]["npoints"]
        lenZ = dic["w2"]["npoints"]
        lenY = dic["w3"]["npoints"]
        lenX = dic["w4"]["npoints"]
        lentA = dic["w1"]["bsize"]
        lentZ = dic["w2"]["bsize"]
        lentY = dic["w3"]["bsize"]
        lentX = dic["w4"]["bsize"]
        data = get_data(f)

        data = untile_data4D(data, (lentA, lentZ, lentY, lentX), (lenA, lenZ, lenY, lenX))

        return dic, data

# TODO: read_lowmem_4D ?

# read_lowmem functions
def read_lowmem_2D(filename):
    """
    Read a 2D Sparky file using minimal memory. See :py:func:`read_lowmem`.
    """
    seek_pos = os.stat(filename).st_size

    # create the sparky_2d file
    data = sparky_2d(filename)
    dic = dict(data.dic)

    # check for file size mismatch
    expected_seek_pos = dic['seek_pos']
    if seek_pos != expected_seek_pos:
        warn(f'Bad file size in header {seek_pos} vs {expected_seek_pos}')
    return dic, data


def read_lowmem_3D(filename):
    """
    Read a 3D sparky file using minimal memory. See :py:func:`read_lowmem`.
    """
    seek_pos = os.stat(filename).st_size

    # create the sparky_3d file
    data = sparky_3d(filename)
    dic = dict(data.dic)

    # check for file size mismatch
    expected_seek_pos = dic['seek_pos']
    if seek_pos != expected_seek_pos:
        warn(f'Bad file size in header {seek_pos} vs {expected_seek_pos}')

    return dic, data


class SparkySaveParser(HTMLParser):
    """
    A parser for Sparky .save files. The file structure is similar to
    simple HTML files, except the use of <end tag> instead of </tag>.
    The following structure is assumed:

    <sparky save file>
    <version ...>
    <user>
    ...
    <end user>
    <spectrum>
    ...
        <view>
        ...
            <params>
            ...
            <end params>
            <params>
            ...
            <end params>
        <end view>
        <view>
        ...
            <params>
            ...
            <end params>
            <params>
            ...
            <end params>
        <end view>
        <ornament>
        ...
        <end ornament>
    <end spectrum>

    TODO: some .save files do not have this exact structure
    They need to be treated differently

    """

    # main dictionaries to parse data into
    user, spectrum, view, ornament = {}, {}, {}, {}

    # tracker if there are multiple views
    viewnum = -1

    curtag, curdict = None, None

    def _parse_info(self, string, dtype=None):
        """
        Reads a list of strings into a dictionary, with the first item of
        the list as the key and the remaining list as the value. In addition,
        it parses all values in the list to do the following: (i) convert the
        values to float wherever possible and (ii) if the list has a single
        item, upack and return that item alone as the value

        """

        dic = {}
        for s in string:
            i = s.split()
            try:
                if dtype is None:
                    key = i[0].replace(".", "_").replace("-", "_")
                    value = i[1:]

                if dtype == "user":
                    key = i[0] + "_" + i[1]
                    value = i[2:]

                parsed_value = []
                for v in value:
                    try:
                        if key == "id":
                            parsed_value.append(int(v))
                        else:
                            parsed_value.append(float(v))
                    except ValueError:
                        parsed_value.append(v)

                if len(value) == 1:
                    dic[key] = parsed_value[0]

                else:
                    dic[key] = parsed_value

            except IndexError:
                pass

        return dic

    def _parse_peak(self, peak):
        """
        Parses a single peak into a dictionary, the input being a list
        that corresponds to a single peak in a sparky save file. In addition,
        it parses all values in the list to do the following: (i) convert the
        values to float wherever possible and (ii) if the list has a single
        item, upack and return that item alone as the value. Currently assumes
        the following structure for a single peak:

        type peak
        ...
        [
        type label
        ...
        ]

        """
        label_pos = [i for i, word in enumerate(peak) if word in "[]"]

        # see if label exists and split things up
        # according to label position
        if label_pos:
            peak_info = peak[:label_pos[0]]
            label_info = peak[label_pos[0] + 1 : label_pos[1]]
        else:
            peak_info = peak
            label_info = []

        # peak_data dict stores all relevant information for a peak,
        # including all information about the label, label location, etc
        peak_data = {}

        # parse and add data under "peak"
        for item in peak_info:
            it = item.split()

            key, *vals = it
            for i, v in enumerate(vals):
                try:
                    vals[i] = float(v)
                except ValueError:
                    vals[i] = str(v)

            if len(vals) == 1:
                vals = vals[0]

            peak_data[key] = vals

        # parse and add data under "label"
        for item in label_info:
            key, *vals = item.split()
            if key not in peak_data.keys():
                peak_data[key] = vals

        try:
            for i, position in enumerate(peak_data["xy"]):
                position = position.split(",")
                peak_data["xy"][i] = [float(pos) for pos in position]
        except KeyError:
            pass

        return peak_data

    def _parse_ornaments(self, data):
        """
        Parses a string containing all ornaments into a dictionary. This
        is for all the data inside the <ornament> tag. The key for each
        peak item is given by the peak ID, which should be unique for each
        peak. The following structure is assumed:

        type peak # peak 1
        ...
        type peak # peak 2
        ...
        type peak # peak 3
        ...

        """

        data = data.split("\n")
        p = [i for i, word in enumerate(data) if word == "type peak"]
        peaklist = [data[p[i]: p[i+1]] for i in range(len(p)-1)]

        dic = {}
        for peak in peaklist:
            d = self._parse_peak(peak)
            dic[int(d["id"])] = d

        return dic


    def handle_starttag(self, tag, attrs):

        if tag in ["sparky", "version"]:
            self.curdict = self.spectrum
            self.spectrum[tag] = attrs[0][0]

        elif tag == "user":
            self.curdict = self.user

        elif tag == "spectrum":
            self.curdict = self.spectrum

        elif tag == "view":
            self.viewnum += 1
            self.view[self.viewnum] = {}
            self.curdict = self.view[self.viewnum]

        elif tag == "ornament":
            self.curdict = self.ornament

        else:
            # params tag for each view
            self.curtag = tag


    def handle_endtag(self, tag):
        self.curtag = None


    def handle_data(self, data):

        # ignore blank lines
        if len(data.strip()) == 0:
            pass

        elif self.curtag not in self.curdict.keys():

            # all the files that are read in a split at the newline character
            if (self.curtag is None) and (self.curdict == self.spectrum):
                dic = self._parse_info(data.split("\n"),)
                for k, v in dic.items():
                    self.curdict[k] = v

            elif (self.curtag is None) and (self.curdict == self.user):
                dic = self._parse_info(data.split("\n"), "user")
                for k, v in dic.items():
                    self.curdict[k] = v

            elif (self.curtag is None) and (self.curdict == self.view[self.viewnum]):
                dic = self._parse_info(data.split("\n"),)
                for k, v in dic.items():
                    self.curdict[k] = v

            elif (self.curtag is None) and (self.curdict == self.ornament):
                self.ornament = self._parse_ornaments(data)

            else:
                # for a param tag inside a view
                dic = self._parse_info(data.split("\n"),)
                self.curdict[self.curtag] = [dic]

        else:
            # this is only executed for multiple params tags in a view
            dic = self._parse_info(data.split("\n"),)
            self.curdict[self.curtag].append(dic)


def read_savefile(filename, spectrum_file=None):
    """
    Reads in a Sparky .save file and the corresponding spectrum (.ucsf)
    file. In addition to the usual dictionary contents that come with
    a .ucsf file, these additional dictionary keys are created with the content
    from .save file: "spectrum", "view", "user" and "ornament". The together
    contain all edits and annotations. By default, it tries to read in
    the spectrum file given in the .save file (but this fails many times due
    to relative paths in .save file)

    Parameters
    ----------
    savefile : str
        Filename of Sparky .save file.
    spectrum_file : str
        Filename of Sparky .ucsf file.


    Returns
    -------
    dic : dict
        Dictionary of Sparky .ucsf and .save parameters.
    data : ndarray
        Array of NMR data.

    """

    with open(filename, "r") as f:
        savefile = f.read().replace("<end ", r"</")

    parser = SparkySaveParser()
    parser.feed(savefile)
    parser.close()

    dic = {
        "user": parser.user,
        "spectrum": parser.spectrum,
        "view": parser.view,
        "ornament": parser.ornament,
    }


    if spectrum_file is None:
        try:
            spectrum_file = dic["spectrum"]["abspathname"]
        except KeyError:
            spectrum_file = os.path.join(os.path.dirname(filename), dic["spectrum"]["pathname"])

    try:
        d, data = read(spectrum_file)
    except FileNotFoundError:
        warn(f"Cannot find this file: {spectrum_file}")
        d, data = {}, None
    except UnicodeDecodeError:
        warn(f"Could not correctly parse the spectrum file {spectrum_file}")
        d, data = {}, None

    dic.update(d)

    return dic, data


# sparky_ low memory objects
class sparky_2d(fileiobase.data_nd):
    """
    Emulates a ndarray object without loading data into memory for low memory
    reading of 2D Sparky files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filename : str
        Filename of 2D Sparky file.
    order : tuple, optional
        Order of axes against file.  None is equivalent to (0, 1).

    """

    def __init__(self, filename, order=None):
        """
        Create and set up object
        """
        # open the file
        self.filename = filename
        f = open(filename, 'rb')

        # read the fileheader
        self.dic = fileheader2dic(get_fileheader(f))

        if self.dic["naxis"] != 2:
            raise Exception("file is not a 2D Sparky file")

        # read in the axisheaders
        self.dic["w1"] = axisheader2dic(get_axisheader(f))
        self.dic["w2"] = axisheader2dic(get_axisheader(f))
        f.close()

        # sizes
        self.lenY = self.dic["w1"]["npoints"]
        self.lenX = self.dic["w2"]["npoints"]

        # tile sizes
        self.lentY = self.dic["w1"]["bsize"]
        self.lentX = self.dic["w2"]["bsize"]

        # check order
        if order is None:
            order = (0, 1)

        # finalize
        self.dtype = np.dtype("float32")
        self.order = order
        self.fshape = (self.lenY, self.lenX)
        self.__setdimandshape__()

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = sparky_2d(self.filename, order)
        return n

    def __fgetitem__(self, slices):
        """
        Returns ndarray of selected values.

        (sY, sX) is a well formatted tuple of slices
        """
        sY, sX = slices
        f = open(self.filename, 'rb')

        # print(sY,sX)
        gY = range(self.lenY)[sY]  # list of values to take in Y
        gX = range(self.lenX)[sX]  # list of values to take in X

        # tiles to get in each dim to read
        # Y tile to read
        gtY = set(int(np.floor(i / self.lentY)) for i in gY)
        # X tile to read
        gtX = set(int(np.floor(i / self.lentX)) for i in gX)

        # create a empty output directory
        out = np.empty((len(gY), len(gX)), dtype=self.dtype)

        for iY in gtY:      # loop over Y tiles to get
            for iX in gtX:  # loop over X tiles to get

                # get the tile and reshape it
                ntile = int(iY * np.ceil(self.lenX / self.lentX) + iX)
                tile = get_tilen(f, ntile, (self.lentX, self.lentY))
                tile = tile.reshape(self.lentY, self.lentX)

                # tile minimum and max values for each dim
                minX = iX * self.lentX
                maxX = (iX + 1) * self.lentX

                minY = iY * self.lentY
                maxY = (iY + 1) * self.lentY

                # determine what elements are needed from this tile
                XinX = [i for i in gX if maxX > i >= minX]  # values in gX
                XinT = [i - minX for i in XinX]  # tile index values
                XinO = [gX.index(i) for i in XinX]  # output indexes

                YinY = [i for i in gY if maxY > i >= minY]  # values in gX
                YinT = [i - minY for i in YinY]  # tile index values
                YinO = [gY.index(i) for i in YinY]  # output indexes

                # take elements from the tile
                ctile = tile.take(XinT, axis=1).take(YinT, axis=0)

                # DEBUGGING info
                # print("-------------------------------")
                # print("iX:",iX,"iY:",iY,"ntile:",ntile)
                # print("tile.shape",tile.shape)
                # print("minX:",minX,"maxX",maxX)
                # print("minY:",minY,"maxY",maxY)
                # print("XinX",XinX)
                # print("XinT",XinT)
                # print("XinO",XinO)
                # print("YinY",YinY)
                # print("YinT",YinT)
                # print("YinO",YinO)

                # put the cut tile to the out array (uses some fancy indexing)
                out[np.ix_(YinO, XinO)] = ctile

        f.close()
        return out


class sparky_3d(fileiobase.data_nd):
    """
    Emulates a ndarray object without loading data into memory for low memory
    read of 3D Sparky files.

    * slicing operations return ndarray objects.
    * can iterate over with expected results.
    * transpose and swapaxes methods create a new objects with correct axes
      ordering.
    * has ndim, shape, and dtype attributes.

    Parameters
    ----------
    filename : str
        Filename of 3D Sparky file.
    order : tuple
        Ordering of axes against file. None is equilent to (0, 1, 2)

    """
    def __init__(self, filename, order=None):
        """
        Create and set up object
        """
        # open the file
        self.filename = filename
        f = open(filename, 'rb')

        # read the fileheader
        self.dic = fileheader2dic(get_fileheader(f))

        if self.dic["naxis"] != 3:
            raise Exception("file not 3D Sparky file")

        # read in the axisheaders
        self.dic["w1"] = axisheader2dic(get_axisheader(f))
        self.dic["w2"] = axisheader2dic(get_axisheader(f))
        self.dic["w3"] = axisheader2dic(get_axisheader(f))
        f.close()

        # sizes
        self.lenZ = self.dic["w1"]["npoints"]
        self.lenY = self.dic["w2"]["npoints"]
        self.lenX = self.dic["w3"]["npoints"]

        # tile sizes
        self.lentZ = self.dic["w1"]["bsize"]
        self.lentY = self.dic["w2"]["bsize"]
        self.lentX = self.dic["w3"]["bsize"]

        # check order
        if order is None:
            order = (0, 1, 2)

        # finalize
        self.dtype = np.dtype("float32")
        self.order = order
        self.fshape = (self.lenZ, self.lenY, self.lenX)
        self.__setdimandshape__()

    def __fcopy__(self, order):
        """
        Create a copy
        """
        n = sparky_3d(self.filename, order)
        return n

    def __fgetitem__(self, slices):
        """
        Returns ndarray of selected values.

        (sZ, sY, sX) is a well formateed tuple of slices
        """
        sZ, sY, sX = slices
        f = open(self.filename, 'rb')

        gZ = range(self.lenZ)[sZ]  # list of values to take in Z
        gY = range(self.lenY)[sY]  # list of values to take in Y
        gX = range(self.lenX)[sX]  # list of values to take in X

        # tiles to get in each dim to read
        # Z tiles
        gtZ = set(int(np.floor(float(i) / self.lentZ)) for i in gZ)
        # Y tiles
        gtY = set(int(np.floor(float(i) / self.lentY)) for i in gY)
        # X tiles
        gtX = set(int(np.floor(float(i) / self.lentX)) for i in gX)

        # total tiles in each dim
        ttX = int(np.ceil(self.lenX / float(self.lentX)))  # total tiles in X
        ttY = int(np.ceil(self.lenY / float(self.lentY)))  # total tiles in Y
        ttZ = int(np.ceil(self.lenZ / float(self.lentZ)))  # total tiles in Z

        tile_tup = (self.lentZ, self.lentY, self.lentX)

        # create a empty output array
        out = np.empty((len(gZ), len(gY), len(gX)), dtype=self.dtype)

        for iZ in gtZ:          # loop over Z tiles to get
            for iY in gtY:      # loop over Y tiles to get
                for iX in gtX:  # loop over X tiles to get

                    # get the tile and reshape it
                    ntile = iZ * ttX * ttY + iY * ttX + iX
                    tile = get_tilen(f, ntile, tile_tup)
                    tile = tile.reshape(tile_tup)

                    # tile minimum and max values for each dim
                    minX = iX * self.lentX
                    maxX = (iX + 1) * self.lentX

                    minY = iY * self.lentY
                    maxY = (iY + 1) * self.lentY

                    minZ = iZ * self.lentZ
                    maxZ = (iZ + 1) * self.lentZ

                    # determine what elements are needed from this tile
                    XinX = [i for i in gX if maxX > i >= minX]  # values in gX
                    XinT = [i - minX for i in XinX]  # tile index values
                    XinO = [gX.index(i) for i in XinX]  # output indexes

                    YinY = [i for i in gY if maxY > i >= minY]  # values in gX
                    YinT = [i - minY for i in YinY]  # tile index values
                    YinO = [gY.index(i) for i in YinY]  # output indexes

                    ZinZ = [i for i in gZ if maxZ > i >= minZ]  # values in gX
                    ZinT = [i - minZ for i in ZinZ]  # tile index values
                    ZinO = [gZ.index(i) for i in ZinZ]  # output indexes

                    # take elements from the tile
                    ctile = tile.take(XinT, axis=2).take(YinT, axis=1)
                    ctile = ctile.take(ZinT, axis=0)

                    # DEBUGGING info
                    # print("-------------------------------")
                    # print("iX:",iX,"iY:",iY,"iZ:",iZ,"ntile:",ntile)
                    # print("ttX:",ttX,"ttY:",ttY,"ttZ",ttZ)
                    # print("tile.shape",tile.shape)
                    # print("minX:",minX,"maxX",maxX)
                    # print("minY:",minY,"maxY",maxY)
                    # print("minZ:",minZ,"maxZ",maxZ)
                    # print("XinX",XinX)
                    # print("XinT",XinT)
                    # print("XinO",XinO)
                    # print("YinY",YinY)
                    # print("YinT",YinT)
                    # print("YinO",YinO)
                    # print("ZinZ",ZinZ)
                    # print("ZinT",ZinT)
                    # print("ZinO",ZinO)

                    # put the cut tile to the out array
                    out[np.ix_(ZinO, YinO, XinO)] = ctile
        f.close()
        return out


# tile and data get/put functions
def get_tilen(f, n_tile, tw_tuple):
    """
    Read a tile from a Sparky file object.

    Parameters
    ----------
    f : file object
        Open file object pointing to a Sparky file.
    n_tile : int
        Tile number to read
    tw_tuple : tuple of ints
        Tile size

    Returns
    -------
    tile : ndarray
        Tile of NMR data. Data is returned as a 1D array.

    Notes
    -----
    Current file position is loss. In can be stored before calling if the
    position is later needed.

    """
    # determine the size of the tile in bytes
    tsize = 4
    for i in tw_tuple:
        tsize = tsize * i

    # seek to the beginning of the tile
    f.seek(int(180 + 128 * len(tw_tuple) + n_tile * tsize))
    return np.frombuffer(f.read(tsize), dtype='>f4')

def untile_data4D(data, tile_size, data_size):
    """
    Rearrange 4D tiled/Sparky formatted data into standard format.

    Parameters
    ----------
    data : 1D ndarray
        Tiled/Sparky formatted 2D NMR data.
    (lentA, lentZ, lentY, lentX) : tuple of ints
        Size of tile
    (lenA, lenZ, lenY, lenX) : tuple of ints
        Size of NMR data.

    Returns
    -------
    sdata : 4D ndarray
        NMR data, untiled/standard format.

    """
    lentA, lentZ, lentY, lentX = tile_size
    lenA, lenZ, lenY, lenX = data_size

    # determine the number of tiles in data
    ttX = int(np.ceil(lenX / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(lenY / float(lentY)))  # total tiles in Y dim
    ttZ = int(np.ceil(lenZ / float(lentZ)))  # total tiles in Z dim
    ttA = int(np.ceil(lenA / float(lentA)))  # total tiles in A dim
    tt = ttX * ttY * ttZ * ttA

    # calc some basic parameter
    tsize = lentX * lentY * lentZ *lentA # number of points in one tile
    t_tup = (lentA, lentZ, lentY, lentX)  # tile size tuple

    # create an empty array to store file data
    out = np.empty((ttA * lentA, ttZ * lentZ, ttY * lentY, ttX * lentX), dtype="float32")
    for iA in range(int(ttA)):
        for iZ in range(int(ttZ)):
            for iY in range(int(ttY)):
                for iX in range(int(ttX)):

                    minX = iX * lentX
                    maxX = (iX + 1) * lentX

                    minY = iY * lentY
                    maxY = (iY + 1) * lentY

                    minZ = iZ * lentZ
                    maxZ = (iZ + 1) * lentZ

                    minA = iA * lentA
                    maxA = (iA + 1) * lentA

                    ntile = iA * ttZ * ttY * ttX + iZ * ttX * ttY + iY * ttX + iX
                    minT = ntile * tsize
                    maxT = (ntile + 1) * tsize

                    out[minA:maxA, minZ:maxZ, minY:maxY, minX:maxX] =  \
                        data[minT:maxT].reshape(t_tup)

    return out[:lenA, :lenZ, :lenY, :lenX]

def get_tile(f, num_points):
    """
    Read the next tile from a Sparky file object.

    Parameters
    ----------
    f : file object
        Open file object pointing to a Sparky file.
    num_points : int
        Number of points in the tile.

    Returns
    -------
    tile : ndarray
        Tile of NMR data. Data is returned as a 1D array.

    """
    bsize = num_points * 4        # size in bytes
    return np.frombuffer(f.read(bsize), dtype='>f4')


def put_tile(f, tile):
    """
    Put a tile to a Sparky file object.

    Parameters
    ----------
    f : file object
        Open file object pointing to a Sparky file, to be written to.
    tile : ndarray
        Tile of NMR data to be written.

    """
    f.write(tile.astype('>f4').tobytes())
    return


def get_data(f):
    """
    Read all data from sparky file object.
    """
    return np.frombuffer(f.read(), dtype='>f4')


def put_data(f, data):
    """
    Put data to a Sparky file object.

    This function does not untile data. This should be done before calling
    this function

    """
    f.write(data.astype('>f4').tobytes())
    return


# tiling/untiling functions
def find_tilen_2d(data, ntile, tile_size):
    """
    Return a tile from a 2D NMR data set.

    Parameters
    ----------
    data : 2D ndarray
        NMR data, untiled/standard format.
    ntile : int
        Tile number to extract.
    (lentY, lentX) : tuple of ints
        Tile size (w1, w2).

    Returns
    -------
    tile : 1D ndarray
        Tile of NMR data, returned as 1D array.

    Notes
    -----
    Edge tiles are zero filled to the indicated tile size.

    """
    lentY, lentX = tile_size
    ttX = int(np.ceil(data.shape[1] / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(data.shape[0] / float(lentY)))  # total tiles in Y dim

    # tile number in each dim
    Xt = ntile % ttX
    Yt = int(np.floor(ntile / ttX))

    # dimension limits
    Xmin = int(Xt * lentX)
    Xmax = int((Xt + 1) * lentX)

    Ymin = int(Yt * lentY)
    Ymax = int((Yt + 1) * lentY)

    tile = data[Ymin:Ymax, Xmin:Xmax]

    # some edge tiles might need zero filling
    # see if this is the case
    if tile.shape == (lentY, lentX):    # well sized tile
        return tile.flatten()
    else:
        new_tile = np.zeros((lentY, lentX), dtype="float32")
        new_tile[:tile.shape[0], :tile.shape[1]] = tile
        return new_tile.flatten()


def tile_data2d(data, tile_size):
    """
    Tile 2D data into a 1D array.

    Parameters
    ----------
    data : 2D ndarray
        NMR data, untiled/standard format.
    (lentY, lentX) : tuple of ints
        Tile size.

    Returns
    -------
    tdata : 1D ndarray
        Tiled/Sparky formatted NMR data, returned as 1D array.

    """
    lentY, lentX = tile_size
    # determine the number of tiles in data
    ttX = int(np.ceil(data.shape[1] / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(data.shape[0] / float(lentY)))  # total tiles in Y dim
    tt = ttX * ttY  # total number of tiles

    # calc some basic parameter
    tsize = lentX * lentY   # number of points in one tile
    t_tup = (lentY, lentX)  # tile size tuple

    # create an empty array to store file data
    out = np.empty((tt * tsize), dtype="float32")

    for i in range(int(tt)):
        out[i * tsize:(i + 1) * tsize] = find_tilen_2d(data, i, t_tup)

    return out


def untile_data2D(data, tile_size, data_size):
    """
    Rearrange 2D Tiled/Sparky formatted data into standard format.

    Parameters
    ----------
    data : 1D ndarray
        Tiled/Sparky formatted 2D NMR data.
    (lentY, lenX) : tuple of ints
        Size of tile.
    (lenY, lenX) : tuple of ints
        Size of NMR data.

    Returns
    -------
    sdata : 2D ndarray
        NMR data, untiled/standard format.

    """
    lentY, lentX = tile_size
    lenY, lenX = data_size
    # determine the number of tiles in data
    ttX = int(np.ceil(lenX / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(lenY / float(lentY)))  # total tiles in Y dim
    tt = ttX * ttY

    # calc some basic parameter
    tsize = lentX * lentY  # number of points in one tile
    t_tup = (lentY, lentX)  # tile size tuple

    # create an empty array to store file data
    out = np.empty((ttY * lentY, ttX * lentX), dtype="float32")

    for iY in range(int(ttY)):
        for iX in range(int(ttX)):
            minX = iX * lentX
            maxX = (iX + 1) * lentX

            minY = iY * lentY
            maxY = (iY + 1) * lentY

            ntile = iY * ttX + iX
            minT = ntile * tsize
            maxT = (ntile + 1) * tsize

            # DEBUG
            # print("ntile",ntile)
            # print("minX",minX,"maxX",maxX)
            # print("minY",minY,"maxY",maxY)
            # print("minT",minT,"maxT",maxT)

            # print(out[minY:maxY,minX:maxX].shape)
            # print(data[minT:maxT].reshape(t_tup).shape)

            out[minY:maxY, minX:maxX] = data[minT:maxT].reshape(t_tup)

    return out[:lenY, :lenX]


def find_tilen_3d(data, ntile, tile_size):
    """
    Return a single tile from a 3D NMR data set.

    Parameters
    ----------
    data : 3D ndarray
        NMR data, untiled/standard format.
    ntile : int
        Tile number to extract.
    (lentZ, lentY, lentX) : tuple of ints
        Tile size (w1, w2, w3).

    Returns
    -------
    tile : 1D ndarray
        Tile of NMR data, returned as 1D array.

    Notes
    -----
    Edge tiles are zero filled to the indicated tile size.

    """
    lentZ, lentY, lentX = tile_size
    ttX = int(np.ceil(data.shape[2] / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(data.shape[1] / float(lentY)))  # total tiles in Y dim
    ttZ = int(np.ceil(data.shape[0] / float(lentZ)))  # total tiles in Z dim

    # tile number in each dim
    Xt = ntile % ttX
    Yt = int(np.floor(ntile / ttX)) % ttY
    Zt = int(np.floor(ntile / (ttX * ttY)))

    # dimension limits
    Xmin = int(Xt * lentX)
    Xmax = int((Xt + 1) * lentX)

    Ymin = int(Yt * lentY)
    Ymax = int((Yt + 1) * lentY)

    Zmin = int(Zt * lentZ)
    Zmax = int((Zt + 1) * lentZ)

    tile = data[Zmin:Zmax, Ymin:Ymax, Xmin:Xmax]

    # some edge tiles might need zero filling
    # see if this is the case
    if tile.shape == (lentZ, lentY, lentX):  # well sized tile
        return tile.flatten()
    else:
        new_tile = np.zeros((lentZ, lentY, lentX), dtype="float32")
        new_tile[:tile.shape[0], :tile.shape[1], :tile.shape[2]] = tile
        return new_tile.flatten()


def tile_data3d(data, tile_size):
    """
    Tile 3D data into a 1D numpy array

    Parameters
    ----------
    data : 3D ndarray
        NMR data, untiled/standard format.
    (lentZ, lentY, lentX) : tuple of ints
        Tile size (w1, w2, w3).

    Returns
    -------
    tile : 1D ndarray
        Tiled/Sparky formatted NMR data, returned as 1D array.

    """
    lentZ, lentY, lentX = tile_size
    # determine the number of tiles in data
    ttX = int(np.ceil(data.shape[2] / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(data.shape[1] / float(lentY)))  # total tiles in Y dim
    ttZ = int(np.ceil(data.shape[0] / float(lentZ)))  # total tiles in Z dim

    tt = ttX * ttY * ttZ    # total number of tiles

    # calc some basic parameter
    tsize = lentX * lentY * lentZ  # number of points in one tile
    t_tup = (lentZ, lentY, lentX)  # tile size tuple

    # create an empty array to store file data
    out = np.empty((tt * tsize), dtype="float32")

    for i in range(int(tt)):
        out[i * tsize:(i + 1) * tsize] = find_tilen_3d(data, i, t_tup)
    return out


def untile_data3D(data, tile_size, data_size):
    """
    Rearrange 3D tiled/Sparky formatted data into standard format.

    Parameters
    ----------
    data : 1D ndarray
        Tiled/Sparky formatted 2D NMR data.
    (lentZ, lentY, lentX) : tuple of ints
        Size of tile
    (lenZ, lenY, lenX) : tuple of ints
        Size of NMR data.

    Returns
    -------
    sdata : 3D ndarray
        NMR data, untiled/standard format.

    """
    lentZ, lentY, lentX = tile_size
    lenZ, lenY, lenX = data_size

    # determine the number of tiles in data
    ttX = int(np.ceil(lenX / float(lentX)))  # total tiles in X dim
    ttY = int(np.ceil(lenY / float(lentY)))  # total tiles in Y dim
    ttZ = int(np.ceil(lenZ / float(lentZ)))  # total tiles in Z dim
    tt = ttX * ttY * ttZ

    # calc some basic parameter
    tsize = lentX * lentY * lentZ  # number of points in one tile
    t_tup = (lentZ, lentY, lentX)  # tile size tuple

    # create an empty array to store file data
    out = np.empty((ttZ * lentZ, ttY * lentY, ttX * lentX), dtype="float32")

    for iZ in range(int(ttZ)):
        for iY in range(int(ttY)):
            for iX in range(int(ttX)):

                minX = iX * lentX
                maxX = (iX + 1) * lentX

                minY = iY * lentY
                maxY = (iY + 1) * lentY

                minZ = iZ * lentZ
                maxZ = (iZ + 1) * lentZ

                ntile = iZ * ttX * ttY + iY * ttX + iX
                minT = ntile * tsize
                maxT = (ntile + 1) * tsize

                out[minZ:maxZ, minY:maxY, minX:maxX] =  \
                    data[minT:maxT].reshape(t_tup)

    return out[:lenZ, :lenY, :lenX]


# fileheader functions
def get_fileheader(f):
    """
    Get fileheader from file and return a list.

    Reads the 180 byte file header of a Sparky file

    """
    # file header as described in ucsffile.cc of sparky source
    # header is packed as follows:
    # ident(10s),naxis(c),ncomponents(c),encoding(c),version(c)
    # owner(9s),date(26s),comment(80s),pad(3x),seek_pos(l),scratch(40s),
    # pad(4x)

    # note that between comment and seek_pos is a 3 byte pad
    # so that the long is @ a multiple of 4
    # also sparky always packs big-endian, hence >
    return struct.unpack('>10s 4c 9s 26s 80s 3x l 40s 4x', f.read(180))


def put_fileheader(f, fl):
    """
    Write fileheader list to file (180-bytes).
    """
    f.write(struct.pack('>10s 4c 9s 26s 80s 3x l 40s 4x', *fl))
    return


def fileheader2dic(header):
    """
    Convert a fileheader list into a Sparky parameter dictionary.
    """
    dic = dict()
    dic["ident"] = str(header[0].decode()).strip('\x00')
    dic["naxis"] = ord(header[1].decode())
    dic["ncomponents"] = ord(header[2].decode())
    dic["encoding"] = ord(header[3].decode())
    dic["version"] = ord(header[4].decode())
    dic["owner"] = str(header[5].decode()).strip('\x00')
    dic["date"] = str(header[6].decode()).strip('\x00')
    dic["comment"] = str(header[7].decode()).strip('\x00')
    dic["seek_pos"] = header[8]     # eof seek position
    dic["scratch"] = str(header[9].decode()).strip('\x00')
    return dic


def dic2fileheader(dic):
    """
    Convert a Sparky parameter dictionary into a fileheader list.
    """
    fl = [0] * 10
    fl[0] = dic["ident"].encode()
    fl[1] = chr(dic["naxis"]).encode()
    fl[2] = chr(dic["ncomponents"]).encode()
    fl[3] = chr(dic["encoding"]).encode()
    fl[4] = chr(dic["version"]).encode()
    fl[5] = dic["owner"].encode()
    fl[6] = dic["date"].encode()
    fl[7] = dic["comment"].encode()
    fl[8] = dic["seek_pos"]
    fl[9] = dic["scratch"].encode()
    return fl


# axisheader functions
def get_axisheader(f):
    """
    Get an axisheader from file and return a list.

    Only the first 44 bytes are examined, the NMR_PROCESSED and other header
    parameters are ignored since the current version of Sparky does not use
    them.

    """
    # axis header is described in ucsffile.cc
    # axis header is packed as follows
    # nucleus(6s),spectral_shift(h),npoints(I),size(I),bsize(I)
    # spectrometer_freq(f),spectral_width(f),xmtr_freq(f),zero_order(f),
    # first_order(f),first_pt_scale(f),ZEROS
    return struct.unpack('>6s h 3I 6f 84s', f.read(128))


def put_axisheader(f, al):
    """
    Write an axisheader list to file (128-bytes written).
    """
    f.write(struct.pack('>6s h 3I 6f 84s', *al))
    return


def axisheader2dic(header):
    """
    Convert an axisheader list into Sparky parameter axis dictionary.
    """
    dic = dict()
    dic["nucleus"] = str(header[0].decode()).strip('\x00')
    dic["spectral_shift"] = header[1]
    dic["npoints"] = header[2]
    dic["size"] = header[3]
    dic["bsize"] = header[4]
    dic["spectrometer_freq"] = header[5]
    dic["spectral_width"] = header[6]
    dic["xmtr_freq"] = header[7]
    dic["zero_order"] = header[8]
    dic["first_order"] = header[9]
    dic["first_pt_scale"] = header[10]
    dic["extended"] = header[11]
    return dic


def dic2axisheader(dic):
    """
    Convert a Sparky parameter axis diction into a axisherder list.
    """
    al = [0] * 12
    al[0] = dic["nucleus"].encode()
    al[1] = dic["spectral_shift"]
    al[2] = dic["npoints"]
    al[3] = dic["size"]
    al[4] = dic["bsize"]
    al[5] = dic["spectrometer_freq"]
    al[6] = dic["spectral_width"]
    al[7] = dic["xmtr_freq"]
    al[8] = dic["zero_order"]
    al[9] = dic["first_order"]
    al[10] = dic["first_pt_scale"]
    al[11] = dic["extended"]
    return al

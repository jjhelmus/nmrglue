import nmrglue as ng
import numpy as np
import pytest


_ONE_D_PEAKS = {
    "shape": (1024,),
    "positions": [100, 200, 300, 400, 500],
    "lws": [10, 20, 10, 20, 10],
    "amps": [100, 200, 300, 300, 150],
    "vparams": [0.2, 0.4, 0.2, 0.5, 1.0],
}

_TWO_D_PEAKS = {
    "shape": (512, 128),
    "positions": [(100, 100), (200, 53), (300, 110), (400, 75)],
    "lws": [(10, 5), (20, 5), (10, 5), (20, 5)],
    "amps": [100, 200, 300, 300],
    "vparams": [(0.2, 0.2), (0.4, 0.4), (0.2, 0.2), (0.5, 0.4)],
}


_THREE_D_PEAKS = {
    "shape": (256, 64, 64),
    "positions": [(150, 24, 22), (200, 10, 50), (210, 50, 10), (77, 30, 15)],
    "lws": [(5, 5, 5), (3, 8, 5), (7, 5, 5), (7, 6, 5)],
    "amps": [100, 200, 300, 300],
    "vparams": [(0.2, 0.2, 0.3), (0.5, 0.4, 0.4), (0.1, 0.2, 0.2), (0.3, 0.5, 0.4)],
}


def _generate_1d_data(shape, positions, lws, amps, vparams):
    """
    Generates a test 1d dataset with multiple peaks

    Parameters
    ----------
    shape : Iterable[int]
        shape of the numpy array to be created
    positions : Iterable[int]
        a list of the positions of the peaks
    lws : Iterable[float|int]
        a list of linewidths for each peak
    amps : Iterable[float|int]
        a list of amplitudes for each peak
    vparams : Iterable[float|int]
        a list of list containing the eta parameter
        for the pseudo voigt lineshape

    Returns
    -------
    numpy.ndarray
        simulated one-d dataset

    """
    data = ng.linesh.sim_NDregion(
        shape=shape,
        lineshapes=["pv"],
        params=[[(pos, lws, vp)] for pos, lws, vp in zip(positions, lws, vparams)],
        amps=amps,
    )
    return data


def _generate_2d_data(dataset):
    """
    Generates a test 2d dataset with multiple peaks

    Parameters
    ----------
    shape : Iterable[Iterable[int, int]]
        shape of the numpy array to be created
    positions : Iterable[Iterable[int, int]]
        a list of list of two positions for each peak
    lws : Iterable[Iterable[float, float]]
        a list of list of two linewidths for each peak
    amps : Iterable[Iterable[float, float]]
        a list of list of two amplitudes for each peak
    vparams : Iterable[Iterable[float, float]]
        a list of list containing 2 values for the
        eta parameter for the pseud-voigt lineshape

    Returns
    -------
    numpy.ndarray
        simulated two-d dataset

    """
    params = []
    for i in range(len(dataset["positions"])):
        d = [[], []]
        for j in range(2):
            d[j] = (
                dataset["positions"][i][j],
                dataset["lws"][i][j],
                dataset["vparams"][i][j],
            )

        params.append(d)

    data = ng.linesh.sim_NDregion(
        shape=(512, 128), lineshapes=["pv", "pv"], params=params, amps=dataset["amps"]
    )

    return data


def _generate_3d_data(dataset):
    """
    Generates a test 3d dataset with multiple peaks

    Parameters
    ----------
    shape : Iterable[Iterable[int, int, int]]
        shape of the numpy array to be created
    positions : Iterable[Iterable[int, int, int]]
        a list of list of three positions for each peak
    lws : Iterable[Iterable[float, float, float]]
        a list of list of three linewidths for each peak
    amps : Iterable[Iterable[float, float, float]]
        a list of list of three amplitudes for each peak
    vparams : Iterable[Iterable[float, float, float]]
        a list of list containing three values for the
        eta parameter for the pseud-voigt lineshape

    Returns
    -------
    numpy.ndarray
        simulated three-d dataset

    """
    params = []
    for i in range(len(dataset["positions"])):
        d = [[], [], []]
        for j in range(3):
            d[j] = (
                dataset["positions"][i][j],
                dataset["lws"][i][j],
                dataset["vparams"][i][j],
            )

        params.append(d)

    data = ng.linesh.sim_NDregion(
        shape=(256, 64, 64),
        lineshapes=["pv", "pv", "pv"],
        params=params,
        amps=dataset["amps"],
    )

    return data


def _test_1d(dataset, algorithm, msep=None, rtol=1):
    """test 1d peak picking"""

    data = _generate_1d_data(**dataset)
    peaks = ng.peakpick.pick(data, pthres=50, algorithm=algorithm, msep=msep)
    assert np.allclose(peaks.X_AXIS, dataset["positions"], rtol=1)


def _test_2d(dataset, algorithm, msep=None, rtol=1):
    """test 2d peak picking"""

    data = _generate_2d_data(dataset)
    peaks = ng.peakpick.pick(data, pthres=50, algorithm=algorithm, msep=msep)
    assert np.allclose(peaks.X_AXIS, [i[1] for i in dataset["positions"]], rtol=1)
    assert np.allclose(peaks.Y_AXIS, [i[0] for i in dataset["positions"]], rtol=1)


def _test_3d(dataset, algorithm, msep=None, rtol=1):
    """test 3d peak picking"""

    data = _generate_3d_data(dataset)
    peaks = ng.peakpick.pick(data, pthres=50, algorithm=algorithm, msep=msep)

    assert np.allclose(
        sorted(peaks.X_AXIS), sorted([i[2] for i in dataset["positions"]]), rtol=rtol
    )
    assert np.allclose(
        sorted(peaks.Y_AXIS), sorted([i[1] for i in dataset["positions"]]), rtol=rtol
    )
    assert np.allclose(
        sorted(peaks.Z_AXIS), sorted([i[0] for i in dataset["positions"]]), rtol=rtol
    )


@pytest.mark.fast
def test_1d_connected():
    _test_1d(dataset=_ONE_D_PEAKS, algorithm="connected")


@pytest.mark.fast
def test_1d_downward():
    _test_1d(dataset=_ONE_D_PEAKS, algorithm="downward")


@pytest.mark.fast
def test_1d_thres():
    _test_1d(dataset=_ONE_D_PEAKS, algorithm="thres", msep=(5,))


@pytest.mark.fast
def test_1d_thres_fast():
    _test_1d(dataset=_ONE_D_PEAKS, algorithm="thres-fast", msep=(5,))


@pytest.mark.fast
def test_2d_connected():
    _test_2d(dataset=_TWO_D_PEAKS, algorithm="connected")


@pytest.mark.fast
def test_2d_downward():
    _test_2d(dataset=_TWO_D_PEAKS, algorithm="downward")


@pytest.mark.fast
def test_2d_thres():
    _test_2d(dataset=_TWO_D_PEAKS, algorithm="thres", msep=(5, 5))


@pytest.mark.fast
def test_2d_thres_fast():
    _test_2d(dataset=_TWO_D_PEAKS, algorithm="thres-fast", msep=(5, 5))


@pytest.mark.fast
def test_3d_connected():
    _test_3d(dataset=_THREE_D_PEAKS, algorithm="connected")


@pytest.mark.fast
def test_3d_downward():
    _test_3d(dataset=_THREE_D_PEAKS, algorithm="downward")


@pytest.mark.fast
def test_3d_thres():
    _test_3d(dataset=_THREE_D_PEAKS, algorithm="thres", msep=(2, 2, 2))


@pytest.mark.fast
def test_3d_thres_fast():
    _test_3d(dataset=_THREE_D_PEAKS, algorithm="thres-fast", msep=(2, 2, 2))

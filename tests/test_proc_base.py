import numpy as np
import nmrglue as ng
import scipy.linalg
import scipy.signal


def test_reorder_nus_data_2d():
    """reorder 2d nus data"""
    nus_data = np.linspace(0, 1, 128 * 4 * 2).reshape(-1, 128)
    nuslist = [(0,), (2,), (5,), (7,)]
    full_shape = (16, 128)
    full_data = ng.proc_base.expand_nus(
        data=nus_data, shape=full_shape, nuslist=nuslist
    )

    assert full_data.shape == full_shape
    assert np.allclose(full_data[0], nus_data[0], atol=1e-7)
    assert np.allclose(full_data[1], nus_data[1], atol=1e-7)
    assert np.allclose(full_data[4], nus_data[2], atol=1e-7)
    assert np.allclose(full_data[5], nus_data[3], atol=1e-7)
    assert np.allclose(full_data[10], nus_data[4], atol=1e-7)
    assert np.allclose(full_data[11], nus_data[5], atol=1e-7)
    assert np.allclose(full_data[14], nus_data[6], atol=1e-7)
    assert np.allclose(full_data[15], nus_data[7], atol=1e-7)


def test_reorder_nus_data_3d():
    """reorder 3d nus data"""
    nus_data = np.linspace(0, 1, 128 * 4 * 2 * 2).reshape(-1, 128)
    nuslist = [(0, 0), (2, 3), (5, 3), (6, 7)]
    full_shape = (16, 16, 128)
    full_data = ng.proc_base.expand_nus(
        data=nus_data, shape=full_shape, nuslist=nuslist
    )

    assert full_data.shape == full_shape
    assert np.allclose(full_data[0, 0], nus_data[0], atol=1e-7)
    assert np.allclose(full_data[0, 1], nus_data[1], atol=1e-7)
    assert np.allclose(full_data[1, 0], nus_data[2], atol=1e-7)
    assert np.allclose(full_data[1, 1], nus_data[3], atol=1e-7)

    assert np.allclose(full_data[6, 4], nus_data[4], atol=1e-7)
    assert np.allclose(full_data[6, 5], nus_data[5], atol=1e-7)
    assert np.allclose(full_data[7, 4], nus_data[6], atol=1e-7)
    assert np.allclose(full_data[7, 5], nus_data[7], atol=1e-7)

    assert np.allclose(full_data[6, 10], nus_data[8], atol=1e-7)
    assert np.allclose(full_data[6, 11], nus_data[9], atol=1e-7)
    assert np.allclose(full_data[7, 10], nus_data[10], atol=1e-7)
    assert np.allclose(full_data[7, 11], nus_data[11], atol=1e-7)

    assert np.allclose(full_data[14, 12], nus_data[12], atol=1e-7)
    assert np.allclose(full_data[14, 13], nus_data[13], atol=1e-7)
    assert np.allclose(full_data[15, 12], nus_data[14], atol=1e-7)
    assert np.allclose(full_data[15, 13], nus_data[15], atol=1e-7)


def test_reorder_nus_data_4d():
    """reorder 4d nus data"""
    nus_data = np.linspace(0, 1, 32 * 2 * 2 * 2 * 2).reshape(-1, 32)
    nuslist = [(0, 0, 0), (2, 3, 1)]
    full_shape = (8, 8, 8, 32)
    full_data = ng.proc_base.expand_nus(
        data=nus_data, shape=full_shape, nuslist=nuslist
    )

    assert full_data.shape == full_shape
    assert np.allclose(full_data[0, 0, 0], nus_data[0], atol=1e-7)
    assert np.allclose(full_data[0, 0, 1], nus_data[1], atol=1e-7)
    assert np.allclose(full_data[0, 1, 0], nus_data[2], atol=1e-7)
    assert np.allclose(full_data[0, 1, 1], nus_data[3], atol=1e-7)
    assert np.allclose(full_data[1, 0, 0], nus_data[4], atol=1e-7)
    assert np.allclose(full_data[1, 0, 1], nus_data[5], atol=1e-7)
    assert np.allclose(full_data[1, 1, 0], nus_data[6], atol=1e-7)
    assert np.allclose(full_data[1, 1, 1], nus_data[7], atol=1e-7)

    assert np.allclose(full_data[2, 6, 4], nus_data[8], atol=1e-7)
    assert np.allclose(full_data[2, 6, 5], nus_data[9], atol=1e-7)
    assert np.allclose(full_data[2, 7, 4], nus_data[10], atol=1e-7)
    assert np.allclose(full_data[2, 7, 5], nus_data[11], atol=1e-7)
    assert np.allclose(full_data[3, 6, 4], nus_data[12], atol=1e-7)
    assert np.allclose(full_data[3, 6, 5], nus_data[13], atol=1e-7)
    assert np.allclose(full_data[3, 7, 4], nus_data[14], atol=1e-7)
    assert np.allclose(full_data[3, 7, 5], nus_data[15], atol=1e-7)


def test_reorder_nus_rev_aqorder():
    """reorder 3d nus data"""
    nus_data = np.linspace(0, 1, 128 * 4 * 2 * 2).reshape(-1, 128)
    nuslist = [(0, 0), (2, 3), (5, 3), (6, 7)]
    full_shape = (16, 16, 128)
    full_data = ng.proc_base.expand_nus(
        data=nus_data, shape=full_shape, nuslist=nuslist, aqorder=[0, 1]
    )

    assert full_data.shape == full_shape
    assert np.allclose(full_data[0, 0], nus_data[0], atol=1e-7)
    assert np.allclose(full_data[0, 1], nus_data[1], atol=1e-7)
    assert np.allclose(full_data[1, 0], nus_data[2], atol=1e-7)
    assert np.allclose(full_data[1, 1], nus_data[3], atol=1e-7)

    assert np.allclose(full_data[4, 6], nus_data[4], atol=1e-7)
    assert np.allclose(full_data[4, 7], nus_data[5], atol=1e-7)
    assert np.allclose(full_data[5, 6], nus_data[6], atol=1e-7)
    assert np.allclose(full_data[5, 7], nus_data[7], atol=1e-7)

    assert np.allclose(full_data[10, 6], nus_data[8], atol=1e-7)
    assert np.allclose(full_data[10, 7], nus_data[9], atol=1e-7)
    assert np.allclose(full_data[11, 6], nus_data[10], atol=1e-7)
    assert np.allclose(full_data[11, 7], nus_data[11], atol=1e-7)

    assert np.allclose(full_data[12, 14], nus_data[12], atol=1e-7)
    assert np.allclose(full_data[12, 15], nus_data[13], atol=1e-7)
    assert np.allclose(full_data[13, 14], nus_data[14], atol=1e-7)
    assert np.allclose(full_data[13, 15], nus_data[15], atol=1e-7)


def test_reorder_nus_3d_quadorder():
    """reorder 3d nus data"""
    nus_data = np.linspace(0, 1, 128 * 4 * 2 * 2).reshape(-1, 128)
    nuslist = [(0, 0), (2, 3), (5, 3), (6, 7)]
    full_shape = (16, 16, 128)
    full_data = ng.proc_base.expand_nus(
        data=nus_data,
        shape=full_shape,
        nuslist=nuslist,
        quadrature_order=[(0, 0), (1, 0), (0, 1), (1, 1)],
    )

    assert full_data.shape == full_shape
    assert np.allclose(full_data[0, 0], nus_data[0], atol=1e-7)
    assert np.allclose(full_data[1, 0], nus_data[1], atol=1e-7)
    assert np.allclose(full_data[0, 1], nus_data[2], atol=1e-7)
    assert np.allclose(full_data[1, 1], nus_data[3], atol=1e-7)

    assert np.allclose(full_data[6, 4], nus_data[4], atol=1e-7)
    assert np.allclose(full_data[7, 4], nus_data[5], atol=1e-7)
    assert np.allclose(full_data[6, 5], nus_data[6], atol=1e-7)
    assert np.allclose(full_data[7, 5], nus_data[7], atol=1e-7)

    assert np.allclose(full_data[6, 10], nus_data[8], atol=1e-7)
    assert np.allclose(full_data[7, 10], nus_data[9], atol=1e-7)
    assert np.allclose(full_data[6, 11], nus_data[10], atol=1e-7)
    assert np.allclose(full_data[7, 11], nus_data[11], atol=1e-7)

    assert np.allclose(full_data[14, 12], nus_data[12], atol=1e-7)
    assert np.allclose(full_data[15, 12], nus_data[13], atol=1e-7)
    assert np.allclose(full_data[14, 13], nus_data[14], atol=1e-7)
    assert np.allclose(full_data[15, 13], nus_data[15], atol=1e-7)


def _reference_hadamard(data):
    """Reference implementation matching the historic Hadamard pipeline."""
    order = int(np.ceil(np.log2(data.shape[-1])))
    max_size = 2 ** order
    pad = max_size - data.shape[-1]
    zdata = ng.proc_base.zf(data, pad)

    nat = np.array(np.dot(zdata, scipy.linalg.hadamard(max_size)),
                   dtype=data.dtype)
    bit_reversal = [ng.proc_base.bin2int(
        ng.proc_base.int2bin(value, digits=order)[::-1]
    ) for value in range(max_size)]
    gray = ng.proc_base.gray(order)
    return np.take(np.take(nat, bit_reversal, axis=-1), gray, axis=-1)


def test_hadamard_matches_reference_for_complex_input():
    """The FWHT path should match the previous Hadamard implementation."""
    real = np.arange(12, dtype=np.float32).reshape(2, 6)
    imag = np.arange(12, 24, dtype=np.float32).reshape(2, 6)
    data = real + 1j * imag

    expected = _reference_hadamard(data)
    result = ng.proc_base.ha(data)

    assert result.dtype == data.dtype
    np.testing.assert_allclose(result, expected, atol=1e-6)


def test_hilbert_default_matches_scipy_last_axis():
    """The default N should use the trace length and preserve real values."""
    data = np.arange(12, dtype=np.float32).reshape(3, 4)

    expected = scipy.signal.hilbert(data.real, N=data.shape[-1], axis=-1)
    expected = np.asarray(expected, dtype=np.complex64)
    expected.real = data.real

    result = ng.proc_base.ht(data)

    assert result.dtype == expected.dtype
    np.testing.assert_allclose(result, expected, atol=1e-6)


def test_hilbert_respects_requested_fourier_components():
    """The Hilbert transform should match SciPy for explicit zero filling."""
    data = np.linspace(0, 1, 10, dtype=np.float32).reshape(2, 5)
    n_components = 8

    expected = scipy.signal.hilbert(data.real, N=n_components, axis=-1)
    expected = np.asarray(expected[..., :data.shape[-1]] *
                          (n_components / data.shape[-1]),
                          dtype=np.complex64)
    expected.real = data.real

    result = ng.proc_base.ht(data, N=n_components)

    np.testing.assert_allclose(result, expected, atol=1e-6)

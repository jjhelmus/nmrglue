import numpy as np

import nmrglue as ng


def test_find_lpc_qr():
    D = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=np.float64,
    )
    d = np.array([[1.0], [2.0], [3.0]], dtype=np.float64)

    result = np.asarray(ng.proc_lp.find_lpc_qr(D, d))

    assert result.shape == (2, 1)
    assert result.dtype == np.dtype("float64")
    np.testing.assert_allclose(result, [[1.0], [2.0]], rtol=1e-14, atol=1e-14)

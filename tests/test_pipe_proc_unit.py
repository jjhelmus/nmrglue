import numpy as np

import nmrglue as ng


def test_ht_ps90_180_even_length():
    data = np.array([1.0, 2.0, 4.0, 8.0, 16.0, 32.0], dtype="float32")
    expected = np.array(
        [
            1.0 - 25.333336j,
            2.0 - 7.333332j,
            4.0 - 14.666667j,
            8.0 - 7.6666665j,
            16.0 - 15.333332j,
            32.0 + 12.666667j,
        ],
        dtype="complex64",
    )

    _, result = ng.pipe_proc.ht(
        ng.pipe.create_empty_dic(), data, mode="ps90-180"
    )

    assert result.shape == data.shape
    assert result.dtype == expected.dtype
    np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-6)


def test_ht_ps90_180_odd_length():
    data = np.array([1.0, 2.0, 4.0, 8.0, 16.0], dtype="float32")
    expected = np.array(
        [
            1.0 - 12.242127j,
            2.0 - 5.3340144j,
            4.0 - 5.8728495j,
            8.0 - 6.950515j,
            16.0 + 5.830447j,
        ],
        dtype="complex64",
    )

    _, result = ng.pipe_proc.ht(
        ng.pipe.create_empty_dic(), data, mode="ps90-180"
    )

    assert result.shape == data.shape
    assert result.dtype == expected.dtype
    np.testing.assert_allclose(result, expected, rtol=1e-6, atol=1e-6)

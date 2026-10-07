"""
tests/test_autophase.py — Tests for nmrglue.process.proc_autophase

Covers autops() with and without p1_bounds, both built-in scoring functions,
callable scoring functions, return_phases flag, and error handling.

The p1-bounded path uses scipy.optimize.minimize (Nelder-Mead); the unbounded
path uses scipy.optimize.fmin. Both are exercised here.
"""
import numpy as np
from numpy.testing import assert_allclose
import pytest

import nmrglue as ng
from nmrglue.process.proc_autophase import (
    autops,
    _ps_acme_score,
    _ps_peak_minima_score,
)
from nmrglue.process.proc_base import ps

# Tolerance (degrees) when comparing recovered phases with the applied error.
PHASE_TOL = 5.0


# =============================================================================
# Helpers
# =============================================================================

def _lorentzian_spectrum(n=512, sw=10000.0, f0=1000.0, t2=0.01,
                         p0=0.0, p1=0.0):
    """Return a phase-shifted Lorentzian spectrum (complex ndarray).

    Parameters
    ----------
    n   : int   — number of points
    sw  : float — spectral width (Hz)
    f0  : float — peak frequency (Hz)
    t2  : float — transverse relaxation time (s)
    p0  : float — zero-order phase error (degrees)
    p1  : float — first-order phase error (degrees)
    """
    t = np.arange(n) / sw
    fid = np.exp(1j * 2 * np.pi * f0 * t) * np.exp(-t / t2)
    # Halve the first point so the DFT does not add a baseline offset that
    # ACME would compensate with a spurious p1.
    fid[0] *= 0.5
    spec = np.fft.fftshift(np.fft.fft(fid))
    return ps(spec, p0=p0, p1=p1)


# =============================================================================
# Backward-compatibility — unbounded (original fmin path)
# =============================================================================

def test_autops_acme_returns_ndarray():
    """autops with 'acme' returns an ndarray of the same shape."""
    data = _lorentzian_spectrum(p0=30.0)
    result = autops(data, 'acme', disp=False)
    assert isinstance(result, np.ndarray)
    assert result.shape == data.shape


def test_autops_peak_minima_returns_ndarray():
    """autops with 'peak_minima' returns an ndarray of the same shape."""
    data = _lorentzian_spectrum(p0=30.0)
    result = autops(data, 'peak_minima', disp=False)
    assert isinstance(result, np.ndarray)
    assert result.shape == data.shape


def test_autops_return_phases_unbounded():
    """return_phases=True returns (ndarray, array-like) with two elements."""
    data = _lorentzian_spectrum(p0=45.0)
    result = autops(data, 'acme', return_phases=True, disp=False)
    assert isinstance(result, tuple) and len(result) == 2
    phased, opt = result
    assert isinstance(phased, np.ndarray)
    assert len(opt) == 2


def test_autops_callable_fn_unbounded():
    """A callable scoring function is accepted in place of a string alias."""
    data = _lorentzian_spectrum(p0=20.0)
    result = autops(data, _ps_acme_score, disp=False)
    assert isinstance(result, np.ndarray)


def test_autops_unknown_fn_raises():
    """An unrecognised string fn raises KeyError with a helpful message."""
    data = _lorentzian_spectrum()
    with pytest.raises(KeyError, match='Unable to find algorithm'):
        autops(data, 'not_a_real_algorithm')


def test_autops_acme_reduces_phase_error():
    """autops 'acme' reduces a known phase error towards zero."""
    p0_true = 60.0
    data = _lorentzian_spectrum(p0=p0_true)
    phased, opt = autops(data, 'acme', return_phases=True, disp=False)
    # The correction must undo the applied error: opt[0] ~ -p0_true, p1 ~ 0.
    assert abs(p0_true + opt[0]) < PHASE_TOL, f"p0={opt[0]}"
    assert abs(opt[1]) < PHASE_TOL, f"p1={opt[1]}"


# =============================================================================
# p1-bounded path — scipy.optimize.minimize(Nelder-Mead)
# =============================================================================

def _wrapped_error(p0_true, p0_found):
    """Residual zero-order phase error wrapped into [-180, 180)."""
    return (p0_true + p0_found + 180) % 360 - 180


def test_autops_p1_bounded_returns_ndarray():
    """autops with p1_bounds= returns an ndarray of the same shape."""
    data = _lorentzian_spectrum(p0=30.0)
    result = autops(data, 'acme', p1_bounds=(-1800, 1800))
    assert isinstance(result, np.ndarray)
    assert result.shape == data.shape


def test_autops_p1_bounded_return_phases():
    """return_phases=True with p1_bounds returns (ndarray, array-like)."""
    data = _lorentzian_spectrum(p0=45.0)
    phased, opt = autops(data, 'acme', p1_bounds=(-1800, 1800),
                         return_phases=True)
    assert isinstance(phased, np.ndarray)
    assert len(opt) == 2


def test_autops_p1_fixed_to_zero():
    """p1_bounds=(0, 0) keeps the first-order phase exactly zero."""
    data = _lorentzian_spectrum(p0=50.0)
    _, opt = autops(data, 'acme', p1_bounds=(0, 0), return_phases=True)
    assert opt[1] == 0.0, f"p1 should be 0.0, got {opt[1]}"
    assert abs(_wrapped_error(50.0, opt[0])) < PHASE_TOL, f"p0={opt[0]}"


def test_autops_p1_bounded_respects_range():
    """Optimised p1 stays within p1_bounds."""
    data = _lorentzian_spectrum(p1=300.0)
    b1 = (-200, 200)
    _, opt = autops(data, 'acme', p1_bounds=b1, return_phases=True)
    assert b1[0] <= opt[1] <= b1[1], f"p1={opt[1]} outside {b1}"


@pytest.mark.parametrize('p1_bounds', [(None, 1800), (-1800, None),
                                       (None, None)])
def test_autops_p1_bounds_open_side(p1_bounds):
    """None leaves that side of the p1 range unbounded."""
    data = _lorentzian_spectrum(p0=60.0, p1=30.0)
    _, opt = autops(data, 'acme', p1_bounds=p1_bounds, return_phases=True)
    assert abs(_wrapped_error(60.0, opt[0])) < PHASE_TOL, f"p0={opt[0]}"
    assert abs(30.0 + opt[1]) < PHASE_TOL, f"p1={opt[1]}"


def test_autops_p1_bounded_acme_reduces_phase_error():
    """p1-bounded autops 'acme' undoes a known phase error."""
    p0_true, p1_true = 60.0, 30.0
    data = _lorentzian_spectrum(p0=p0_true, p1=p1_true)
    _, opt = autops(data, 'acme', p1_bounds=(-1800, 1800),
                    return_phases=True)
    assert abs(p0_true + opt[0]) < PHASE_TOL, f"p0={opt[0]}"
    assert abs(p1_true + opt[1]) < PHASE_TOL, f"p1={opt[1]}"


@pytest.mark.parametrize('p0_true', range(-175, 180, 5))
def test_autops_p1_bounded_recovers_any_p0(p0_true):
    """p0 is recovered over the whole circle, including near +-180."""
    data = _lorentzian_spectrum(p0=p0_true)
    _, opt = autops(data, 'acme', p1_bounds=(-1800, 1800),
                    return_phases=True)
    assert abs(_wrapped_error(p0_true, opt[0])) < PHASE_TOL, f"p0={opt[0]}"
    assert abs(opt[1]) < PHASE_TOL, f"p1={opt[1]}"


@pytest.mark.parametrize('p0_start', [190.0, -530.0])
def test_autops_p1_bounded_p0_wrapped(p0_start):
    """The returned p0 is wrapped into [-180, 180)."""
    data = _lorentzian_spectrum(p0=170.0)
    _, opt = autops(data, 'acme', p0=p0_start, p1_bounds=(-1800, 1800),
                    return_phases=True)
    assert -180 <= opt[0] < 180, f"p0={opt[0]}"
    assert abs(_wrapped_error(170.0, opt[0])) < PHASE_TOL, f"p0={opt[0]}"


def test_autops_p1_bounded_peak_minima():
    """p1-bounded autops works with 'peak_minima' scoring function."""
    data = _lorentzian_spectrum(p0=40.0)
    result = autops(data, 'peak_minima', p1_bounds=(-1800, 1800))
    assert isinstance(result, np.ndarray)
    assert result.shape == data.shape


def test_autops_p1_bounded_callable_fn():
    """p1-bounded autops accepts a callable scoring function."""
    data = _lorentzian_spectrum(p0=25.0)
    phased, opt = autops(data, _ps_acme_score, p1_bounds=(-1800, 1800),
                         return_phases=True)
    assert isinstance(phased, np.ndarray)
    assert len(opt) == 2


def test_autops_p1_bounded_unknown_fn_raises():
    """An unrecognised string fn raises KeyError even with p1_bounds."""
    data = _lorentzian_spectrum()
    with pytest.raises(KeyError, match='Unable to find algorithm'):
        autops(data, 'not_a_real_algorithm', p1_bounds=(-1800, 1800))


@pytest.mark.parametrize('p1_bounds', [
    (0,),
    (-1800, 0, 1800),
    [(-180, 180), (-1800, 1800)],
])
def test_autops_invalid_p1_bounds_raises(p1_bounds):
    """p1_bounds must be a single (min, max) pair."""
    data = _lorentzian_spectrum(p0=30.0)
    with pytest.raises(ValueError, match='p1_bounds must be a'):
        autops(data, 'acme', p1_bounds=p1_bounds)


def test_autops_p1_bounded_kwargs_passed_to_minimize():
    """Extra kwargs are forwarded to scipy.optimize.minimize."""
    data = _lorentzian_spectrum(p0=30.0)
    # options= is a minimize-specific kwarg; should not raise
    phased, opt = autops(data, 'acme', p1_bounds=(-1800, 1800),
                         return_phases=True,
                         options={'xatol': 1e-2, 'fatol': 1e-2, 'disp': False})
    assert isinstance(phased, np.ndarray)


def test_autops_p1_bounded_accepts_fmin_kwargs():
    """fmin-style kwargs (disp, xtol, ftol, ...) also work with p1_bounds."""
    data = _lorentzian_spectrum(p0=30.0)
    _, opt = autops(data, 'acme', p1_bounds=(-1800, 1800),
                    return_phases=True,
                    disp=False, xtol=1e-4, ftol=1e-4, maxiter=2000)
    assert abs(30.0 + opt[0]) < PHASE_TOL, f"p0={opt[0]}"


def test_autops_p1_bounded_vs_unbounded_agreement():
    """p1-bounded and unbounded results agree on a clean spectrum."""
    data = _lorentzian_spectrum(p0=30.0, p1=20.0)
    _, opt_ub = autops(data, 'acme', return_phases=True, disp=False)
    _, opt_b = autops(data, 'acme', p1_bounds=(-1800, 1800),
                      return_phases=True)
    assert_allclose(opt_b, opt_ub, atol=1.0)


# =============================================================================
# Score function unit tests
# =============================================================================

def test_acme_score_returns_scalar():
    """_ps_acme_score returns a finite scalar."""
    data = _lorentzian_spectrum()
    score = _ps_acme_score((0.0, 0.0), data)
    assert np.isfinite(score)
    assert np.isscalar(score) or score.ndim == 0


def test_acme_score_lower_for_phased():
    """_ps_acme_score is lower for correctly phased than misphased data."""
    data = _lorentzian_spectrum()
    score_good = _ps_acme_score((0.0, 0.0), data)
    score_bad  = _ps_acme_score((90.0, 0.0), data)
    assert score_good < score_bad


def test_peak_minima_score_returns_scalar():
    """_ps_peak_minima_score returns a finite scalar."""
    data = _lorentzian_spectrum()
    score = _ps_peak_minima_score((0.0, 0.0), data, peak_width=50)
    assert np.isfinite(score)


if __name__ == '__main__':
    import pytest as _pytest
    _pytest.main([__file__, '-v'])

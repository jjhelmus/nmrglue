"""
tests/test_autophase.py — Tests for nmrglue.process.proc_autophase

Covers autops() with and without bounds, both built-in scoring functions,
callable scoring functions, return_phases flag, and error handling.

The bounded path uses scipy.optimize.minimize (Nelder-Mead); the unbounded
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


# =============================================================================
# Helpers
# =============================================================================

def _lorentzian_spectrum(n=512, sw=10000.0, f0=1000.0, t2=0.1, p0=0.0, p1=0.0):
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
    # After correction the residual phase should be much smaller than the
    # original error (we don't require exact recovery — ACME is heuristic).
    assert abs(opt[0]) < abs(p0_true)


# =============================================================================
# Bounded path — scipy.optimize.minimize(Nelder-Mead)
# =============================================================================

def test_autops_bounded_returns_ndarray():
    """autops with bounds= returns an ndarray of the same shape."""
    data = _lorentzian_spectrum(p0=30.0)
    result = autops(data, 'acme', bounds=[(-180, 180), (-1800, 1800)])
    assert isinstance(result, np.ndarray)
    assert result.shape == data.shape


def test_autops_bounded_return_phases():
    """return_phases=True with bounds returns (ndarray, array-like)."""
    data = _lorentzian_spectrum(p0=45.0)
    phased, opt = autops(data, 'acme',
                         bounds=[(-180, 180), (-1800, 1800)],
                         return_phases=True)
    assert isinstance(phased, np.ndarray)
    assert len(opt) == 2


def test_autops_p1_fixed_to_zero():
    """Bounds (0, 0) on p1 keep the first-order phase exactly zero."""
    data = _lorentzian_spectrum(p0=50.0, p1=0.0)
    _, opt = autops(data, 'acme',
                    bounds=[(-180, 180), (0, 0)],
                    return_phases=True)
    assert abs(opt[1]) < 1e-6, f"p1 should be 0.0, got {opt[1]}"


def test_autops_p0_fixed_to_zero():
    """Bounds (0, 0) on p0 keep the zero-order phase exactly zero."""
    data = _lorentzian_spectrum(p0=0.0, p1=200.0)
    _, opt = autops(data, 'acme',
                    bounds=[(0, 0), (-1800, 1800)],
                    return_phases=True)
    assert abs(opt[0]) < 1e-6, f"p0 should be 0.0, got {opt[0]}"


def test_autops_bounded_respects_p0_range():
    """Optimised p0 stays within the specified bounds."""
    data = _lorentzian_spectrum(p0=30.0)
    b0 = (-90, 90)
    _, opt = autops(data, 'acme',
                    bounds=[b0, (-1800, 1800)],
                    return_phases=True)
    assert b0[0] <= opt[0] <= b0[1], f"p0={opt[0]} outside {b0}"


def test_autops_bounded_respects_p1_range():
    """Optimised p1 stays within the specified bounds."""
    data = _lorentzian_spectrum(p0=0.0, p1=300.0)
    b1 = (-500, 500)
    _, opt = autops(data, 'acme',
                    bounds=[(-180, 180), b1],
                    return_phases=True)
    assert b1[0] <= opt[1] <= b1[1], f"p1={opt[1]} outside {b1}"


def test_autops_bounded_acme_reduces_phase_error():
    """Bounded autops 'acme' reduces a known phase error towards zero."""
    p0_true = 60.0
    data = _lorentzian_spectrum(p0=p0_true)
    _, opt = autops(data, 'acme',
                    bounds=[(-180, 180), (-1800, 1800)],
                    return_phases=True)
    assert abs(opt[0]) < abs(p0_true)


def test_autops_bounded_peak_minima():
    """Bounded autops works with 'peak_minima' scoring function."""
    data = _lorentzian_spectrum(p0=40.0)
    result = autops(data, 'peak_minima',
                    bounds=[(-180, 180), (-1800, 1800)])
    assert isinstance(result, np.ndarray)
    assert result.shape == data.shape


def test_autops_bounded_callable_fn():
    """Bounded autops accepts a callable scoring function."""
    data = _lorentzian_spectrum(p0=25.0)
    phased, opt = autops(data, _ps_acme_score,
                         bounds=[(-180, 180), (-1800, 1800)],
                         return_phases=True)
    assert isinstance(phased, np.ndarray)
    assert len(opt) == 2


def test_autops_bounded_unknown_fn_raises():
    """An unrecognised string fn raises KeyError even when bounds are set."""
    data = _lorentzian_spectrum()
    with pytest.raises(KeyError, match='Unable to find algorithm'):
        autops(data, 'not_a_real_algorithm',
               bounds=[(-180, 180), (-1800, 1800)])


def test_autops_bounded_kwargs_passed_to_minimize():
    """Extra kwargs are forwarded to scipy.optimize.minimize."""
    data = _lorentzian_spectrum(p0=30.0)
    # options= is a minimize-specific kwarg; should not raise
    phased, opt = autops(data, 'acme',
                         bounds=[(-180, 180), (-1800, 1800)],
                         return_phases=True,
                         options={'xatol': 1e-2, 'fatol': 1e-2, 'disp': False})
    assert isinstance(phased, np.ndarray)


def test_autops_bounded_vs_unbounded_agreement():
    """Bounded and unbounded results are consistent on a well-conditioned spectrum."""
    data = _lorentzian_spectrum(p0=30.0, p1=0.0)
    _, opt_ub = autops(data, 'acme', return_phases=True, disp=False)
    _, opt_b  = autops(data, 'acme',
                       bounds=[(-180, 180), (-1800, 1800)],
                       return_phases=True)
    # Both should find similar p0 (within solver tolerance + ACME heuristic noise)
    assert abs(opt_ub[0] - opt_b[0]) < 10.0, (
        f"Unbounded p0={opt_ub[0]:.2f} and bounded p0={opt_b[0]:.2f} disagree")


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

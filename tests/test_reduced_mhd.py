"""Validation suite for plasmaplay.reduced_mhd (2-D Strauss reduced MHD).

Falsifiable checks for the linear tearing phase of the nonlinear solver:
  * the elliptic inversion lap(phi) = U is exact (to round-off);
  * a seeded mode GROWS for k a < 1 (Delta' > 0) and DECAYS for k a > 1;
  * the linear growth rate obeys the Furth-Killeen-Rosenbluth gamma ~ S^(-3/5) law
    (the same resistive-layer scaling validated by the slab eigenvalue in T4).

The nonlinear island *saturation* is exercised separately in the gif/experiment; the
absolute growth coefficient differs from the T4 eigenvalue by an O(1) convention
factor (documented), so the tests assert the *scaling*, not the absolute value.
"""

import numpy as np
import pytest

from plasmaplay import reduced_mhd as rm


def _growth_rate(k, S, *, nx=224, ny=16, Lx=4.0, t_end=55.0, dt=0.008):
    """Median late-time growth rate of the m=1 reconnected flux (Alfven units)."""
    s = rm.ReducedMHD(k, S=S, Pm=0.0, nx=nx, ny=ny, Lx=Lx).seed(1e-6)
    ts, amps = [], []
    for i in range(int(t_end / dt)):
        s.step(dt)
        if i % 40 == 0:
            ts.append(s.t)
            amps.append(s.mode_amplitude())
    ts, amps = np.array(ts), np.array(amps)
    g = np.gradient(np.log(amps), ts)
    return float(np.median(g[len(g) // 2:]))


@pytest.fixture(scope="module")
def growth():
    """The three growth rates the linear tests share, computed once per module."""
    return {
        ("u", 0.5, 400.0): _growth_rate(0.5, 400.0),
        ("u", 0.5, 1600.0): _growth_rate(0.5, 1600.0),
        ("s", 1.4, 400.0): _growth_rate(1.4, 400.0, t_end=35.0),
    }


def test_elliptic_inversion_is_exact():
    """lap(phi) = U recovers phi (phi = 0 at the x-walls, periodic in y) to round-off."""
    s = rm.ReducedMHD(0.5, nx=128, ny=32, Lx=6.0)
    phi_true = (np.sin(np.pi * (s.x[:, None] + s.Lx) / (2 * s.Lx))
                * np.cos(s.k * s.y)[None, :])
    U = s._lap(phi_true)
    phi_rec = s.invert_phi(U)
    err = np.max(np.abs(phi_rec - phi_true)) / np.max(np.abs(phi_true))
    assert err < 1e-10


def test_tearing_grows_below_threshold_decays_above(growth):
    """A seeded mode grows for k a < 1 (Delta' > 0) and decays for k a > 1."""
    assert growth[("u", 0.5, 400.0)] > 0.0         # k a = 0.5 < 1: unstable
    assert growth[("s", 1.4, 400.0)] < 0.0         # k a = 1.4 > 1: stable


def test_linear_growth_scales_as_S_minus_three_fifths(growth):
    """gamma ~ S^(-3/5): a factor-4 rise in Lundquist number drops gamma by ~4^(-0.6)
    (the FKR resistive-layer law, reproduced by direct simulation)."""
    g_lo = growth[("u", 0.5, 400.0)]
    g_hi = growth[("u", 0.5, 1600.0)]
    exponent = np.log(g_hi / g_lo) / np.log(1600.0 / 400.0)
    assert -0.75 < exponent < -0.45                # FKR -0.6, with discretization slack

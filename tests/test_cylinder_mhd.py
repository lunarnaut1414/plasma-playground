"""Validation suite for plasmaplay.cylinder_mhd (straight-tokamak linear MHD).

Falsifiable checks:
  * the screw-pinch q-profile has the analytic on-axis / edge values and is monotone;
  * the q = m/n rational surface is located where q actually equals m/n;
  * the m=1 internal kink is unstable exactly when q(0) < 1 (the sawtooth trigger);
  * the sign of the cylindrical Delta' predicts tearing stability, and Delta'
    decreases (more stable) with mode number m;
  * the resistive growth rate follows the FKR gamma ~ S^(-3/5) layer law.
"""

import numpy as np
import pytest

from plasmaplay import cylinder_mhd as cm


# --- screw-pinch q-profile -------------------------------------------------
@pytest.mark.parametrize("q0, nu", [(0.8, 1.0), (1.1, 1.0), (0.9, 2.0)])
def test_q_profile_axis_and_edge(q0, nu):
    """q(0) = q0 and q(a) = (nu+1) q0 for the (1-r^2)^nu current profile."""
    assert cm.screw_pinch_q(0.0, q0, nu) == pytest.approx(q0, rel=1e-6)
    assert cm.screw_pinch_q(1.0, q0, nu) == pytest.approx((nu + 1.0) * q0, rel=1e-6)


def test_q_profile_monotonic():
    r = np.linspace(0.0, 1.0, 200)
    q = cm.screw_pinch_q(r, 0.9, 1.0)
    assert np.all(np.diff(q) > 0)               # q rises monotonically outward


# --- rational surfaces -----------------------------------------------------
def test_rational_surface_is_resonant():
    """q at the located surface equals m/n; absent when m/n is out of range."""
    r_s = cm.rational_surface(2, 1, 1.2, 2.0)   # q from 1.2 to 3.6 -> q=2 exists
    assert r_s is not None
    assert cm.screw_pinch_q(r_s, 1.2, 2.0) == pytest.approx(2.0, abs=1e-4)
    assert cm.rational_surface(1, 1, 1.5, 1.0) is None   # q0=1.5 > 1: no q=1 surface


# --- the m=1 internal kink: the sawtooth trigger ---------------------------
@pytest.mark.parametrize("q0, unstable", [(0.7, True), (0.85, True), (0.95, True),
                                          (1.05, False), (1.3, False)])
def test_internal_kink_threshold_at_q0_unity(q0, unstable):
    """The m=1/n=1 internal kink is unstable iff a q=1 surface exists, i.e. q(0)<1."""
    assert cm.internal_kink_unstable(q0, nu=1.0) is unstable


# --- tearing Delta': sign predicts stability, falls with m -----------------
def test_delta_prime_sign_predicts_stability():
    """Delta' > 0 (unstable) for a broad-current low-m surface; Delta' < 0 (stable)
    for a high-m surface near the edge. The sign is robust to the layer-skip gap."""
    d_unstable = cm.delta_prime_cylinder(2, 1, 1.2, 2.0)
    d_stable = cm.delta_prime_cylinder(3, 1, 0.9, 3.0)
    assert d_unstable > 0
    assert d_stable < 0


def test_delta_prime_decreases_with_mode_number():
    """At a fixed equilibrium, higher-m modes are more wall-stabilized (smaller
    Delta'): Delta'(m=2) > Delta'(m=3) on a profile carrying both surfaces."""
    q0, nu = 1.2, 2.0                            # q: 1.2 -> 3.6, so q=2 and q=3 exist
    assert cm.delta_prime_cylinder(2, 1, q0, nu) > cm.delta_prime_cylinder(3, 1, q0, nu)


# --- resistive growth rate: the FKR S^(-3/5) law ---------------------------
def test_fkr_growth_rate_scales_as_S_minus_three_fifths():
    """gamma ~ S^(-3/5): a decade in Lundquist number drops gamma by 10^(-0.6)."""
    g_lo = cm.fkr_growth_rate(1.5, 1e5)
    g_hi = cm.fkr_growth_rate(1.5, 1e6)
    assert g_hi / g_lo == pytest.approx(10 ** (-0.6), rel=1e-6)


def test_fkr_growth_rate_zero_when_stable():
    """A tearing-stable surface (Delta' <= 0) has no resistive growth."""
    assert cm.fkr_growth_rate(-3.0, 1e5) == 0.0
    assert cm.fkr_growth_rate(2.0, 1e5) > 0.0

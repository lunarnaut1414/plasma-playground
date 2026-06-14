"""Validation suite for plasmaplay.sawtooth (Kadomtsev reconnection — rung B3a).

Falsifiable checks for the reconnection physics (the periodic cycle / period-vs-tau_R
scaling is the follow-on rung B3b):
  * the helical flux psi* peaks exactly at the q = 1 surface;
  * the Kadomtsev mixing radius sits outside the q = 1 surface (and is absent when
    q(0) >= 1, i.e. no crash);
  * the flatten conserves the cross-section integral (thermal energy) exactly and
    leaves a flat core;
  * a single crash reconnects the core: T flattens, the helical flux in the core goes
    to zero, q(0) is reset to ~1, and thermal energy is conserved.
"""

import numpy as np
import pytest

from plasmaplay import sawtooth as st


def test_helical_flux_peaks_at_q1_surface():
    """psi*(r) = int B_theta (1-q) dr' has its maximum exactly at the q = 1 surface
    (where the integrand changes sign)."""
    r = np.linspace(0.0, 1.0, 401)
    q = st.screw_pinch_q(r, 0.8, 1.0)
    psi = st.helical_flux(r, q)
    r_at_max = r[int(np.argmax(psi))]
    i_q1 = int(np.argmin(np.abs(q - 1.0)))
    assert r_at_max == pytest.approx(r[i_q1], abs=2 * (r[1] - r[0]))


def test_mixing_radius_outside_q1_and_absent_when_stable():
    """r_mix lies outside the q = 1 surface and inside the wall; it is None (no crash)
    when q(0) >= 1."""
    r = np.linspace(0.0, 1.0, 401)
    q = st.screw_pinch_q(r, 0.8, 1.0)
    r1 = r[int(np.argmin(np.abs(q - 1.0)))]
    r_mix = st.mixing_radius(r, q)
    assert r1 < r_mix <= 1.0
    assert st.mixing_radius(r, st.screw_pinch_q(r, 1.2, 1.0)) is None


def test_kadomtsev_flatten_conserves_energy_and_flattens():
    """Flattening inside r_mix conserves int field * r dr over [0, r_mix] exactly and
    leaves the core constant."""
    r = np.linspace(0.0, 1.0, 257)
    T = 1.0 - 0.8 * r ** 2
    r_mix = 0.6
    inside = r <= r_mix
    e0 = np.trapezoid(T[inside] * r[inside], r[inside])
    Tf = st.kadomtsev_flatten(r, T, r_mix)
    e1 = np.trapezoid(Tf[inside] * r[inside], r[inside])
    assert e1 == pytest.approx(e0, rel=1e-12)          # energy conserved (exact)
    assert Tf[inside].std() < 1e-12                     # core is flat
    assert np.allclose(Tf[~inside], T[~inside])         # untouched outside


def test_single_crash_reconnects_core():
    """A Kadomtsev crash on a q(0) < 1 state: T flattens inside r_mix, the helical flux
    in the core is reconnected to ~0, q(0) resets to ~1, energy is conserved."""
    sim = st.SawtoothCycle(eta=0.3, nr=129, q0_init=0.82, chi=0.03, heat=0.0,
                           q_reset=1.05)
    r_mix = st.mixing_radius(sim.r, sim.q())
    inside = sim.r <= r_mix
    psi_max_before = st.helical_flux(sim.r, sim.q()).max()
    e0 = np.trapezoid(sim.T[inside] * sim.r[inside], sim.r[inside])
    assert sim.q0() < 1.0                                # unstable before
    sim._crash()
    e1 = np.trapezoid(sim.T[inside] * sim.r[inside], sim.r[inside])
    psi_max_after = st.helical_flux(sim.r, sim.q()).max()
    assert sim.q0() == pytest.approx(1.05, rel=1e-6)    # reconnected toward q = 1
    assert sim.T[inside].std() < 1e-10                  # core temperature flattened
    assert psi_max_after < 0.05 * psi_max_before        # helical flux reconnected
    assert e1 == pytest.approx(e0, rel=1e-10)           # thermal energy conserved

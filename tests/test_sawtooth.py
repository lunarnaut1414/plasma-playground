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


# --- Track C: coupling the sawtooth into a transport burn -------------------
def test_q_from_temperature_peaks_lower_q0():
    """A peaked (hot-core) temperature gives a lower q(0) than a flat one (Spitzer
    current peaks on axis), and q rises monotonically to q_edge at the edge."""
    rho = np.linspace(0.0, 1.0, 129)
    q_peaked = st.q_from_temperature(rho, 1.0 - 0.95 * rho ** 2, q_edge=3.0)
    q_flat = st.q_from_temperature(rho, np.ones_like(rho), q_edge=3.0)
    assert q_peaked[0] < q_flat[0]                       # peaking lowers q(0)
    assert q_peaked[-1] == pytest.approx(3.0, rel=1e-6)  # normalised to q_edge
    assert np.all(np.diff(q_peaked) > -1e-9)             # monotone rising


def test_resistive_relaxation_time_scales_as_T_three_halves():
    """The sawtooth recovery time follows the Spitzer current-diffusion law tau ~ T^1.5:
    it equals tau_ref at T_ref, grows monotonically with T, and a 4x hotter core gives an
    8x longer period (so a burning reactor sawtooths slowly -> 'monster' sawteeth)."""
    assert st.resistive_relaxation_time(20.0, tau_ref=2.0, T_ref_keV=20.0) == pytest.approx(2.0)
    assert st.resistive_relaxation_time(80.0, tau_ref=2.0, T_ref_keV=20.0) == pytest.approx(16.0)
    taus = [st.resistive_relaxation_time(T) for T in (5.0, 10.0, 20.0, 30.0)]
    assert np.all(np.diff(taus) > 0)                     # hotter core -> longer period
    # exact 3/2 power law: doubling T multiplies the period by 2^1.5
    ratio = st.resistive_relaxation_time(40.0) / st.resistive_relaxation_time(20.0)
    assert ratio == pytest.approx(2.0 ** 1.5)


def test_crash_profiles_conserve_particles_and_energy():
    """A crash flattens n and T inside r_mix while conserving the particle content and
    the thermal energy exactly."""
    rho = np.linspace(0.0, 1.0, 129)
    n = 1e20 * (1.0 - 0.5 * rho ** 2)
    T = 20.0 * (1.0 - 0.9 * rho ** 2)
    r_mix = 0.5
    inside = rho <= r_mix
    P0 = np.trapezoid(n[inside] * rho[inside], rho[inside])
    E0 = np.trapezoid(3 * n[inside] * T[inside] * rho[inside], rho[inside])
    n2, T2 = st.crash_profiles(rho, n, T, r_mix)
    P1 = np.trapezoid(n2[inside] * rho[inside], rho[inside])
    E1 = np.trapezoid(3 * n2[inside] * T2[inside] * rho[inside], rho[inside])
    assert P1 == pytest.approx(P0, rel=1e-10)            # particles conserved
    assert E1 == pytest.approx(E0, rel=1e-10)            # thermal energy conserved
    assert T2[inside].std() / T2[inside].mean() < 1e-9   # core flattened


def test_sawtooth_event_fires_only_on_unstable_core():
    """sawtooth_event crashes a hot, kink-unstable core (q(0) < trigger) and leaves a
    cool, kink-stable one untouched (the events-off / Track-A regression)."""
    rho = np.linspace(0.0, 1.0, 129)
    n = np.full_like(rho, 1e20)
    T_hot = 1.0 - 0.97 * rho ** 2                         # strongly peaked -> q(0) low
    _, _, crashed_hot = st.sawtooth_event(rho, n, T_hot, q_edge=2.0)
    assert crashed_hot                                    # the core reconnects
    T_flat = np.full_like(rho, 1.0)                       # flat -> kink-stable
    n_out, T_out, crashed_flat = st.sawtooth_event(rho, n, T_flat, q_edge=2.0)
    assert not crashed_flat                               # no event...
    assert np.array_equal(T_out, T_flat) and np.array_equal(n_out, n)  # ...state intact


def test_stellarator_burn_is_sawtooth_free():
    """The tokamak contrast (Track E2): a hot peaked core that triggers a sawtooth in a
    current-driven tokamak (q from T crosses 1) is sawtooth-FREE in a stellarator, whose
    q is set by the external coils (q > 1 everywhere) and does not respond to the burn —
    so there is no q = 1 surface and no kink, ever (inherently steady-state)."""
    rho = np.linspace(0.0, 1.0, 129)
    n = np.full_like(rho, 1e20)
    T_hot = 1.0 - 0.97 * rho ** 2                         # the same hot, peaked core
    assert st.sawtooth_event(rho, n, T_hot, q_edge=2.0)[2]   # tokamak: crashes
    q_stell = st.external_q_profile(rho)                  # stellarator: q from coils
    assert q_stell.min() > 1.0                            # no q = 1 surface...
    assert st.mixing_radius(rho, q_stell) is None         # ...so no kink / reconnection

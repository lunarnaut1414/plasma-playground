"""Validation suite for plasmaplay.transport (burning-plasma transport).

Falsifiable checks, in the spirit of FUNDAMENTALS.md:
  * Bosch-Hale D-T reactivity reproduces published <sigma v> values.
  * Fusion power partitions 3.5/14.1 MeV correctly; bremsstrahlung scales n^2 sqrt(T).
  * 0-D burn IGNITES above the Lawson triple product and DIES below it.
  * At a burning steady state, alpha self-heating balances the losses.
  * The 1-D solver heats a peaked core and fuelling raises the density.
"""

import numpy as np
import pytest

from plasmaplay import transport as tr


# --- Bosch-Hale reactivity ------------------------------------------------
@pytest.mark.parametrize("T_keV, expected", [
    (10.0, 1.1e-22),    # published <sigma v>_DT ~ 1.1e-22 m^3/s
    (20.0, 4.2e-22),
    (64.0, 8.7e-22),    # near the peak
])
def test_reactivity_published(T_keV, expected):
    val = tr.reactivity_dt(T_keV)
    assert val == pytest.approx(expected, rel=0.05)


def test_reactivity_monotonic_below_peak():
    T = np.linspace(2, 50, 50)
    sv = tr.reactivity_dt(T)
    assert np.all(np.diff(sv) > 0)            # rises monotonically up to ~64 keV


def test_reactivity_peaks_near_64keV():
    T = np.linspace(10, 120, 400)
    Tpeak = T[np.argmax(tr.reactivity_dt(T))]
    assert 55.0 < Tpeak < 75.0


# --- power densities ------------------------------------------------------
def test_alpha_is_one_fifth_of_fusion():
    n, T = 1e20, 15.0
    ratio = tr.fusion_power_density(n, T, "alpha") / tr.fusion_power_density(n, T, "total")
    assert ratio == pytest.approx(3.5 / 17.6, rel=1e-6)


def test_fusion_power_scales_density_squared():
    p1 = tr.fusion_power_density(1e20, 15.0, "total")
    p2 = tr.fusion_power_density(2e20, 15.0, "total")
    assert p2 / p1 == pytest.approx(4.0, rel=1e-6)


def test_bremsstrahlung_scaling():
    base = tr.bremsstrahlung_density(1e20, 10.0, z_eff=1.0)
    assert tr.bremsstrahlung_density(2e20, 10.0) / base == pytest.approx(4.0, rel=1e-9)
    assert tr.bremsstrahlung_density(1e20, 40.0) / base == pytest.approx(2.0, rel=1e-9)
    assert tr.bremsstrahlung_density(1e20, 10.0, z_eff=2.0) / base == pytest.approx(2.0)


# --- 0-D burn: ignition is a threshold in the triple product --------------
def _kick(t):
    """A gentle auxiliary-heating kick that switches off at t = 4 s."""
    return 3e5 if t < 4.0 else 0.0


def test_0d_ignites_above_lawson():
    """Density held at 1e20, modest heating off at t=4s -> stays lit (Q->inf)."""
    n0 = 1e20
    tau_E = 2.0
    fuel = n0 / (3 * tau_E)                    # hold density
    r = tr.burn_0d(n0, 2.0, tau_E=tau_E, p_aux=_kick, fuel_rate=fuel, t_end=20.0)
    assert r["p_aux"][-1] == 0.0               # heating is off at the end
    assert r["T"][-1] > 10.0                   # and the plasma is still hot -> ignited


def test_0d_dies_below_lawson():
    """Low density + poor confinement: even with the same kick it cannot stay lit."""
    r = tr.burn_0d(2e19, 2.0, tau_E=0.3, p_aux=_kick, fuel_rate=0.0, t_end=20.0)
    assert r["T"][-1] < 3.0                     # collapses once the heating stops


def test_0d_steady_state_power_balance():
    """At the burning steady state alpha heating balances transport + radiation."""
    n0 = 1e20
    tau_E = 2.0
    fuel = n0 / (3 * tau_E)
    r = tr.burn_0d(n0, 2.0, tau_E=tau_E, p_aux=_kick, fuel_rate=fuel, t_end=40.0)
    p_alpha = r["p_alpha"][-1]
    p_loss = r["W"][-1] / tau_E + r["p_brem"][-1]
    assert p_alpha == pytest.approx(p_loss, rel=0.02)


def test_0d_triple_product_at_ignition_near_lawson():
    """The triple product when self-heating first matches losses is order 1e21."""
    n0 = 1e20
    tau_E = 2.0
    fuel = n0 / (3 * tau_E)
    r = tr.burn_0d(n0, 2.0, tau_E=tau_E, p_aux=_kick, fuel_rate=fuel, t_end=20.0)
    # find where alpha first overtakes the loss term (ignition crossover)
    p_loss = r["W"] / tau_E + r["p_brem"]
    cross = np.where(r["p_alpha"] > p_loss)[0]
    assert cross.size > 0
    triple = r["triple"][cross[0]]
    assert 1e21 < triple < 2e22                # Lawson scale (~3e21 for DT)


# --- 1-D transport --------------------------------------------------------
def test_gaussian_deposition_normalized():
    rho = np.linspace(0, 1, 257)
    g = tr.gaussian_deposition(rho, 0.0, 0.3)
    assert np.trapezoid(g * 2.0 * rho, rho) == pytest.approx(1.0, rel=1e-3)


def test_1d_heating_builds_peaked_core():
    sim = tr.Transport1D(a=1.0, n_grid=65, chi=1.5, D=0.5)
    sim.set_state(T=2.0, n=1e20)
    for _ in range(1000):
        sim.step(2e-3, p_aux_total=3e6)
    assert sim.T[0] > 20.0                      # core got hot
    assert sim.T[0] > sim.T[32] > sim.T[-1]     # monotonically peaked profile


def test_1d_fuelling_raises_density():
    sim = tr.Transport1D(a=1.0, n_grid=65, chi=1.5, D=0.5, n_edge=1e19)
    sim.set_state(T=5.0, n=2e19)
    n_before = sim.diagnostics()["n_avg"]
    for _ in range(500):
        sim.step(2e-3, fuel_total=5e20,
                 fuel_profile=tr.gaussian_deposition(sim.rho, 0.0, 0.2))
    assert sim.diagnostics()["n_avg"] > n_before


def test_1d_diffusion_conserves_without_sources():
    """No sources, no fuelling, insulating-ish: energy only leaks out the edge,
    never grows. A blunt check that the implicit diffusion is stable & dissipative."""
    sim = tr.Transport1D(a=1.0, n_grid=65, chi=1.0, D=0.3, T_edge=1.0, n_edge=1e20)
    sim.set_state(T=10.0, n=1e20)
    W0 = sim.diagnostics()["W"]
    for _ in range(300):
        sim.step(2e-3)
    assert sim.diagnostics()["W"] <= W0 + 1e-9  # energy does not spontaneously appear

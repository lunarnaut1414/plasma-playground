"""Validation suite for plasmaplay.operating_limits + the operating-mode burns.

Falsifiable checks:
  * Greenwald density limit reproduces n_G = Ip/(pi a^2) and its scalings.
  * The Martin 2008 L->H power threshold reproduces the ITER value (~50 MW).
  * The L->H and Greenwald confinement multipliers bifurcate / collapse correctly.
  * A 0-D burn reaches H-mode (hot) above the power threshold and stays in L-mode
    (cool) below it.
  * Over-fuelling past the Greenwald limit collapses the burn — and the collapse is
    *reversible* when the fuelling is backed off.
"""

import numpy as np
import pytest

from plasmaplay import operating_limits as ol
from plasmaplay import transport as tr

# A small ITER-like toy device shared by the scenario tests.
R0, A_MIN, B0, IP, KAPPA = 3.0, 1.0, 5.3, 7.0, 1.5
S_PLASMA = 4 * np.pi ** 2 * R0 * A_MIN * np.sqrt((1 + KAPPA ** 2) / 2)
VOL = 2 * np.pi ** 2 * R0 * A_MIN ** 2
N_G = ol.greenwald_density(IP, A_MIN)
BETA_LIM = 0.04


# --- Greenwald density limit ----------------------------------------------
def test_greenwald_density_iter_value():
    """n_G = Ip/(pi a^2): for the ITER baseline (15 MA, a=2 m) ~1.19e20 m^-3."""
    assert ol.greenwald_density(15.0, 2.0) == pytest.approx(1.19e20, rel=0.02)


def test_greenwald_scaling():
    """n_G ∝ Ip and ∝ 1/a^2."""
    base = ol.greenwald_density(10.0, 1.0)
    assert ol.greenwald_density(20.0, 1.0) / base == pytest.approx(2.0)
    assert ol.greenwald_density(10.0, 2.0) / base == pytest.approx(0.25)


# --- L->H power threshold (Martin 2008) ------------------------------------
def test_lh_threshold_iter_scale():
    """The Martin 2008 scaling gives ~50 MW at the ITER operating point
    (n20~0.5, B=5.3 T, S~680 m^2)."""
    p_lh = ol.lh_power_threshold(0.5, 5.3, 680.0)
    assert 35.0 < p_lh < 75.0


def test_lh_threshold_density_monotonic():
    """Higher density raises the L->H threshold (P_LH ∝ n^0.717)."""
    assert ol.lh_power_threshold(1.0, 5.3, 680.0) > ol.lh_power_threshold(0.5, 5.3, 680.0)


# --- confinement multipliers ------------------------------------------------
def test_lh_factor_bifurcates():
    """The L->H multiplier is ~1 well below threshold and ~h_factor well above."""
    assert ol.confinement_factor_lh(2.0, 10.0, h_factor=2.0) == pytest.approx(1.0, abs=0.02)
    assert ol.confinement_factor_lh(40.0, 10.0, h_factor=2.0) == pytest.approx(2.0, abs=0.02)


def test_greenwald_factor_collapses():
    """The density-limit multiplier is ~1 well below n_G and collapses to the floor
    once the density crosses toward n_G."""
    assert ol.confinement_factor_greenwald(0.4 * N_G, N_G) == pytest.approx(1.0, abs=0.02)
    assert ol.confinement_factor_greenwald(1.2 * N_G, N_G, floor=0.08) == pytest.approx(0.08, abs=0.02)


# --- 0-D operating modes ----------------------------------------------------
def _tau_factor(n0):
    """Confinement multiplier = L->H bifurcation x Greenwald collapse for this device."""
    p_lh = ol.lh_power_threshold(n0 / 1e20, B0, S_PLASMA)

    def factor(t, n_e, T, p_heat_density):
        p_heat_MW = p_heat_density * VOL / 1e6
        return (ol.confinement_factor_lh(p_heat_MW, p_lh)
                * ol.confinement_factor_greenwald(n_e, N_G))
    return factor


def test_hmode_hotter_than_lmode():
    """Above the L->H power threshold the burn reaches H-mode (hot); below it the
    same device stays in L-mode (cool)."""
    lmode = tr.burn_0d_ash(5e19, 3.0, tau_E=1.0, p_aux=lambda t: 1.2e5, B=B0,
                           tau_p=6.0, tau_he=10.0, fuel_rate=5e19 / 6,
                           beta_limit=BETA_LIM, tau_factor=_tau_factor(5e19), t_end=45.0)
    hmode = tr.burn_0d_ash(7e19, 3.0, tau_E=1.8, p_aux=lambda t: 3e5, B=B0,
                           tau_p=6.0, tau_he=10.0, fuel_rate=7e19 / 6,
                           beta_limit=BETA_LIM, tau_factor=_tau_factor(7e19), t_end=45.0)
    assert lmode["T"][-1] < 10.0                 # L-mode stays cool
    assert hmode["T"][-1] > 15.0                 # H-mode burns hot
    assert hmode["T"][-1] > 2 * lmode["T"][-1]   # a real confinement jump


def test_density_limit_collapse_is_reversible():
    """Over-fuelling past n_G collapses the burn (the density-limit disruption); back
    the fuelling off and the temperature recovers — the collapse is reversible."""
    def fuel(t):
        return 4.5e19 if 20 <= t < 33 else 7e19 / 6     # over-fuel window, then off

    r = tr.burn_0d_ash(7e19, 3.0, tau_E=1.8, p_aux=lambda t: 3e5, B=B0, tau_p=6.0,
                       tau_he=10.0, fuel_rate=fuel, beta_limit=BETA_LIM,
                       tau_factor=_tau_factor(7e19), t_end=60.0, dt=1e-3)
    t = r["t"]

    def T_at(tt):
        return r["T"][np.argmin(np.abs(t - tt))]

    def ng_at(tt):
        return r["n_e"][np.argmin(np.abs(t - tt))] / N_G

    assert T_at(18) > 15.0                       # burning before over-fuelling
    assert ng_at(31) > 1.0                       # pushed past the Greenwald limit
    assert T_at(31) < 5.0                        # and the burn has collapsed
    assert T_at(58) > 15.0                       # backing off recovers the burn

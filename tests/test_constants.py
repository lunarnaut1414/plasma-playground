"""V0 — formulary sanity checks.

Validate the `plasmaplay.constants` helpers two ways:
  1. against a first-principles SI computation (catches wrong constant / wrong
     power / transcription errors), and
  2. against canonical reference numbers a plasma physicist would recognize.
A third, optional cross-check against PlasmaPy runs only if it is installed.
"""

import math

import pytest
from scipy import constants as sc

from plasmaplay import constants as k

# A representative plasma: dense-ish lab plasma, 1 T field, a few eV.
N = 1.0e18          # number density [m^-3]
B = 1.0             # magnetic field [T]
T_EV = 1.0          # temperature [eV]
V_PERP = 1.0e5      # perpendicular speed [m/s]


# --- 1. First-principles SI agreement -------------------------------------

def test_gyrofrequency_first_principles():
    expected = k.e * B / k.m_e            # electron cyclotron frequency
    assert math.isclose(k.gyrofrequency(k.e, B, k.m_e), expected, rel_tol=1e-12)


def test_gyroradius_first_principles():
    expected = k.m_e * V_PERP / (k.e * B)
    assert math.isclose(k.gyroradius(V_PERP, k.e, B, k.m_e), expected, rel_tol=1e-12)


def test_plasma_frequency_first_principles():
    expected = math.sqrt(N * sc.e**2 / (sc.epsilon_0 * sc.electron_mass))
    assert math.isclose(k.plasma_frequency(N), expected, rel_tol=1e-12)


def test_debye_length_first_principles():
    expected = math.sqrt(sc.epsilon_0 * (T_EV * sc.e) / (N * sc.e**2))
    assert math.isclose(k.debye_length(T_EV, N), expected, rel_tol=1e-12)


def test_alfven_speed_first_principles():
    expected = B / math.sqrt(sc.mu_0 * N * sc.proton_mass)
    assert math.isclose(k.alfven_speed(B, N), expected, rel_tol=1e-12)


def test_thermal_velocity_first_principles():
    expected = math.sqrt(2 * T_EV * sc.e / sc.electron_mass)
    assert math.isclose(k.thermal_velocity(T_EV, sc.electron_mass), expected, rel_tol=1e-12)


# --- 2. Canonical reference numbers ---------------------------------------
# Values every plasma physicist carries around; loose tolerance on purpose.

def test_electron_plasma_frequency_reference():
    # f_pe ~ 8980 * sqrt(n[cm^-3]) Hz ; n = 1e18 m^-3 = 1e12 cm^-3
    f_pe = k.plasma_frequency(1.0e18) / (2 * math.pi)
    assert math.isclose(f_pe, 8.98e9, rel_tol=2e-2)


def test_debye_length_reference():
    # lambda_D ~ 7.43 um for T_e = 1 eV, n = 1e18 m^-3
    assert math.isclose(k.debye_length(1.0, 1.0e18), 7.43e-6, rel_tol=2e-2)


def test_alfven_speed_reference():
    # v_A ~ 6.9e6 m/s for B = 1 T, n = 1e19 m^-3 hydrogen
    assert math.isclose(k.alfven_speed(1.0, 1.0e19), 6.9e6, rel_tol=2e-2)


# --- 3. Optional cross-check vs PlasmaPy ----------------------------------

def test_against_plasmapy():
    plasmapy = pytest.importorskip("plasmapy")
    from astropy import units as u
    from plasmapy.formulary import plasma_frequency, Debye_length, Alfven_speed

    pp_wpe = plasma_frequency(N * u.m**-3, particle="e-").to_value(u.rad / u.s)
    assert math.isclose(k.plasma_frequency(N), pp_wpe, rel_tol=1e-3)

    pp_ld = Debye_length(T_EV * u.eV, N * u.m**-3).to_value(u.m)
    assert math.isclose(k.debye_length(T_EV, N), pp_ld, rel_tol=1e-3)

    pp_va = Alfven_speed(B * u.T, N * u.m**-3, ion="p+").to_value(u.m / u.s)
    assert math.isclose(k.alfven_speed(B, N), pp_va, rel_tol=1e-3)

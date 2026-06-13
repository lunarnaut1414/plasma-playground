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


# --- 2b. Unit-level behavior: aliases, scaling laws, cross-relations ------

def test_constant_aliases():
    assert k.e == k.ELEMENTARY_CHARGE
    assert k.c == k.SPEED_OF_LIGHT
    assert k.m_e == k.ELECTRON_MASS
    assert k.m_p == k.PROTON_MASS
    assert k.m_p > k.m_e        # protons are heavier (sanity)


def test_gyrofrequency_scalings():
    f1 = k.gyrofrequency(k.e, 1.0, k.m_e)
    assert math.isclose(k.gyrofrequency(k.e, 2.0, k.m_e), 2 * f1)   # linear in B
    assert math.isclose(k.gyrofrequency(k.e, 1.0, 2 * k.m_e), f1 / 2)  # inverse in m


def test_gyrofrequency_charge_sign_independent():
    # uses |q|, so electrons and positrons share a gyrofrequency
    assert k.gyrofrequency(k.e, 1.0, k.m_e) == k.gyrofrequency(-k.e, 1.0, k.m_e)


def test_gyroradius_scalings_and_charge_sign():
    r1 = k.gyroradius(1.0e5, k.e, 1.0, k.m_p)
    assert math.isclose(k.gyroradius(2.0e5, k.e, 1.0, k.m_p), 2 * r1)   # linear in v
    assert math.isclose(k.gyroradius(1.0e5, k.e, 2.0, k.m_p), r1 / 2)   # inverse in B
    assert k.gyroradius(1.0e5, -k.e, 1.0, k.m_p) == r1                  # |q|


def test_gyroradius_equals_vperp_over_gyrofrequency():
    v_perp = 3.0e5
    r_L = k.gyroradius(v_perp, k.e, 1.0, k.m_p)
    omega_c = k.gyrofrequency(k.e, 1.0, k.m_p)
    assert math.isclose(r_L, v_perp / omega_c, rel_tol=1e-12)


def test_plasma_frequency_scaling_and_species():
    f1 = k.plasma_frequency(1.0e18)
    assert math.isclose(k.plasma_frequency(4.0e18), 2 * f1)            # sqrt(n)
    # ion plasma frequency is smaller by sqrt(m_e/m_p)
    f_ion = k.plasma_frequency(1.0e18, charge=k.e, mass=k.m_p)
    assert math.isclose(f_ion / f1, math.sqrt(k.m_e / k.m_p), rel_tol=1e-12)


def test_debye_length_scalings():
    d1 = k.debye_length(1.0, 1.0e18)
    assert math.isclose(k.debye_length(4.0, 1.0e18), 2 * d1)          # sqrt(T)
    assert math.isclose(k.debye_length(1.0, 4.0e18), d1 / 2)          # 1/sqrt(n)


def test_alfven_speed_scalings():
    a1 = k.alfven_speed(1.0, 1.0e19)
    assert math.isclose(k.alfven_speed(2.0, 1.0e19), 2 * a1)          # linear in B
    assert math.isclose(k.alfven_speed(1.0, 4.0e19), a1 / 2)          # 1/sqrt(n)


def test_debye_thermal_plasmafreq_relation():
    # lambda_D = v_th / (sqrt(2) * omega_pe) for a single Maxwellian species
    n, T_eV = 1.0e18, 5.0
    lhs = k.debye_length(T_eV, n)
    rhs = k.thermal_velocity(T_eV, k.m_e) / (math.sqrt(2) * k.plasma_frequency(n))
    assert math.isclose(lhs, rhs, rel_tol=1e-12)


# --- 3. Optional cross-check vs PlasmaPy ----------------------------------

def test_against_plasmapy():
    pytest.importorskip("plasmapy")
    from astropy import units as u
    from plasmapy.formulary import plasma_frequency, Debye_length, Alfven_speed

    pp_wpe = plasma_frequency(N * u.m**-3, particle="e-").to_value(u.rad / u.s)
    assert math.isclose(k.plasma_frequency(N), pp_wpe, rel_tol=1e-3)

    pp_ld = Debye_length(T_EV * u.eV, N * u.m**-3).to_value(u.m)
    assert math.isclose(k.debye_length(T_EV, N), pp_ld, rel_tol=1e-3)

    pp_va = Alfven_speed(B * u.T, N * u.m**-3, ion="p+").to_value(u.m / u.s)
    assert math.isclose(k.alfven_speed(B, N), pp_va, rel_tol=1e-3)

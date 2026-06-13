"""Validation for the propulsion kernel (experiment 07).

The channel model obeys two exact relations that make good falsifiable checks:
impulse-momentum (thrust = integrated body force) and energy balance (the work
done by the Lorentz force equals the kinetic-energy gain). Plus the headline MPD
scaling T ∝ I².
"""

import numpy as np

from plasmaplay import propulsion as prop


def test_maecker_scales_as_current_squared():
    t1 = prop.maecker_thrust(1000.0, 0.05, 0.01)
    t2 = prop.maecker_thrust(2000.0, 0.05, 0.01)
    assert np.isclose(t2 / t1, 4.0)              # T ∝ I²


def test_maecker_ballpark():
    # ~24 N at 10 kA with r_a/r_c = 5 — a high-power MPD operating point
    T = prop.maecker_thrust(10_000.0, 0.05, 0.01)
    assert 20.0 < T < 30.0


def test_specific_impulse_and_jet_power():
    T, mdot = 2.0, 1e-4
    isp = prop.specific_impulse(T, mdot)
    assert np.isclose(isp, T / (mdot * prop.G0))
    # jet power = T²/(2 ṁ) = ½ ṁ u_e², with u_e = T/ṁ
    u_e = T / mdot
    assert np.isclose(prop.jet_power(T, mdot), 0.5 * mdot * u_e**2)


def test_channel_thrust_equals_integrated_force():
    # impulse-momentum: ṁ Δu = A ∫ f dx
    n, L = 401, 0.2
    dx = L / (n - 1)
    f = np.full(n, 6.0e3)                        # uniform force density [N/m³]
    area, mdot = 1e-3, 1e-3
    G = mdot / area
    u = prop.channel_velocity(f, dx, G, u_inlet=1000.0)
    thrust = mdot * (u[-1] - u[0])
    integrated_force = area * np.trapezoid(f, dx=dx)
    assert np.isclose(thrust, integrated_force, rtol=1e-6)


def test_channel_energy_balance():
    # Lorentz work ∫ f u A dx  ==  KE gain ½ ṁ (u_exit² − u_inlet²)
    n, L = 801, 0.25
    dx = L / (n - 1)
    x = np.linspace(0, L, n)
    f = 5.0e3 * (1 + np.sin(np.pi * x / L))      # a non-uniform profile
    area, mdot = 1e-3, 1e-3
    G = mdot / area
    u = prop.channel_velocity(f, dx, G, u_inlet=500.0)
    lorentz_work = area * np.trapezoid(f * u, dx=dx)
    ke_gain = 0.5 * mdot * (u[-1] ** 2 - u[0] ** 2)
    assert np.isclose(lorentz_work, ke_gain, rtol=1e-3)


def test_channel_uniform_force_gives_linear_velocity():
    n, L = 201, 0.1
    dx = L / (n - 1)
    f = np.full(n, 1.0e4)
    G = 0.5
    u = prop.channel_velocity(f, dx, G, u_inlet=0.0)
    # du/dx = f/G constant -> u linear with slope f/G
    np.testing.assert_allclose(u, (1.0e4 / G) * np.linspace(0, L, n), rtol=1e-6)

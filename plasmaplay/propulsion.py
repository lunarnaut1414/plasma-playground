"""Electric-propulsion scaling and a 1-D Lorentz-accelerator channel model.

Background physics for experiment 07 (an MHD drive for space propulsion). The
thrust mechanism is the J×B body force from experiment 06 applied to make
exhaust: cross a current J with a magnetic field B and the plasma is pushed.

Two views:
  * `maecker_thrust` — the self-field MPD scaling, where the current's *own*
    azimuthal field provides B, giving the famous T ∝ I² law.
  * `channel_velocity` — a steady 1-D channel where a prescribed Lorentz force
    density accelerates the flow; integrate the cold momentum equation.
"""

from __future__ import annotations

import numpy as np

from .constants import MU_0

G0 = 9.80665   # standard gravity [m/s²], for specific impulse


def maecker_thrust(current, r_anode, r_cathode, mu0=MU_0):
    """Self-field magnetoplasmadynamic (MPD) thrust — Maecker's formula:

        T = (μ0 / 4π) I² [ ln(r_a / r_c) + 3/4 ].

    The J×B force of the discharge current interacting with its own azimuthal
    magnetic field. The hallmark is T ∝ I² (independent of mass flow), which is
    why MPD thrusters want very high currents (and megawatts of power).
    """
    return (mu0 / (4.0 * np.pi)) * current**2 * (np.log(r_anode / r_cathode) + 0.75)


def specific_impulse(thrust, mass_flow, g0=G0):
    """Specific impulse Isp = exhaust velocity / g0 = T / (ṁ g0)  [s]."""
    return thrust / (mass_flow * g0)


def jet_power(thrust, mass_flow):
    """Ideal jet (kinetic) power = T² / (2 ṁ) = ½ ṁ u_e²  [W]."""
    return thrust**2 / (2.0 * mass_flow)


def channel_velocity(force_density, dx, mass_flux_density, u_inlet=0.0):
    """Velocity along a steady 1-D accelerator channel.

    Cold single-fluid momentum balance at constant mass-flux density
    G = ρu = ṁ/A:   G du/dx = f(x),  with f the axial Lorentz force per volume.
    Integrating gives u(x) = u_inlet + (1/G) ∫₀ˣ f dx'.

    Parameters
    ----------
    force_density : (n,) axial Lorentz force per unit volume f(x) [N/m³].
    dx : grid spacing [m].
    mass_flux_density : G = ṁ/A [kg/m²/s].
    u_inlet : inflow speed [m/s].

    Returns
    -------
    u : (n,) velocity profile [m/s].
    """
    f = np.asarray(force_density, dtype=float)
    integral = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * dx)])
    return u_inlet + integral / mass_flux_density

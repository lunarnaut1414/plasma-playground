"""Particle pushers — integrators that advance a charged particle in E & B.

The Boris pusher is the workhorse of plasma PIC codes: it is second-order
accurate, time-reversible, and (crucially) energy-conserving over long runs
because it exactly rotates the velocity about B without changing its magnitude.
A naive RK4 integrator, by contrast, slowly gains or loses energy and the
gyro-orbit spirals open or closed — try it and see.
"""

from __future__ import annotations

import numpy as np


def boris_push(position, velocity, charge, mass, E_func, B_func, dt, n_steps):
    """Integrate one particle with the Boris algorithm.

    Parameters
    ----------
    position, velocity : array-like, shape (3,)
        Initial state in SI units (m, m/s).
    charge, mass : float
        Particle charge (C) and mass (kg).
    E_func, B_func : callable
        Field models, position -> (3,) vector. See plasmaplay.fields.
    dt : float
        Time step (s). Resolve the gyro-period: dt ~ T_c / 20 or smaller.
    n_steps : int
        Number of steps to take.

    Returns
    -------
    t : (n_steps + 1,) ndarray
    positions : (n_steps + 1, 3) ndarray
    velocities : (n_steps + 1, 3) ndarray
    """
    x = np.array(position, dtype=float)
    v = np.array(velocity, dtype=float)

    positions = np.empty((n_steps + 1, 3))
    velocities = np.empty((n_steps + 1, 3))
    positions[0] = x
    velocities[0] = v

    qmdt2 = (charge / mass) * 0.5 * dt  # (q/m)(dt/2), reused each step

    for i in range(n_steps):
        E = np.asarray(E_func(x), dtype=float)
        B = np.asarray(B_func(x), dtype=float)

        # Half electric acceleration
        v_minus = v + qmdt2 * E

        # Magnetic rotation (the exact rotation that conserves |v|)
        t_vec = qmdt2 * B
        t2 = t_vec @ t_vec
        s_vec = 2.0 * t_vec / (1.0 + t2)
        v_prime = v_minus + np.cross(v_minus, t_vec)
        v_plus = v_minus + np.cross(v_prime, s_vec)

        # Second half electric acceleration
        v = v_plus + qmdt2 * E

        # Drift position with the updated velocity (leapfrog)
        x = x + v * dt

        positions[i + 1] = x
        velocities[i + 1] = v

    t = np.arange(n_steps + 1) * dt
    return t, positions, velocities

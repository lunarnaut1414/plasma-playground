"""General-purpose ODE integrators.

The Boris pusher (`pushers.py`) is special-cased for the Lorentz force because
energy conservation matters there. For everything else — field-line tracing,
guiding-center drift equations — a classic 4th-order Runge-Kutta is the right
default: simple, accurate, and (unlike Boris) it works for an arbitrary
right-hand side f(t, y).
"""

from __future__ import annotations

import numpy as np


def rk4_step(f, t, y, dt):
    """One classic 4th-order Runge-Kutta step of dy/dt = f(t, y)."""
    k1 = np.asarray(f(t, y), dtype=float)
    k2 = np.asarray(f(t + 0.5 * dt, y + 0.5 * dt * k1), dtype=float)
    k3 = np.asarray(f(t + 0.5 * dt, y + 0.5 * dt * k2), dtype=float)
    k4 = np.asarray(f(t + dt, y + dt * k3), dtype=float)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def integrate_ode(f, y0, t0, dt, n_steps):
    """Integrate dy/dt = f(t, y) with fixed-step RK4.

    Parameters
    ----------
    f : callable
        Right-hand side f(t, y) -> dy/dt, where y is a (dim,) array.
    y0 : array-like
        Initial state.
    t0, dt : float
        Initial time and (fixed) step size.
    n_steps : int
        Number of steps to take.

    Returns
    -------
    t : (n_steps + 1,) ndarray
    y : (n_steps + 1, dim) ndarray
    """
    y = np.array(y0, dtype=float)
    dim = y.size
    ys = np.empty((n_steps + 1, dim))
    ys[0] = y
    t = t0
    for i in range(n_steps):
        y = rk4_step(f, t, y, dt)
        t += dt
        ys[i + 1] = y
    ts = t0 + np.arange(n_steps + 1) * dt
    return ts, ys

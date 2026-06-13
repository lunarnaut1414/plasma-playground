"""V3 — RK4 integrator accuracy.

A 4th-order method's global error should scale like dt^4: halving the step cuts
the error by ~16x. We verify both the convergence order (slope of log-error vs
log-dt) and exactness on a trivially integrable case.
"""

import numpy as np

from plasmaplay.integrators import integrate_ode


def test_v3_fourth_order_convergence():
    # Harmonic oscillator y'' = -w^2 y  ->  state [y, y'].
    w = 2 * np.pi                      # period 1
    T = 1.0

    def f(t, s):
        return np.array([s[1], -w**2 * s[0]])

    y0 = np.array([1.0, 0.0])          # exact: y(t) = cos(w t)

    dts, errs = [], []
    for n in (100, 200, 400, 800):
        dt = T / n
        t, y = integrate_ode(f, y0, 0.0, dt, n)
        err = np.max(np.abs(y[:, 0] - np.cos(w * t)))
        dts.append(dt)
        errs.append(err)

    # slope of log(err) vs log(dt) should be ~4
    slope = np.polyfit(np.log(dts), np.log(errs), 1)[0]
    assert 3.8 < slope < 4.2


def test_v3_constant_derivative_is_exact():
    # y' = c  ->  y(t) = y0 + c t, which RK4 integrates exactly.
    def f(t, y):
        return np.array([3.0, -1.0])

    t, y = integrate_ode(f, [0.0, 5.0], 0.0, 0.1, 50)
    np.testing.assert_allclose(y[:, 0], 3.0 * t, atol=1e-12)
    np.testing.assert_allclose(y[:, 1], 5.0 - t, atol=1e-12)

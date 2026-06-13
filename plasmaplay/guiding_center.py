"""Guiding-center motion.

Instead of resolving the fast gyration (experiment 01's Boris orbit), follow only
the slowly-drifting *guiding center*. The perpendicular velocity is the sum of
the standard drifts, and the parallel motion feels the mirror force:

    dX/dt = v∥ b̂ + v_E + v_∇B + v_curv
    m dv∥/dt = q E∥ − μ (b̂·∇)|B|

with the drifts
    v_E     = (E × B) / B²                      (E×B; charge-independent)
    v_∇B    = (μ/q) (b̂ × ∇|B|) / |B|            (grad-B)
    v_curv  = (m v∥²/q) (b̂ × κ) / |B|,  κ=(b̂·∇)b̂ (curvature)

and the magnetic moment μ = m v⊥²/(2|B|) conserved as an adiabatic invariant.
Field gradients are taken numerically, so any field from `plasmaplay.fields`
(or a real equilibrium) drops straight in.
"""

from __future__ import annotations

import numpy as np

from .integrators import integrate_ode


def _field_geometry(B_func, x, h=1e-6):
    """Return |B|, b̂, ∇|B|, and curvature κ=(b̂·∇)b̂ at point x (finite diff)."""
    x = np.asarray(x, dtype=float)
    B = np.asarray(B_func(x), dtype=float)
    Bmag = np.linalg.norm(B)
    bhat = B / Bmag

    gradB = np.empty(3)        # ∇|B|
    dbhat = np.empty((3, 3))   # dbhat[i, j] = ∂ b̂_i / ∂ x_j
    for j in range(3):
        xp = x.copy(); xp[j] += h
        xm = x.copy(); xm[j] -= h
        Bp = np.asarray(B_func(xp), dtype=float)
        Bm = np.asarray(B_func(xm), dtype=float)
        Bp_mag, Bm_mag = np.linalg.norm(Bp), np.linalg.norm(Bm)
        gradB[j] = (Bp_mag - Bm_mag) / (2 * h)
        dbhat[:, j] = (Bp / Bp_mag - Bm / Bm_mag) / (2 * h)

    kappa = dbhat @ bhat       # (b̂·∇)b̂
    return Bmag, bhat, gradB, kappa


def drift_velocity(x, v_par, mu, charge, mass, E_func, B_func):
    """Total guiding-center velocity dX/dt at (x, v_par) for invariant mu."""
    E = np.asarray(E_func(x), dtype=float)
    Bmag, bhat, gradB, kappa = _field_geometry(B_func, x)
    v_E = np.cross(E, bhat) / Bmag                      # (E×B)/B² = (E×b̂)/B
    v_gradB = (mu / charge) * np.cross(bhat, gradB) / Bmag
    v_curv = (mass * v_par**2 / charge) * np.cross(bhat, kappa) / Bmag
    return v_par * bhat + v_E + v_gradB + v_curv


def gc_push(position, v_par, mu, charge, mass, E_func, B_func, dt, n_steps):
    """Integrate the guiding-center equations with RK4.

    Returns t (n_steps+1,), positions (n_steps+1, 3), v_par (n_steps+1,).
    `mu` is the magnetic moment (constant); set it from the initial v_perp via
    mu = mass * v_perp**2 / (2 * |B(position)|).
    """
    def rhs(t, y):
        x, u = y[:3], y[3]
        E = np.asarray(E_func(x), dtype=float)
        Bmag, bhat, gradB, kappa = _field_geometry(B_func, x)
        E_par = E @ bhat
        v_E = np.cross(E, bhat) / Bmag
        v_gradB = (mu / charge) * np.cross(bhat, gradB) / Bmag
        v_curv = (mass * u**2 / charge) * np.cross(bhat, kappa) / Bmag
        dX = u * bhat + v_E + v_gradB + v_curv
        du = (charge / mass) * E_par - (mu / mass) * (gradB @ bhat)
        return np.array([dX[0], dX[1], dX[2], du])

    y0 = np.array([position[0], position[1], position[2], v_par])
    t, y = integrate_ode(rhs, y0, 0.0, dt, n_steps)
    return t, y[:, :3], y[:, 3]


def magnetic_moment(v_perp, mass, Bmag):
    """μ = m v⊥² / (2 |B|)."""
    return mass * v_perp**2 / (2.0 * Bmag)

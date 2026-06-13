"""Diagnostics shared across experiments.

Currently: magnetic field-line tracing and the Poincaré-section / rotational-
transform tools that turn a 3-D field into the classic 2-D "puncture plot" used
to diagnose flux surfaces and islands (experiments 04, 05).
"""

from __future__ import annotations

import numpy as np

from .integrators import rk4_step


def _unit_field(B_func):
    """Wrap B(x) as an arc-length RHS dx/ds = B/|B| for the ODE integrator."""
    def rhs(s, x):
        B = np.asarray(B_func(x), dtype=float)
        norm = np.linalg.norm(B)
        return B / norm if norm > 0 else B
    return rhs


def trace_field_line(B_func, x0, ds, n_steps):
    """Trace a magnetic field line by integrating dx/ds = B/|B| (arc length).

    Returns positions of shape (n_steps + 1, 3). Arc-length parameterization
    keeps the step size uniform along the line regardless of |B|.
    """
    rhs = _unit_field(B_func)
    x = np.array(x0, dtype=float)
    pts = np.empty((n_steps + 1, 3))
    pts[0] = x
    s = 0.0
    for i in range(n_steps):
        x = rk4_step(rhs, s, x, ds)
        s += ds
        pts[i + 1] = x
    return pts


def poincare_section(B_func, x0, period, n_crossings, ds=0.01, axis=2,
                     max_steps=2_000_000):
    """Record where a field line punctures planes spaced `period` apart along `axis`.

    Starting from x0 (assumed to sit on a plane, i.e. x0[axis] ≈ 0), trace the
    line and linearly interpolate the crossing point each time the `axis`
    coordinate passes the next multiple of `period`. Returns the crossing
    points, shape (n_crossings, 3) — the dots of a Poincaré plot.
    """
    rhs = _unit_field(B_func)
    x = np.array(x0, dtype=float)
    crossings = []
    target = x[axis] + period
    s = 0.0
    steps = 0
    while len(crossings) < n_crossings and steps < max_steps:
        x_new = rk4_step(rhs, s, x, ds)
        if x[axis] < target <= x_new[axis]:
            frac = (target - x[axis]) / (x_new[axis] - x[axis])
            crossings.append(x + frac * (x_new - x))
            target += period
        x = x_new
        s += ds
        steps += 1
    return np.array(crossings)


def rotational_transform(crossings, axis=2):
    """Rotational transform ι from a sequence of Poincaré crossings.

    ι is the average poloidal angle advanced per toroidal period, divided by 2π.
    Angles are measured in the plane perpendicular to `axis` about its centroid.
    """
    perp = [i for i in range(3) if i != axis]
    pts = crossings[:, perp]
    center = pts.mean(axis=0)
    rel = pts - center
    theta = np.unwrap(np.arctan2(rel[:, 1], rel[:, 0]))
    return np.mean(np.diff(theta)) / (2.0 * np.pi)

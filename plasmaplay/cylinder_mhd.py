"""Linear MHD stability of a periodic cylinder — the "straight tokamak" (rung B1).

The slab tearing mode (`tearing.py`, T4) is the cleanest possible reconnection
problem, but it has no *safety factor* and no *rational surfaces* — the two things
that organise tokamak stability. This module lifts the same physics to a periodic
cylinder (a torus cut and straightened, radius a, length 2 pi R), where a current
profile J_z(r) gives a poloidal field B_theta(r), a strong axial guide field B_z,
and the safety factor

    q(r) = r B_z / (R B_theta(r)).

A perturbation ~ exp(i(m theta + k z)),  k = -n/R, resonates on the **rational
surface** q(r_s) = m/n. Two instabilities live here:

  * the **m/n tearing mode** — reconnection at r_s, governed by the outer
    (marginal-ideal) Newcomb equation and its stability index Delta';
  * the **m=1 internal kink** — the rigid sideways shift of the core inside the
    q=1 surface, which exists exactly when q(0) < 1. This is the **sawtooth
    trigger**, the gateway to the Track-B/C coupling.

Outer (constant-psi) tearing equation for the perturbed helical flux psi(r):

    psi'' + psi'/r - (m^2/r^2) psi - [ mu0 J_z'(r) / (B_theta(r) (1 - n q/m)) ] psi = 0,

singular at r_s. Delta' = [psi'/psi]_{r_s^+} - [psi'/psi]_{r_s^-} is the jump in the
logarithmic derivative across the rational surface: Delta' > 0 is tearing-unstable.
The resistive growth rate then follows the same Furth-Killeen-Rosenbluth layer law
validated in the slab, gamma tau_A ~ Delta'^{4/5} S^{-3/5}.

Normalisation: a = 1, B_z = 1, R = 1 (the instability depends only on the q-profile
and the mode numbers — the overall B_theta scale R cancels in the singular term).
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp


# ---------------------------------------------------------------------------
# Screw-pinch equilibrium: a peaked current profile and its safety factor
# ---------------------------------------------------------------------------
def screw_pinch_q(r, q0, nu=1.0):
    """Safety factor q(r) of the J_z(r) = J0 (1 - (r/a)^2)^nu screw pinch (a = 1).

    Integrating Ampere's law gives the closed form

        q(r) = q0 * (nu+1) r^2 / (1 - (1 - r^2)^(nu+1)),

    with q(0) = q0 on axis and q(1) = (nu+1) q0 at the edge — a monotonically
    rising q, peaked current for larger nu. Vectorised; the r -> 0 limit is q0.
    """
    r = np.asarray(r, dtype=float)
    out = np.full_like(r, float(q0))
    nz = r > 1e-9
    rr = r[nz]
    out[nz] = q0 * (nu + 1.0) * rr ** 2 / (1.0 - (1.0 - rr ** 2) ** (nu + 1.0))
    return out


def b_theta(r, q0, nu=1.0):
    """Poloidal field B_theta(r) = r B_z / (R q(r)) with B_z = R = 1."""
    r = np.asarray(r, dtype=float)
    q = screw_pinch_q(r, q0, nu)
    return np.where(r > 1e-9, r / np.maximum(q, 1e-12), 0.0)


def rational_surface(m, n, q0, nu=1.0, a=1.0):
    """Radius r_s where q(r_s) = m/n, or None if no such surface is in (0, a].

    q rises monotonically from q0 to (nu+1) q0, so the surface exists iff
    q0 <= m/n <= (nu+1) q0. Found by bisection on the monotone q.
    """
    qs = m / n
    if not (screw_pinch_q(0.0, q0, nu) <= qs <= screw_pinch_q(a, q0, nu)):
        return None
    lo, hi = 1e-6, a
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if screw_pinch_q(mid, q0, nu) < qs:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


# ---------------------------------------------------------------------------
# Outer Newcomb equation and the tearing stability index Delta'
# ---------------------------------------------------------------------------
def _newcomb_rhs_factory(m, n, q0, nu):
    """Return f(r, y) for y = [psi, psi'] of the cylindrical Newcomb equation."""
    h = 1e-5

    def mu0_Jz(r):
        # mu0 J_z = (1/r) d(r B_theta)/dr, by a centred difference
        rp, rm = r + h, max(r - h, 1e-9)
        return (rp * b_theta(rp, q0, nu) - rm * b_theta(rm, q0, nu)) / (2 * h) / max(r, 1e-9)

    def rhs(r, y):
        psi, dpsi = y
        bt = float(b_theta(r, q0, nu))
        q = float(screw_pinch_q(r, q0, nu))
        denom = bt * (1.0 - n * q / m)
        jzp = (mu0_Jz(r + h) - mu0_Jz(r - h)) / (2 * h)   # d(mu0 J_z)/dr
        g = m ** 2 / r ** 2 + jzp / denom
        return [dpsi, g * psi - dpsi / r]

    return rhs


def delta_prime_cylinder(m, n, q0, nu=1.0, a=1.0, gap=2e-3, n_eval=600):
    """Tearing stability index Delta' for the m/n mode of the screw pinch.

    Integrates the outer Newcomb equation from the axis outward (regular solution
    psi ~ r^m) to just inside the rational surface, and from the conducting wall at
    r = a inward (psi(a) = 0) to just outside it, then returns

        Delta' = (psi'/psi)|_{r_s^+} - (psi'/psi)|_{r_s^-}.

    `gap` is the small offset from r_s on each side (the constant-psi layer is
    skipped). Returns np.nan if the m/n surface is absent. Delta' > 0 => unstable.
    """
    r_s = rational_surface(m, n, q0, nu, a)
    if r_s is None or r_s <= gap or r_s >= a - gap:
        return np.nan
    rhs = _newcomb_rhs_factory(m, n, q0, nu)

    r0 = 1e-4
    inner = solve_ivp(rhs, (r0, r_s - gap), [r0 ** m, m * r0 ** (m - 1)],
                      t_eval=np.linspace(r0, r_s - gap, n_eval), rtol=1e-8, atol=1e-10)
    li = inner.y[1, -1] / inner.y[0, -1]            # (psi'/psi)_-

    outer = solve_ivp(rhs, (a, r_s + gap), [1e-6 * a, -1.0],
                      t_eval=np.linspace(a, r_s + gap, n_eval), rtol=1e-8, atol=1e-10)
    lo = outer.y[1, -1] / outer.y[0, -1]            # (psi'/psi)_+
    return lo - li


def newcomb_eigenfunction(m, n, q0, nu=1.0, a=1.0, gap=2e-3, n_eval=400):
    """The (constant-psi) outer eigenfunction psi(r) of the m/n mode, normalised to
    psi(r_s) = 1 and continuous across the rational surface. Returns (r, psi).

    Used to *show* the mode: the inner branch (axis -> r_s) and the outer branch
    (r_s -> wall) are each normalised to 1 at r_s, so the kink/tearing structure —
    a core displacement that peaks at the resonant surface and falls to the wall —
    is visible. (The reconnected layer itself is sub-grid.)
    """
    r_s = rational_surface(m, n, q0, nu, a)
    if r_s is None:
        # no resonant surface (e.g. m=1 with q0>1): a smooth r^m core mode to the wall
        r = np.linspace(1e-4, a, 2 * n_eval)
        psi = (r / a) ** m * (1.0 - (r / a) ** 2)
        return r, psi / psi.max()
    rhs = _newcomb_rhs_factory(m, n, q0, nu)
    r0 = 1e-4
    inner = solve_ivp(rhs, (r0, r_s - gap), [r0 ** m, m * r0 ** (m - 1)],
                      t_eval=np.linspace(r0, r_s - gap, n_eval), rtol=1e-8, atol=1e-10)
    outer = solve_ivp(rhs, (a, r_s + gap), [1e-6 * a, -1.0],
                      t_eval=np.linspace(a, r_s + gap, n_eval), rtol=1e-8, atol=1e-10)
    ri, pi = inner.t, inner.y[0] / inner.y[0, -1]
    ro, po = outer.t[::-1], outer.y[0][::-1] / outer.y[0][-1]
    r = np.concatenate([ri, [r_s], ro])
    psi = np.concatenate([pi, [1.0], po])
    return r, psi


# ---------------------------------------------------------------------------
# The m = 1 internal kink (the sawtooth trigger) and the FKR growth rate
# ---------------------------------------------------------------------------
def internal_kink_unstable(q0, nu=1.0, a=1.0):
    """True if the m=1/n=1 internal kink is unstable: a q=1 surface exists, i.e.
    q(0) < 1 < q(a). The cylindrical sawtooth trigger."""
    return rational_surface(1, 1, q0, nu, a) is not None


def internal_kink_xi(r, r1, sharpness=30.0):
    """Radial displacement profile xi_r(r) of the ideal m=1 internal kink.

    In the cylinder the ideal m=1 eigenfunction is a *rigid* core shift: xi_r is
    nearly constant inside the q=1 surface r1 and falls sharply to zero outside it
    (a top hat, discontinuous in the ideal limit; here smoothed by `sharpness`).
    Normalised to xi_r(0) = 1. This is the shape that, displacing the core sideways,
    makes the characteristic internal-kink crescent. r1 is `rational_surface(1,1,...)`.
    """
    r = np.asarray(r, dtype=float)
    return 0.5 * (1.0 - np.tanh(sharpness * (r - r1)))


def fkr_growth_rate(delta_prime, S, k_factor=0.55):
    """Furth-Killeen-Rosenbluth resistive-tearing growth rate gamma*tau_A.

        gamma tau_A = k_factor * Delta'^{4/5} S^{-3/5}    (constant-psi regime)

    The resistive-layer law is local, so it carries over unchanged from the slab to
    the cylinder — only Delta' becomes the cylindrical value. Returns 0 for a stable
    (Delta' <= 0) surface. The S^{-3/5} scaling is the falsifiable signature.
    """
    dp = np.asarray(delta_prime, dtype=float)
    return np.where(dp > 0, k_factor * np.abs(dp) ** 0.8 * S ** (-0.6), 0.0)

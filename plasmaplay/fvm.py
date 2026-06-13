"""Finite-volume ideal MHD, 1-D.

Treat the plasma as a single magnetized fluid and solve the ideal MHD
conservation laws in one dimension with a shock-capturing Godunov scheme:

    MUSCL (minmod) reconstruction  ->  HLL Riemann flux  ->  SSP-RK2 in time

State (primitive) is W = [ρ, u, v, w, p, By, Bz]; Bx is constant in 1-D because
∂Bx/∂x = 0 follows from ∇·B = 0. Conserved variables are
U = [ρ, ρu, ρv, ρw, By, Bz, E] with E = p/(γ-1) + ½ρ|v|² + ½|B|².

The HLL solver is robust and positivity-friendly but diffusive (it smears
contacts and the slow compound wave); HLLD would sharpen those. That trade is
fine for learning — the shock/rarefaction structure comes out clearly.

Indices into the last axis:
    W: 0=ρ 1=u 2=v 3=w 4=p 5=By 6=Bz
    U: 0=ρ 1=ρu 2=ρv 3=ρw 4=By 5=Bz 6=E
"""

from __future__ import annotations

import numpy as np


def prim_to_cons(W, Bx, gamma):
    rho, u, v, w, p, By, Bz = W.T
    B2 = Bx**2 + By**2 + Bz**2
    E = p / (gamma - 1) + 0.5 * rho * (u**2 + v**2 + w**2) + 0.5 * B2
    return np.stack([rho, rho * u, rho * v, rho * w, By, Bz, E], axis=-1)


def cons_to_prim(U, Bx, gamma):
    rho = U[..., 0]
    u = U[..., 1] / rho
    v = U[..., 2] / rho
    w = U[..., 3] / rho
    By = U[..., 4]
    Bz = U[..., 5]
    E = U[..., 6]
    B2 = Bx**2 + By**2 + Bz**2
    p = (gamma - 1) * (E - 0.5 * rho * (u**2 + v**2 + w**2) - 0.5 * B2)
    return np.stack([rho, u, v, w, p, By, Bz], axis=-1)


def flux(W, Bx, gamma):
    """Physical MHD flux in the x-direction from primitive state W."""
    rho, u, v, w, p, By, Bz = W.T
    B2 = Bx**2 + By**2 + Bz**2
    ptot = p + 0.5 * B2
    vdotB = u * Bx + v * By + w * Bz
    E = p / (gamma - 1) + 0.5 * rho * (u**2 + v**2 + w**2) + 0.5 * B2
    return np.stack([
        rho * u,
        rho * u**2 + ptot - Bx**2,
        rho * u * v - Bx * By,
        rho * u * w - Bx * Bz,
        u * By - v * Bx,
        u * Bz - w * Bx,
        (E + ptot) * u - Bx * vdotB,
    ], axis=-1)


def fast_speed(W, Bx, gamma):
    """Fast magnetosonic speed in the x-direction."""
    rho, u, v, w, p, By, Bz = W.T
    a2 = gamma * p / rho
    b2 = (Bx**2 + By**2 + Bz**2) / rho
    bx2 = Bx**2 / rho
    disc = np.maximum((a2 + b2) ** 2 - 4.0 * a2 * bx2, 0.0)
    return np.sqrt(0.5 * (a2 + b2 + np.sqrt(disc)))


def hll_flux(WL, WR, Bx, gamma):
    """HLL approximate Riemann flux between left/right interface states."""
    FL = flux(WL, Bx, gamma)
    FR = flux(WR, Bx, gamma)
    UL = prim_to_cons(WL, Bx, gamma)
    UR = prim_to_cons(WR, Bx, gamma)
    uL, uR = WL[..., 1], WR[..., 1]
    cfL = fast_speed(WL, Bx, gamma)
    cfR = fast_speed(WR, Bx, gamma)
    SL = np.minimum(uL, uR) - np.maximum(cfL, cfR)
    SR = np.maximum(uL, uR) + np.maximum(cfL, cfR)
    SLc, SRc = SL[:, None], SR[:, None]
    F = (SRc * FL - SLc * FR + (SLc * SRc) * (UR - UL)) / (SRc - SLc)
    F = np.where(SLc >= 0, FL, F)
    F = np.where(SRc <= 0, FR, F)
    return F


def _minmod(a, b):
    return np.where(a * b > 0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)


def _pad(W, bc):
    mode = "wrap" if bc == "periodic" else "edge"   # 'edge' = transmissive
    return np.pad(W, ((2, 2), (0, 0)), mode=mode)


def _rhs(U, dx, Bx, gamma, bc):
    """Conserved-variable time derivative dU/dt for the interior cells."""
    W = cons_to_prim(U, Bx, gamma)
    Wg = _pad(W, bc)                                   # (n+4, 7)
    slope = np.zeros_like(Wg)
    slope[1:-1] = _minmod(Wg[1:-1] - Wg[:-2], Wg[2:] - Wg[1:-1])
    WL = Wg[:-1] + 0.5 * slope[:-1]                    # left state at each interface
    WR = Wg[1:] - 0.5 * slope[1:]                      # right state at each interface
    Fint = hll_flux(WL, WR, Bx, gamma)                 # (n+3, 7)
    dU = -(Fint[1:] - Fint[:-1]) / dx                  # (n+2, 7), cells 1..n+2
    return dU[1:-1]                                    # interior physical cells


def solve_mhd_1d(W0, dx, Bx, gamma=2.0, t_end=0.1, cfl=0.4, bc="transmissive"):
    """Evolve the 1-D ideal MHD equations to time t_end.

    Parameters
    ----------
    W0 : (n, 7) primitive initial state.
    dx : cell size. Bx : constant normal field. gamma : adiabatic index.
    t_end : final time. cfl : Courant number. bc : 'transmissive' or 'periodic'.

    Returns
    -------
    W : (n, 7) primitive state at t_end.
    """
    U = prim_to_cons(np.asarray(W0, dtype=float), Bx, gamma)
    t = 0.0
    while t < t_end - 1e-14:
        W = cons_to_prim(U, Bx, gamma)
        smax = np.max(np.abs(W[:, 1]) + fast_speed(W, Bx, gamma))
        dt = min(cfl * dx / smax, t_end - t)
        # SSP-RK2 (Heun)
        U1 = U + dt * _rhs(U, dx, Bx, gamma, bc)
        U = 0.5 * U + 0.5 * (U1 + dt * _rhs(U1, dx, Bx, gamma, bc))
        t += dt
    return cons_to_prim(U, Bx, gamma)

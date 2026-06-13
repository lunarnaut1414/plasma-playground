"""The 3-D tokamak field — bridging the experiment-04 equilibrium into 3-D space.

A tokamak equilibrium is *axisymmetric*: experiment 04 solved Grad-Shafranov for
the poloidal flux ψ(R, Z) on a 2-D (R, Z) grid. The **full magnetic field** of
that equilibrium follows directly, with no new physics:

    B = ∇φ × ∇ψ + F(ψ) ∇φ            (∇φ = φ̂ / R)

      B_R = -(1/R) ∂ψ/∂Z
      B_Z =  (1/R) ∂ψ/∂R
      B_φ =  F(ψ) / R                 (F = R B_φ, the toroidal-field function)

`equilibrium_field` differentiates a ψ(R, Z) grid, adds the toroidal field, and
wraps the result as the usual callable ``B(position) -> (3,) ndarray`` in
**Cartesian** (x, y, z) — so it drops straight into the pushers (`boris_push`,
`gc_push`) and the field-line tracer (`trace_field_line`, `poincare_section`)
you already built and validated. That one bridge is the whole of rung T0 in
`docs/3D_TOKAMAK_GUIDE.md`; T1 (q-profile/Poincaré) and T2 (banana orbits) are
then almost pure reuse.

Two facts make this field trustworthy, and both are tested:
  * **∇·B = 0 analytically.** With B_R = -(1/R)∂ψ/∂Z and B_Z = (1/R)∂ψ/∂R the
    poloidal divergence (1/R)∂(R B_R)/∂R + ∂B_Z/∂Z cancels term-by-term (mixed
    partials of ψ), and ∂B_φ/∂φ = 0 by axisymmetry. The numeric check confirms it.
  * **|B| ∝ 1/R on axis.** A vacuum toroidal field is F/R with F constant, the
    defining 1/R falloff of a tokamak; on the magnetic axis (a ψ extremum, where
    ∇ψ = 0) the field is purely toroidal.
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import RegularGridInterpolator


def vacuum_F(R0: float, B0: float):
    """Constant toroidal-field function F(ψ) = R₀ B₀ (the vacuum / current-free case).

    Gives B_φ = R₀ B₀ / R — a pure 1/R toroidal field, the field an external
    toroidal-field coil set would produce with no plasma current contributing to
    F. The simplest valid choice, and the one whose |B| ∝ 1/R is easiest to check.
    """
    val = R0 * B0

    def F(psi):
        return np.full_like(np.asarray(psi, dtype=float), val)

    return F


def solovev_F(F0: float, FFprime: float, psi_b: float = 0.0):
    """Solov'ev toroidal-field function with F F' = const.

        F(ψ) = sign(F0) · √( F0² + 2 F F' (ψ − ψ_b) )

    The diamagnetic/paramagnetic correction the plasma current adds to the vacuum
    F0 = R₀ B₀: F F' = const integrates to this square-root profile. ψ_b is the
    flux at the reference surface where F = F0 (usually the boundary). The
    argument is clipped at 0 so the root stays real on grids that overshoot.
    """
    sign = 1.0 if F0 >= 0 else -1.0

    def F(psi):
        arg = F0**2 + 2.0 * FFprime * (np.asarray(psi, dtype=float) - psi_b)
        return sign * np.sqrt(np.clip(arg, 0.0, None))

    return F


def equilibrium_field(R, Z, psi, F_of_psi):
    """Build the 3-D tokamak B-field callable from a poloidal-flux grid ψ(R, Z).

    Parameters
    ----------
    R, Z : 1-D arrays of grid coordinates (uniform spacing; R > 0). The same
        arrays passed to `solvers.grad_shafranov_solve`.
    psi : (nR, nZ) array — the poloidal flux, ``psi[i, j] = ψ(R[i], Z[j])``
        (i.e. ``meshgrid(R, Z, indexing="ij")`` layout, as experiment 04 returns).
    F_of_psi : callable ``F(ψ) -> F`` mapping flux to the toroidal-field function
        F = R B_φ. Use `vacuum_F(R0, B0)` for a pure 1/R field or `solovev_F(...)`
        for the FF' = const profile. Receives and returns arrays.

    Returns
    -------
    field : callable ``field(position) -> (3,) ndarray``
        Cartesian B = (B_x, B_y, B_z) at ``position = (x, y, z)``. Drops into any
        pusher or tracer unchanged. B_R, B_Z come from a centered finite
        difference of ψ on the grid; B_φ = F(ψ)/R; all three are bilinearly
        interpolated onto the requested (R, Z), then rotated into Cartesian.
    """
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    psi = np.asarray(psi, dtype=float)
    if psi.shape != (R.size, Z.size):
        raise ValueError(f"psi shape {psi.shape} != (nR, nZ) = ({R.size}, {Z.size})")

    RR = R[:, None]                                   # broadcast 1/R over Z
    dpsi_dR, dpsi_dZ = np.gradient(psi, R, Z)         # ∂ψ/∂R, ∂ψ/∂Z on the grid
    BR_grid = -dpsi_dZ / RR                           # B_R = -(1/R) ∂ψ/∂Z
    BZ_grid = dpsi_dR / RR                            # B_Z =  (1/R) ∂ψ/∂R
    Bphi_grid = F_of_psi(psi) / RR                    # B_φ =  F(ψ) / R

    # Bilinear interpolation; extrapolate linearly past the grid so a field line
    # or orbit that briefly steps outside still gets a finite field.
    kw = dict(method="linear", bounds_error=False, fill_value=None)
    interp_BR = RegularGridInterpolator((R, Z), BR_grid, **kw)
    interp_BZ = RegularGridInterpolator((R, Z), BZ_grid, **kw)
    interp_Bphi = RegularGridInterpolator((R, Z), Bphi_grid, **kw)

    def field(position):
        x, y, z = position
        r = np.hypot(x, y)
        pt = (r, z)
        BR = float(interp_BR(pt))
        BZ = float(interp_BZ(pt))
        Bphi = float(interp_Bphi(pt))
        if r > 0.0:
            cos_p, sin_p = x / r, y / r               # cos φ, sin φ
        else:
            cos_p, sin_p = 1.0, 0.0
        # rotate (B_R, B_φ) from cylindrical into Cartesian
        Bx = BR * cos_p - Bphi * sin_p
        By = BR * sin_p + Bphi * cos_p
        return np.array([Bx, By, BZ])

    return field


def to_cylindrical(position):
    """Cartesian (x, y, z) -> (R, φ, Z). φ in (-π, π]."""
    x, y, z = position
    return np.hypot(x, y), np.arctan2(y, x), z


def to_cartesian(R, phi, Z):
    """Cylindrical (R, φ, Z) -> Cartesian (x, y, z) position vector."""
    return np.array([R * np.cos(phi), R * np.sin(phi), Z])


def divergence(B_func, position, h=1e-5):
    """Numeric ∇·B at a Cartesian point via centered differences.

    The headline T0 validation: a magnetic field must be divergence-free, and the
    equilibrium field is so by construction. Returns a scalar that should be ≈ 0
    (to the finite-difference truncation error) anywhere the interpolation is
    smooth — i.e. away from the grid edges.
    """
    x = np.asarray(position, dtype=float)
    div = 0.0
    for j in range(3):
        xp = x.copy(); xp[j] += h
        xm = x.copy(); xm[j] -= h
        div += (B_func(xp)[j] - B_func(xm)[j]) / (2.0 * h)
    return div

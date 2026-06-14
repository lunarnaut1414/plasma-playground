"""Flux-surface-averaged geometry from a Grad-Shafranov equilibrium.

The 1-D transport rung (exp 09 F2) runs on a circular column: the only geometry it
knows is the cylindrical volume element V' ~ rho. A *real* tokamak plasma sits on
the nested, shaped, outboard-shifted flux surfaces of a Grad-Shafranov equilibrium
(exp 04). To run transport on that equilibrium you do not solve a 2-D transport
problem — you keep it 1-D but carry the geometry through two flux-surface-averaged
metric coefficients:

    V'(rho)        = dV/drho, the volume between neighbouring flux surfaces
    <|grad rho|^2> = the flux-surface average of |grad rho|^2

so the transport divergence becomes  (1/V') d/drho( V' <|grad rho|^2> n chi dT/drho ).
When the surfaces are circular this reduces exactly to the cylindrical operator.

This module extracts those metrics from a gridded poloidal flux psi(R,Z) using the
**volume-derivative identity** for a flux-surface average,

    <X>(rho) = d/dV integral_{rho'<rho} X dV  =  (dI_X/drho) / (dV/drho),

evaluated by binning grid cells into flux-label shells (no fragile contour tracing).
The flux label is rho = sqrt(psi_n), psi_n the normalized poloidal flux (0 on axis,
1 on the boundary). dV = 2 pi R dR dZ is the toroidal volume element.
"""

from __future__ import annotations

import numpy as np


def normalized_flux(psi, psi_axis=None, psi_bnd=0.0):
    """Normalized poloidal flux psi_n = (psi_axis - psi)/(psi_axis - psi_bnd).

    psi_n = 0 at the magnetic axis (the psi extremum) and 1 on the plasma boundary.
    If `psi_axis` is None it is taken as psi.max() (the Solov'ev case here peaks at
    the axis). Returns an array the shape of `psi`.
    """
    psi = np.asarray(psi, dtype=float)
    psi_axis = float(psi.max()) if psi_axis is None else float(psi_axis)
    return (psi_axis - psi) / (psi_axis - psi_bnd)


def flux_surface_metrics(R, Z, psi, *, n_rho=64, psi_axis=None, psi_bnd=0.0):
    """Flux-surface-averaged geometry from a gridded poloidal flux psi(R, Z).

    Parameters
    ----------
    R, Z : 1-D uniform coordinate arrays; psi has shape (R.size, Z.size) with
        ``indexing="ij"`` (psi[i, j] at (R[i], Z[j])) — the convention of
        `solvers.grad_shafranov_solve`.
    n_rho : number of flux-label shells (the returned radial resolution).
    psi_axis, psi_bnd : flux at the axis / boundary (defaults: psi.max(), 0).

    Returns
    -------
    dict with
        rho        : (n_rho,) shell-centre flux label sqrt(psi_n) in (0, 1)
        Vprime     : (n_rho,) dV/drho [m^3] — the volume metric
        grad_rho2  : (n_rho,) <|grad rho|^2> [m^-2] — the gradient metric
        V          : (n_rho,) enclosed volume [m^3]
        psi_n      : (nR, nZ) normalized flux on the grid (for rendering)
        rho_grid   : (nR, nZ) sqrt(psi_n) clipped to [0, 1] (for rendering)

    The bin average is the discrete form of <X> = (dI_X/drho)/(dV/drho): each shell
    sums X dV and dV over the cells whose rho falls in it, then divides.
    """
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    psi = np.asarray(psi, dtype=float)
    psi_n = normalized_flux(psi, psi_axis, psi_bnd)
    rho_grid = np.sqrt(np.clip(psi_n, 0.0, 1.0))

    # |grad rho|^2 = |grad psi_n|^2 / (4 rho^2), with rho = sqrt(psi_n).
    dpsidR, dpsidZ = np.gradient(psi_n, R, Z)
    grad_psn2 = dpsidR ** 2 + dpsidZ ** 2
    rho_safe = np.maximum(rho_grid, 1e-6)
    gradrho2_cell = grad_psn2 / (4.0 * rho_safe ** 2)

    # toroidal volume element of each grid cell: dV = 2 pi R dR dZ
    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]
    dV_cell = 2.0 * np.pi * R[:, None] * dR * dZ * np.ones_like(psi)

    inside = (psi_n >= 0.0) & (psi_n <= 1.0)
    edges = np.linspace(0.0, 1.0, n_rho + 1)
    drho = edges[1] - edges[0]
    centers = 0.5 * (edges[:-1] + edges[1:])

    rg = rho_grid[inside]
    dv = dV_cell[inside]
    V_shell, _ = np.histogram(rg, bins=edges, weights=dv)
    IX_shell, _ = np.histogram(rg, bins=edges, weights=gradrho2_cell[inside] * dv)

    with np.errstate(divide="ignore", invalid="ignore"):
        grad_rho2 = np.where(V_shell > 0, IX_shell / V_shell, 0.0)
    Vprime = V_shell / drho
    V = np.cumsum(V_shell)

    # fill any empty inner shells (coarse grid near the axis) by nearest-valid value
    good = V_shell > 0
    if not good.all() and good.any():
        idx = np.arange(n_rho)
        grad_rho2 = np.interp(idx, idx[good], grad_rho2[good])
        Vprime = np.interp(idx, idx[good], Vprime[good])

    return {
        "rho": centers, "Vprime": Vprime, "grad_rho2": grad_rho2, "V": V,
        "psi_n": psi_n, "rho_grid": rho_grid,
    }


def confinement_time_ipb98(Ip_MA, B, n19, P_MW, R, a, kappa, M=2.5):
    """ITER IPB98(y,2) H-mode energy-confinement scaling, tau_E [s].

        tau_E = 0.0562 Ip^0.93 B^0.15 n19^0.41 P^-0.69 R^1.97 kappa^0.78 eps^0.58 M^0.19

    with Ip [MA], B [T], n19 line-averaged density [1e19 m^-3], P [MW] loss power,
    R major radius [m], a minor radius [m], kappa elongation, eps = a/R, M effective
    ion mass [amu, ~2.5 for D-T]. The standard empirical tokamak-confinement law
    (ITER Physics Basis, Nucl. Fusion 39, 2175 (1999)); it anchors what chi a
    reduced transport model should produce. For the ITER baseline it gives ~3.7 s.
    """
    eps = a / R
    return (0.0562 * Ip_MA ** 0.93 * B ** 0.15 * n19 ** 0.41 * P_MW ** -0.69
            * R ** 1.97 * kappa ** 0.78 * eps ** 0.58 * M ** 0.19)

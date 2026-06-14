"""Operational limits and confinement modes of a tokamak discharge.

A burning plasma does not live on a single happy path — it lives inside an
*operating window* bounded by hard limits and split by confinement bifurcations.
This module adds the three that decide whether and how a discharge runs:

  * the **Greenwald density limit** n_G = Ip/(pi a^2): fuel past it and the edge
    cools, confinement collapses, and the plasma disrupts;
  * the **L->H transition**: above a heating-power threshold the edge transport
    bifurcates and confinement roughly doubles (the H-mode pedestal) — real
    machines live in H-mode;
  * the **density-limit / radiative collapse**: driving the density toward n_G
    degrades confinement, which lowers the temperature, which (at fixed heating)
    cannot be sustained — a downward thermal runaway that the burn does not survive.

These are returned as **confinement multipliers** on tau_E so they drop straight
into the 0-D burn (`transport.burn_0d_ash`'s `tau_factor` hook) — and as the two
empirical scaling numbers (n_G, P_LH) that the tests pin against published values.
SI units: densities m^-3, Ip in MA, a/R in m, B in T, powers in MW, area in m^2.
"""

from __future__ import annotations

import numpy as np


def greenwald_density(Ip_MA, a):
    """Greenwald density limit n_G [m^-3] = (Ip[MA] / (pi a[m]^2)) x 1e20.

    The empirical line-averaged density a tokamak can hold before the edge cools and
    it disrupts. For the ITER baseline (Ip=15 MA, a=2 m) it is ~1.2e20 m^-3, and
    machines run at ~0.5-0.85 n_G.
    """
    return Ip_MA / (np.pi * a ** 2) * 1.0e20


def lh_power_threshold(n20, B, S):
    """L->H power threshold P_LH [MW] (Martin 2008 scaling, ITPA).

        P_LH = 0.0488 * n20^0.717 * B^0.803 * S^0.941

    with n20 the line-averaged density [1e20 m^-3], B the toroidal field [T], and S
    the plasma surface area [m^2]. For ITER (n20~0.5, B=5.3, S~680) it gives ~50 MW;
    the heating must exceed it to access H-mode.
    """
    return 0.0488 * n20 ** 0.717 * B ** 0.803 * S ** 0.941


def confinement_factor_lh(p_heat_MW, p_lh_MW, h_factor=2.0, width=0.08):
    """Confinement multiplier across the L->H transition (a smooth bifurcation).

    Returns ~1 in L-mode (heating below threshold) and ~`h_factor` in H-mode (above):

        f = 1 + (h_factor - 1) * 0.5 * (1 + tanh((P/P_LH - 1)/width)).

    The sharp-but-finite `width` models the bifurcation as a fast-but-continuous
    switch so the ODE stays well-behaved. h_factor ~ 2 is the canonical H-mode jump.
    """
    x = (np.asarray(p_heat_MW, float) / p_lh_MW - 1.0) / width
    return 1.0 + (h_factor - 1.0) * 0.5 * (1.0 + np.tanh(x))


def confinement_factor_greenwald(n, n_G, *, onset=0.8, floor=0.08, width=0.06):
    """Confinement multiplier for the density limit: ~1 below the limit, collapsing
    toward `floor` as the density crosses `onset`*n_G toward n_G.

        f = floor + (1 - floor) * 0.5 * (1 - tanh((n/n_G - onset)/width)).

    Pushing the density toward the Greenwald value cools the edge and wrecks
    confinement; this is the knob that turns an over-fuelled burn into a disruption.
    The collapse is reversible — back the density off and f returns to 1.
    """
    x = (np.asarray(n, float) / n_G - onset) / width
    return floor + (1.0 - floor) * 0.5 * (1.0 - np.tanh(x))

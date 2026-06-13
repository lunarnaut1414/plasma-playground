"""Plasma wave dispersion relations (experiment 08).

The plasma dispersion function Z(ζ) and the kinetic Langmuir-wave dispersion
relation. Solving the *kinetic* relation gives a complex frequency: the real
part is the Bohm-Gross oscillation frequency and the imaginary part is the
collisionless **Landau damping** rate — both falling out of the same equation,
which is the content fluid theory misses.

Conventions: frequencies in units of ω_pe, wavenumbers via κ = k λ_D, thermal
speed v_th = √(2 k_B T/m) so that ζ = ω / (k v_th) = ω̃ / (√2 κ).
"""

from __future__ import annotations

import numpy as np
from scipy.special import wofz


def plasma_dispersion_function(zeta):
    """Z(ζ) = i√π · w(ζ), with w the Faddeeva function (scipy.special.wofz)."""
    return 1j * np.sqrt(np.pi) * wofz(zeta)


def _Zprime(zeta, Z):
    """Z'(ζ) = -2 (1 + ζ Z(ζ))."""
    return -2.0 * (1.0 + zeta * Z)


def bohm_gross(k_lambda_D):
    """Bohm-Gross (fluid) frequency ω/ω_pe = √(1 + 3 κ²) — the F0/warm-fluid limit."""
    return np.sqrt(1.0 + 3.0 * k_lambda_D**2)


def langmuir_dispersion(k_lambda_D, max_iter=100, tol=1e-12):
    """Solve the kinetic Langmuir dispersion relation for complex ω/ω_pe.

    Dielectric (electrons, Maxwellian): ε = 1 + (1/κ²)(1 + ζ Z(ζ)) = 0, i.e.
    1 + ζ Z(ζ) = -κ², solved for ζ by complex Newton iteration seeded at the
    Bohm-Gross root. Returns (omega_r, gamma) in units of ω_pe; gamma < 0 is
    damping.
    """
    kappa = float(k_lambda_D)
    # seed ζ from the Bohm-Gross frequency
    omega_seed = bohm_gross(kappa)
    zeta = omega_seed / (np.sqrt(2.0) * kappa) - 0.01j

    for _ in range(max_iter):
        Z = plasma_dispersion_function(zeta)
        g = 1.0 + zeta * Z + kappa**2            # = ε·κ²
        gp = Z + zeta * _Zprime(zeta, Z)         # dg/dζ
        step = g / gp
        zeta -= step
        if abs(step) < tol:
            break

    omega = np.sqrt(2.0) * kappa * zeta          # ω/ω_pe
    return omega.real, omega.imag

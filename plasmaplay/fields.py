"""Analytic electric and magnetic field models.

Each field is a callable f(position) -> 3-vector, so they compose and can be
swapped freely into the particle pushers. `position` is an (x, y, z) array.
"""

from __future__ import annotations

import numpy as np


def uniform_B(Bz: float = 1.0):
    """Uniform magnetic field along +z. The simplest confinement: particles
    gyrate in the x-y plane and stream freely along z."""
    B = np.array([0.0, 0.0, Bz])

    def field(position):
        return B

    return field


def uniform_E(Ex: float = 0.0, Ey: float = 0.0, Ez: float = 0.0):
    """Uniform electric field. Crossed with a uniform B it produces the classic
    E x B drift, v_drift = (E x B) / B^2 — independent of charge and mass."""
    E = np.array([Ex, Ey, Ez])

    def field(position):
        return E

    return field


def zero_field():
    """No field — useful as a default E when only B matters."""
    zero = np.zeros(3)

    def field(position):
        return zero

    return field


def magnetic_mirror(B0: float = 1.0, mirror_ratio: float = 2.0, length: float = 1.0):
    """Simple axisymmetric magnetic-mirror field (paraxial approximation).

    B_z grows toward both ends (|z| -> length/2), so the field strengthens at
    the throats and traps particles whose pitch angle is large enough. The
    radial component follows from div(B) = 0. This is the confinement principle
    behind magnetic-mirror machines and, qualitatively, the tokamak/stellarator
    trapping of banana orbits.
    """
    k = 2.0 * np.pi / length

    def field(position):
        x, y, z = position
        # On-axis field: B0 modulated so it is minimum at z=0, larger at ends.
        amp = (mirror_ratio - 1.0) / 2.0
        Bz = B0 * (1.0 + amp * (1.0 - np.cos(k * z)))
        dBz_dz = B0 * amp * k * np.sin(k * z)
        # Radial field from div(B)=0: B_r ≈ -(r/2) dBz/dz
        Bx = -0.5 * x * dBz_dz
        By = -0.5 * y * dBz_dz
        return np.array([Bx, By, Bz])

    return field

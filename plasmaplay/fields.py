"""Analytic electric and magnetic field models.

Each field is a callable f(position) -> 3-vector, so they compose and can be
swapped freely into the particle pushers. `position` is an (x, y, z) array.
"""

from __future__ import annotations

import numpy as np
from scipy.special import iv, ivp

from .constants import MU_0


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


def gradient_B(B0: float = 1.0, L: float = 1.0):
    """Straight field with a transverse gradient: B = B0 (1 + x/L) ẑ.

    |B| increases in +x, so a gyrating particle feels a grad-B drift
    v_∇B = (m v⊥²/2qB) (b̂×∇B)/B in ±y (sign set by charge). The simplest field
    for isolating the grad-B drift (the field line is straight, so there is no
    curvature drift).
    """
    def field(position):
        x, y, z = position
        return np.array([0.0, 0.0, B0 * (1.0 + x / L)])

    return field


def screw_pinch(Bz: float = 1.0, twist: float = 0.2, b_theta=None):
    """Idealized 'straight stellarator' / screw-pinch field.

        B = Bz ẑ + B_θ(r) θ̂          (cylindrical; θ̂ = (-y, x)/r)

    Field lines lie on nested cylinders r = const and wind helically. With the
    default linear profile B_θ = twist · r the winding is rigid: every surface
    has the same rotational transform per axial length L,

        ι = twist · L / (2π Bz),

    which makes it the clean analytic check for field-line tracing (test V13).
    Pass a callable ``b_theta(r)`` for a sheared ι(r) profile (used by the
    experiment to make non-trivial flux surfaces). The field is divergence-free
    for any B_θ(r) because B_θ depends only on r.
    """
    def field(position):
        x, y, z = position
        r = np.hypot(x, y)
        Bth = twist * r if b_theta is None else b_theta(r)
        if r == 0.0:
            return np.array([0.0, 0.0, Bz])
        # azimuthal unit vector θ̂ = (-y, x)/r
        return np.array([-Bth * y / r, Bth * x / r, Bz])

    return field


def helical_stellarator(B0: float = 1.0, eps: float = 0.4, l: int = 2, h: float = 1.0):
    """A genuine **vacuum** stellarator field — rotational transform from geometry,
    NOT from plasma current.

    The defining difference from a tokamak (and from `screw_pinch`, whose twist comes
    from an axial current): a stellarator's field lines twist purely from the 3-D
    shape of the *external* coils, with ~zero net plasma current. We model the
    straight-stellarator limit as a strong axial guide field plus a single helical
    harmonic, written as the gradient of a harmonic scalar potential (so the field is
    curl-free — current-free — in the plasma region by construction):

        Phi = B0 z + eps I_l(h r) cos(l theta - h z),    grad^2 Phi = 0,

    giving, in cylindrical (r, theta, z),

        B_r     =  eps h I_l'(h r) cos(l theta - h z)
        B_theta = -eps (l/r) I_l(h r) sin(l theta - h z)
        B_z     =  B0 + eps h I_l(h r) sin(l theta - h z)

    A single helicity averages to zero twist at first order; the rotational transform
    appears at **second order**, iota ∝ eps^2 — the hallmark current-free stellarator
    transform. `l` is the field periodicity (l=2 is the classic), `h = 2 pi / L` the
    axial wavenumber. Returns the usual callable ``field(position) -> (3,) ndarray``
    with position in Cartesian (x, y, z).
    """
    def field(position):
        x, y, z = position
        r = np.hypot(x, y)
        if r < 1e-9:
            return np.array([0.0, 0.0, B0])
        theta = np.arctan2(y, x)
        u = l * theta - h * z
        Il = iv(l, h * r)
        Ilp = ivp(l, h * r)
        Br = eps * h * Ilp * np.cos(u)
        Bth = -eps * (l / r) * Il * np.sin(u)
        Bz = B0 + eps * h * Il * np.sin(u)
        return np.array([Br * x / r - Bth * y / r, Br * y / r + Bth * x / r, Bz])

    return field


def circular_loop_coil(radius: float, n_segments: int = 200,
                       center=(0.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)):
    """Return (n_segments, 3) points sampling a circular current loop.

    The points are vertices of a closed polygon (the segment from the last
    vertex back to the first closes the loop) — feed straight into
    `biot_savart`. `normal` sets the loop's orientation.
    """
    n = np.asarray(normal, dtype=float)
    n /= np.linalg.norm(n)
    # pick any vector not parallel to n, then build an orthonormal in-plane basis
    seed = np.array([1.0, 0.0, 0.0]) if abs(n[0]) < 0.9 else np.array([0.0, 1.0, 0.0])
    u = np.cross(seed, n)
    u /= np.linalg.norm(u)
    v = np.cross(n, u)
    phi = np.linspace(0.0, 2.0 * np.pi, n_segments, endpoint=False)
    pts = (np.asarray(center, dtype=float)
           + radius * (np.outer(np.cos(phi), u) + np.outer(np.sin(phi), v)))
    return pts


def biot_savart(coil, current: float):
    """Magnetic field of a closed current loop via the Biot-Savart law.

        B(r) = (μ0 I / 4π) Σ  dl × (r - r_mid) / |r - r_mid|³

    `coil` is an (N, 3) array of vertices of a closed polygon (e.g. from
    `circular_loop_coil`); each consecutive pair is a straight current segment,
    and the loop closes from the last vertex back to the first. Returns the usual
    callable ``field(position) -> (3,) ndarray``.
    """
    coil = np.asarray(coil, dtype=float)
    starts = coil
    ends = np.roll(coil, -1, axis=0)
    dl = ends - starts                       # (N, 3) segment vectors
    mid = 0.5 * (starts + ends)              # (N, 3) segment midpoints
    const = MU_0 * current / (4.0 * np.pi)

    def field(position):
        r = np.asarray(position, dtype=float) - mid          # (N, 3)
        r_norm = np.linalg.norm(r, axis=1)                   # (N,)
        dB = const * np.cross(dl, r) / r_norm[:, None] ** 3  # (N, 3)
        return dB.sum(axis=0)

    return field


def circular_loop(radius: float, current: float, n_segments: int = 200,
                  center=(0.0, 0.0, 0.0), normal=(0.0, 0.0, 1.0)):
    """Convenience: Biot-Savart field of a single circular current loop."""
    coil = circular_loop_coil(radius, n_segments, center, normal)
    return biot_savart(coil, current)

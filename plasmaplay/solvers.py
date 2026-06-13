"""Field solvers.

Spectral (FFT) Poisson solvers on periodic grids — the electrostatic heart of a
PIC code. Given a charge density ρ, solve

    ∇²φ = -ρ/ε0,      E = -∇φ

In Fourier space the Laplacian is just multiplication by -k², so the solve is a
forward FFT, a divide by k², and an inverse FFT — exact (to round-off) for every
resolved mode. The k = 0 component is set to zero, which is equivalent to adding
a uniform neutralizing background so the net charge integrates to zero (required
for a periodic solution to exist).
"""

from __future__ import annotations

import numpy as np

from .constants import EPSILON_0


def solve_poisson_1d(rho, dx, eps0=EPSILON_0):
    """Periodic 1-D Poisson solve: return potential φ with ∇²φ = -ρ/ε0.

    `rho` is charge density on a uniform periodic grid of spacing `dx`. φ has zero
    spatial mean (the k=0 / neutralizing-background gauge).
    """
    rho = np.asarray(rho, dtype=float)
    n = rho.size
    k = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    rho_k = np.fft.fft(rho)
    phi_k = np.zeros_like(rho_k)
    nz = k != 0.0
    phi_k[nz] = rho_k[nz] / (eps0 * k[nz] ** 2)
    return np.fft.ifft(phi_k).real


def solve_efield_1d(rho, dx, eps0=EPSILON_0):
    """Periodic 1-D electric field E = -dφ/dx directly from ρ (spectral).

    Equivalent to solving Poisson then differentiating, but done in one pass:
    E_k = -i k φ_k = -i ρ_k / (ε0 k).
    """
    rho = np.asarray(rho, dtype=float)
    n = rho.size
    k = 2.0 * np.pi * np.fft.fftfreq(n, d=dx)
    rho_k = np.fft.fft(rho)
    E_k = np.zeros_like(rho_k)
    nz = k != 0.0
    E_k[nz] = -1j * rho_k[nz] / (eps0 * k[nz])
    return np.fft.ifft(E_k).real


def solve_poisson_2d(rho, dx, dy, eps0=EPSILON_0):
    """Periodic 2-D Poisson solve. `rho` has shape (nx, ny); returns φ, same shape."""
    rho = np.asarray(rho, dtype=float)
    nx, ny = rho.shape
    kx = 2.0 * np.pi * np.fft.fftfreq(nx, d=dx)
    ky = 2.0 * np.pi * np.fft.fftfreq(ny, d=dy)
    KX, KY = np.meshgrid(kx, ky, indexing="ij")
    k2 = KX ** 2 + KY ** 2
    rho_k = np.fft.fft2(rho)
    phi_k = np.zeros_like(rho_k)
    nz = k2 != 0.0
    phi_k[nz] = rho_k[nz] / (eps0 * k2[nz])
    return np.fft.ifft2(phi_k).real

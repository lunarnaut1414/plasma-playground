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
import scipy.sparse as sp
import scipy.sparse.linalg as spla

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


def grad_shafranov_solve(R, Z, source, boundary=0.0):
    """Solve the fixed-boundary Grad-Shafranov / toroidal elliptic problem

        Δ*ψ = source,     Δ* = ∂²/∂R² - (1/R) ∂/∂R + ∂²/∂Z²

    on a rectangular (R, Z) grid with Dirichlet boundary ψ = `boundary` on the
    grid edges. Δ* is the operator at the heart of tokamak equilibrium: contours
    of the solution ψ are the magnetic flux surfaces.

    Parameters
    ----------
    R, Z : 1-D arrays of grid coordinates (uniform spacing). R must be > 0.
    source : (nR, nZ) array — the right-hand side Δ*ψ. For a real equilibrium
        this is -μ0 R² p'(ψ) - F F'(ψ); for the linear Solov'ev case it is a
        simple function of (R, Z).
    boundary : scalar, or (nR, nZ) array of Dirichlet values (only the edge
        entries are used). Default 0.

    Returns
    -------
    psi : (nR, nZ) ndarray.
    """
    R = np.asarray(R, dtype=float)
    Z = np.asarray(Z, dtype=float)
    source = np.asarray(source, dtype=float)
    nR, nZ = R.size, Z.size
    bc = np.broadcast_to(boundary, (nR, nZ))
    dR = R[1] - R[0]
    dZ = Z[1] - Z[0]
    N = nR * nZ

    def idx(i, j):
        return i * nZ + j

    A = sp.lil_matrix((N, N))
    b = np.zeros(N)
    inv_dR2 = 1.0 / dR**2
    inv_dZ2 = 1.0 / dZ**2

    for i in range(nR):
        for j in range(nZ):
            k = idx(i, j)
            if i == 0 or i == nR - 1 or j == 0 or j == nZ - 1:
                A[k, k] = 1.0          # Dirichlet edge
                b[k] = bc[i, j]
                continue
            grad = 1.0 / (R[i] * 2.0 * dR)        # from the -(1/R) ∂/∂R term
            A[k, idx(i + 1, j)] = inv_dR2 - grad
            A[k, idx(i - 1, j)] = inv_dR2 + grad
            A[k, idx(i, j + 1)] = inv_dZ2
            A[k, idx(i, j - 1)] = inv_dZ2
            A[k, k] = -2.0 * inv_dR2 - 2.0 * inv_dZ2
            b[k] = source[i, j]

    psi = spla.spsolve(A.tocsr(), b)
    return psi.reshape(nR, nZ)

"""V4 & V10 — spectral Poisson solvers.

A spectral solver handles every resolved Fourier mode exactly, so for a source
built from grid-resolved sinusoids the numerical potential matches the analytic
one to round-off — a much stronger statement than the 2nd-order convergence a
finite-difference solver would give.
"""

import numpy as np

from plasmaplay.solvers import (
    grad_shafranov_solve,
    solve_efield_1d,
    solve_poisson_1d,
    solve_poisson_2d,
)

EPS = 1.0  # normalized units keep the test about the operator, not constants


# --- V4: 1-D --------------------------------------------------------------

def test_v4_single_mode_potential():
    # rho = rho0 sin(k x)  ->  phi = rho0/(eps0 k^2) sin(k x)
    L, n = 2 * np.pi, 64
    x = np.linspace(0.0, L, n, endpoint=False)
    dx = L / n
    for m in (1, 2, 5):
        k = 2 * np.pi * m / L
        rho = np.sin(k * x)
        phi = solve_poisson_1d(rho, dx, eps0=EPS)
        phi_exact = np.sin(k * x) / (EPS * k**2)
        np.testing.assert_allclose(phi, phi_exact, atol=1e-12)


def test_v4_multimode_potential():
    # superposition of resolved modes is still exact
    L, n = 2 * np.pi, 128
    x = np.linspace(0.0, L, n, endpoint=False)
    dx = L / n
    k1, k2, k3 = (2 * np.pi * m / L for m in (1, 3, 7))
    rho = 2.0 * np.sin(k1 * x) - 1.5 * np.cos(k2 * x) + 0.5 * np.sin(k3 * x)
    phi = solve_poisson_1d(rho, dx, eps0=EPS)
    phi_exact = (2.0 * np.sin(k1 * x) / k1**2
                 - 1.5 * np.cos(k2 * x) / k2**2
                 + 0.5 * np.sin(k3 * x) / k3**2) / EPS
    np.testing.assert_allclose(phi, phi_exact, atol=1e-12)


def test_v4_efield_matches_minus_grad_phi():
    # E = -dphi/dx ; for rho = sin(kx), E = -cos(kx)*rho0/(eps0 k)
    L, n = 2 * np.pi, 64
    x = np.linspace(0.0, L, n, endpoint=False)
    dx = L / n
    k = 2 * np.pi * 3 / L
    rho = np.sin(k * x)
    E = solve_efield_1d(rho, dx, eps0=EPS)
    E_exact = -np.cos(k * x) / (EPS * k)
    np.testing.assert_allclose(E, E_exact, atol=1e-12)


def test_v4_zero_mean_potential():
    L, n = 2 * np.pi, 64
    x = np.linspace(0.0, L, n, endpoint=False)
    rho = np.sin(2 * np.pi * x / L) + 0.3      # nonzero-mean source
    phi = solve_poisson_1d(rho, L / n, eps0=EPS)
    assert abs(phi.mean()) < 1e-12             # gauge: zero spatial mean


# --- V10: 2-D -------------------------------------------------------------

def test_v10_single_mode_potential_2d():
    # rho = sin(kx x) sin(ky y) -> phi = rho/(eps0 (kx^2+ky^2))
    L, n = 2 * np.pi, 48
    dx = L / n
    x = np.linspace(0.0, L, n, endpoint=False)
    X, Y = np.meshgrid(x, x, indexing="ij")
    for (mx, my) in ((1, 1), (2, 3), (4, 1)):
        kx, ky = 2 * np.pi * mx / L, 2 * np.pi * my / L
        rho = np.sin(kx * X) * np.sin(ky * Y)
        phi = solve_poisson_2d(rho, dx, dx, eps0=EPS)
        phi_exact = rho / (EPS * (kx**2 + ky**2))
        np.testing.assert_allclose(phi, phi_exact, atol=1e-12)


# --- V12: Grad-Shafranov solver (method of manufactured solutions) --------

# ψ = R^4 + Z^4 + R^2 Z^2 is an exact, smooth function whose Δ* is known:
#   Δ*(R^4)   = 8 R^2,   Δ*(Z^4) = 12 Z^2,   Δ*(R^2 Z^2) = 2 R^2
# so  Δ*ψ = 10 R^2 + 12 Z^2.  (A Solov'ev equilibrium is the physical instance
# of the same idea; this manufactured solution checks the operator cleanly.)

def _psi_exact(R, Z):
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    return RR**4 + ZZ**4 + RR**2 * ZZ**2


def _gs_source(R, Z):
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    return 10.0 * RR**2 + 12.0 * ZZ**2


def test_v12_grad_shafranov_manufactured_convergence():
    def max_error(n):
        R = np.linspace(1.0, 2.0, n)      # R > 0 (avoid the 1/R singularity)
        Z = np.linspace(-1.0, 1.0, n)
        exact = _psi_exact(R, Z)
        psi = grad_shafranov_solve(R, Z, _gs_source(R, Z), boundary=exact)
        return np.max(np.abs(psi - exact)), 1.0 / (n - 1)

    errs, hs = zip(*[max_error(n) for n in (21, 41, 81)])
    slope = np.polyfit(np.log(hs), np.log(errs), 1)[0]
    assert 1.8 < slope < 2.2              # 2nd-order accurate Δ* discretization


def test_v12_recovers_exact_solution_on_fine_grid():
    n = 81
    R = np.linspace(1.0, 2.0, n)
    Z = np.linspace(-1.0, 1.0, n)
    exact = _psi_exact(R, Z)
    psi = grad_shafranov_solve(R, Z, _gs_source(R, Z), boundary=exact)
    rel = np.max(np.abs(psi - exact)) / np.max(np.abs(exact))
    assert rel < 1e-3

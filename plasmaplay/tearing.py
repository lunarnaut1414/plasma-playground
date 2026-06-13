"""Linear resistive tearing mode — a reduced-MHD instability (rung T4).

T0–T3 used *prescribed* fields (an equilibrium, a perturbation). T4 is the first
**self-consistent instability**: an equilibrium current sheet that is linearly
unstable to resistive reconnection, growing a magnetic island on its own. We do
the canonical, laptop-scale version — the **slab tearing mode** in reduced MHD —
which has clean analytic checks the full toroidal problem (JOREK/NIMROD) does not.

Geometry: a Harris current sheet, equilibrium field ``B_y(x) = tanh(x/a)`` with a
strong constant guide field ``B_z`` (the reduced-MHD ordering). A perturbation
∝ exp(i k y + γ t) reconnects the field at the neutral line x = 0.

Two pieces, two validations:

1. **Outer (ideal) region & Δ'** — `delta_prime_slab`. Away from the resistive
   layer the marginal ideal equation is the Newcomb equation
        ψ'' − [k² + B_y''/B_y] ψ = 0,   here  B_y''/B_y = −(2/a²) sech²(x/a),
   a Pöschl–Teller potential with the exact decaying solution
        ψ(x) = e^{−k|x|} (1 + tanh(|x|/a)/(k a)).
   Its jump in logarithmic derivative across the sheet is the **tearing stability
   index** Δ' = [ψ'(0⁺) − ψ'(0⁻)]/ψ(0) = (2/a)(1/(k a) − k a). Unstable (Δ' > 0)
   only for **k a < 1**. The numeric Newcomb integrator must reproduce this.

2. **Resistive growth rate & the S^(−3/5) law** — `tearing_growth_rate`. The full
   linear resistive-RMHD eigenvalue problem on an x-grid gives the growth rate γ.
   In the constant-ψ regime Furth–Killeen–Rosenbluth predict
        γ τ_A ∝ Δ'^{4/5} S^{−3/5},   S = τ_R/τ_A (Lundquist number),
   so at fixed k the growth rate must scale as **γ ∝ S^{−3/5}**.

All quantities are normalised: lengths to the sheet width a, B to its asymptotic
value, time to the Alfvén time τ_A = a/v_A, and resistivity η = 1/S.
"""

from __future__ import annotations

import numpy as np
from scipy.sparse import diags, eye
from scipy.sparse import bmat as sp_bmat
from scipy.sparse.linalg import eigs


def harris_By(x, a=1.0):
    """Equilibrium sheet field B_y(x) = tanh(x/a)."""
    return np.tanh(np.asarray(x, dtype=float) / a)


def harris_By_pp(x, a=1.0):
    """B_y''(x) = -(2/a²) sech²(x/a) tanh(x/a)  (the equilibrium current gradient)."""
    x = np.asarray(x, dtype=float)
    return -(2.0 / a**2) * (1.0 / np.cosh(x / a) ** 2) * np.tanh(x / a)


def delta_prime_analytic(k, a=1.0):
    """Exact slab Δ' = (2/a)(1/(k a) − k a) for the tanh sheet (FKR)."""
    return (2.0 / a) * (1.0 / (k * a) - k * a)


def _laplacian_1d(N, dx, k):
    """Sparse d²/dx² − k² with Dirichlet ends (interior N pts)."""
    main = np.full(N, -2.0 / dx**2 - k**2)
    off = np.full(N - 1, 1.0 / dx**2)
    return diags([off, main, off], [-1, 0, 1], format="csr")


def delta_prime_slab(k, a=1.0, L=12.0, N=4001):
    """Numeric Δ' by integrating the Newcomb equation in from both edges.

    Solves ψ'' = [k² + B_y''/B_y] ψ on (0, L] inward to the sheet, with the
    decaying boundary condition ψ ~ e^{−k x} at x = L, then forms the jump in
    ψ'/ψ across x = 0 using the even symmetry of the tearing eigenfunction
    (ψ(−x) = ψ(x) ⇒ Δ' = 2 ψ'(0⁺)/ψ(0)). Returns Δ'·… in 1/length units.
    """
    x = np.linspace(0.0, L, N)
    dx = x[1] - x[0]
    By = harris_By(x, a)
    Bypp = harris_By_pp(x, a)
    # potential Q(x) = k² + B_y''/B_y; at x→0 both →0 so Q→ k² − 2/a² (finite)
    with np.errstate(divide="ignore", invalid="ignore"):
        Q = k**2 + np.where(np.abs(By) > 1e-12, Bypp / By, -2.0 / a**2)
    # integrate inward from the edge (decaying solution) using the ψ'' = Q ψ ODE
    psi = np.zeros(N)
    psi[-1] = np.exp(-k * x[-1])
    psi[-2] = np.exp(-k * x[-2])           # seed the e^{-kx} decay at the edge
    for i in range(N - 2, 0, -1):          # Numerov-free centered scheme
        psi[i - 1] = 2 * psi[i] - psi[i + 1] + dx**2 * Q[i] * psi[i]
    # ψ'(0⁺) by a one-sided derivative; Δ' = 2 ψ'(0⁺)/ψ(0) (even mode)
    dpsi0 = (-3 * psi[0] + 4 * psi[1] - psi[2]) / (2 * dx)
    return 2.0 * dpsi0 / psi[0]


def tearing_growth_rate(k, S, a=1.0, L=8.0, N=2400, return_mode=False):
    """Most-unstable linear resistive tearing growth rate γ (Alfvén units).

    Solves the linearised reduced-MHD eigenvalue problem for a mode
    ∝ exp(i k y + γ t) on x ∈ [−L, L] (Dirichlet ends). With φ = i·φ̂ the system
    is real:

        γ ψ   = −k B_y φ̂ + (1/S) (ψ'' − k² ψ)
        γ Δ*φ̂ =  k B_y Δ*ψ − k B_y'' ψ ,      Δ* ≡ d²/dx² − k²

    posed as a generalised eigenproblem γ M v = A v. Built sparse so N can be
    large enough to resolve the resistive layer, whose width δ ~ a S^{−2/5}
    shrinks with S — under-resolving δ flattens the γ(S) slope, so a fine grid is
    essential to see the S^{−3/5} law. The tearing branch is the eigenvalue of
    largest real part (real, purely growing); `eigs(..., which='LR')` targets it.

    Returns γ (float), or (γ, x, ψ, φ̂) if ``return_mode``.
    """
    x = np.linspace(-L, L, N + 2)[1:-1]    # interior points (Dirichlet at ±L)
    dx = x[1] - x[0]
    Db = diags(harris_By(x, a))
    Dbpp = diags(harris_By_pp(x, a))
    Lm = _laplacian_1d(N, dx, k)           # Δ* = d²/dx² − k²
    Ieye = eye(N, format="csr")

    # state v = [ψ ; φ̂];  γ M v = A v
    A = sp_bmat([[(1.0 / S) * Lm,           -k * Db],
                 [k * (Db @ Lm) - k * Dbpp,  None]], format="csr")
    M = sp_bmat([[Ieye, None], [None, Lm]], format="csr")

    # The tearing eigenvalue is small, real, positive — *interior* to the spectrum
    # (the fine grid carries large spurious modes), so 'largest real part' grabs
    # the wrong mode. Target it by shift-invert near the FKR estimate
    # γ ≈ 0.55 (Δ')^{4/5} S^{−3/5}; among the localized eigenvalues take the
    # most-unstable nearly-real one.
    dprime = delta_prime_analytic(k, a)
    sigma = 0.55 * abs(dprime) ** 0.8 * S ** (-0.6) if dprime > 0 else 1e-6
    w, V = eigs(A, k=6, M=M, sigma=sigma, which="LM", maxiter=10000)
    real_mask = np.abs(w.imag) < 1e-3 * (np.abs(w.real) + 1e-12)
    cand = w[real_mask] if real_mask.any() else w
    idx_local = int(np.argmax(cand.real))
    gamma = float(cand[idx_local].real)
    if not return_mode:
        return gamma
    idx = int(np.argmin(np.abs(w - cand[idx_local])))
    vec = V[:, idx].real
    return gamma, x, vec[:N], vec[N:]

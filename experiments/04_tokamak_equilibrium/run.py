"""Experiment 04 (F1) — Tokamak equilibrium via Grad-Shafranov.

What shape does a confined plasma settle into? Force balance ∇p = J×B in an
axisymmetric torus collapses to a single PDE for the poloidal flux function ψ,
the Grad-Shafranov equation:

    Δ*ψ = -μ0 R² p'(ψ) - F F'(ψ),     Δ* = ∂²/∂R² - (1/R)∂/∂R + ∂²/∂Z²

Contours of ψ are the magnetic flux surfaces. Here we take the *Solov'ev* case —
linear profiles, so the right-hand side is a simple function of (R, Z) and the
solve is one linear elliptic problem (no iteration). We solve it on a fixed
rectangular boundary and look at the equilibrium that emerges.

The headline feature: the **Shafranov shift** — the magnetic axis (the ψ
extremum) sits *outboard* of the geometric center, pushed there by the toroidal
1/R term. Real tokamaks show exactly this.

This is F1 (fixed-boundary, linear profiles). F2 uses FreeGS for free-boundary,
shaped, diverted plasmas. See PLAN.md.

Run:
    python run.py [--save]
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import plotting
from plasmaplay.solvers import grad_shafranov_solve

# Geometry (normalized: R0 = 1 is the geometric center of the box)
R0 = 1.0
A_MINOR = 0.45


def main(save: bool = False, n=141):
    print("=" * 64)
    print("Experiment 04 (F1) — tokamak equilibrium (Grad-Shafranov / Solov'ev)")
    print("=" * 64)

    R = np.linspace(R0 - 0.55, R0 + 0.55, n)
    Z = np.linspace(-0.7, 0.7, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")

    # Solov'ev linear profiles -> source = -μ0 R² p' - F F' = -(c_p R² + c_0).
    # Both coefficients positive => a peaked ψ (a confined current channel).
    c_p, c_0 = 1.0, 1.0
    source = -(c_p * RR**2 + c_0)

    psi = grad_shafranov_solve(R, Z, source, boundary=0.0)

    # Magnetic axis = ψ extremum (here a maximum).
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    R_axis, Z_axis = R[i], Z[j]
    shift = R_axis - R0
    print(f"  magnetic axis at R = {R_axis:.4f}  (geometric center R0 = {R0:.3f})")
    print(f"  Shafranov shift  ΔR = {shift:+.4f}  ({100*shift/A_MINOR:.1f}% of a)")

    fig, ax = plt.subplots(figsize=(6.5, 7.5))
    levels = np.linspace(0, psi.max(), 14)
    cf = ax.contourf(RR, ZZ, psi, levels=levels, cmap="viridis")
    ax.contour(RR, ZZ, psi, levels=levels, colors="white", linewidths=0.5)
    ax.contour(RR, ZZ, psi, levels=[0], colors="k", linewidths=1.2)  # boundary ψ=0
    ax.plot(R_axis, Z_axis, "rx", ms=12, mew=2, label="magnetic axis")
    ax.axvline(R0, color="orange", ls="--", lw=1, label="geometric center R₀")
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
    ax.set_title(f"Grad-Shafranov flux surfaces\nShafranov shift ΔR = {shift:+.3f} m")
    ax.legend(loc="upper right", fontsize=8)
    fig.colorbar(cf, ax=ax, label="poloidal flux ψ", shrink=0.8)
    fig.tight_layout()

    if save:
        out = plotting.ensure_outputs_dir(__file__)
        path = out / "grad_shafranov_equilibrium.png"
        fig.savefig(path, dpi=150)
        print(f"\n  saved -> {path}")
    plt.show()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--save", action="store_true")
    args = p.parse_args()
    main(save=args.save)

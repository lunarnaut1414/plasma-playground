"""T1 — field-line topology: the q-profile & a toroidal Poincaré section.

Rung T1 of `docs/3D_TOKAMAK_GUIDE.md`. With the T0 field in hand, tracing field
lines gives the two central pieces of tokamak topology, almost for free:

  1. **Nested flux surfaces** — a toroidal Poincaré section (puncture the φ = 0
     half-plane once per toroidal turn). Lines on good surfaces draw closed,
     nested curves: the poloidal cross-sections of the flux surfaces.
  2. **The safety factor q(r)** — toroidal turns per poloidal turn of a field
     line. q rises outward (magnetic shear); rational q = m/n surfaces are where
     islands will form in T3.

Two equilibria are shown side by side:
  * left  — a large-aspect *circular* equilibrium with a **sheared** q-profile
    q(r) ≈ q₀(1+(r/a)²) (from ψ ∝ ln(1+(r/a)²)). q is generically irrational, so
    each field line densely fills its circle — clean nested circular surfaces —
    and it crosses the rational q = 1, 3/2, 2 surfaces where T3 will grow islands.
  * right — the experiment-04 *Solov'ev* equilibrium, which has real **shear**:
    q climbs from the axis to the edge.

  (A *uniform* q = 2 circular case — the analytic benchmark in the tests — would
  put every surface on a rational q, so each line closes after 2 toroidal turns
  and punctures φ=0 at only 2 points: correct, but it draws dots, not circles.
  See `docs/T1_QPROFILE_POINCARE.md` for why the viz uses a sheared profile.)

Run:
    MPLBACKEND=Agg python tokamak_t1_viz.py     # headless -> outputs/tokamak_t1.png
    python tokamak_t1_viz.py                     # also show the window
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay.solvers import grad_shafranov_solve
from plasmaplay.tokamak import (
    equilibrium_field,
    safety_factor,
    toroidal_poincare,
    vacuum_F,
)

R0, B0 = 1.0, 2.0


def sheared_circular_equilibrium(Rc=10.0, B0c=1.0, a=1.0, q0=1.0, n=261):
    """Large-aspect circular equilibrium with a sheared q-profile.

    Pick B_θ(r) so that q(r) = r B_φ/(R B_θ) ≈ q₀(1+(r/a)²). With B_φ ≈ B0c
    (large aspect) that needs dψ/dr = R B_θ ≈ B0c r /(q₀(1+(r/a)²)), which
    integrates to ψ(r) = (B0c a²/2q₀) ln(1+(r/a)²). q is irrational on almost
    every surface, so field lines fill their circles densely.
    """
    R = np.linspace(Rc - 1.2, Rc + 1.2, n)
    Z = np.linspace(-1.2, 1.2, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    r2 = (RR - Rc) ** 2 + ZZ ** 2
    psi = (B0c * a**2 / (2.0 * q0)) * np.log1p(r2 / a**2)
    field = equilibrium_field(R, Z, psi, vacuum_F(Rc, B0c))

    def q_analytic(r):
        return q0 * (1.0 + (r / a) ** 2)

    return field, (Rc, 0.0), q_analytic


def solovev_equilibrium(n=161):
    """The experiment-04 Solov'ev equilibrium — real magnetic shear."""
    R = np.linspace(R0 - 0.55, R0 + 0.55, n)
    Z = np.linspace(-0.7, 0.7, n)
    RR, _ = np.meshgrid(R, Z, indexing="ij")
    psi = grad_shafranov_solve(R, Z, -(RR**2 + 1.0), boundary=0.0)
    field = equilibrium_field(R, Z, psi, vacuum_F(R0, B0))
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    return field, (R[i], Z[j]), (R, Z, psi)


def main(save=True):
    fig = plt.figure(figsize=(13, 9))

    # ============ LEFT: sheared circular equilibrium ======================
    cfield, caxis, q_of_r = sheared_circular_equilibrium()

    ax = fig.add_subplot(2, 2, 1)
    radii = [0.25, 0.45, 0.65, 0.85, 1.05]
    for r in radii:
        pc = toroidal_poincare(cfield, (caxis[0] + r, 0.0), n_punctures=130, ds=0.06)
        ax.plot(pc[:, 0], pc[:, 1], ".", ms=1.6)
    ax.plot(*caxis, "rx", ms=10, mew=2)
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
    ax.set_title("Poincaré (φ=0) — nested circular flux surfaces\n(sheared circular equilibrium)")

    ax = fig.add_subplot(2, 2, 3)
    rr = np.linspace(0.2, 1.1, 11)
    qc = [abs(safety_factor(cfield, (caxis[0] + r, 0.0), caxis, n_poloidal=6, ds=0.05))
          for r in rr]
    ax.plot(rr, qc, "o-", label="measured q", zorder=3)
    rfine = np.linspace(0.15, 1.15, 100)
    ax.plot(rfine, q_of_r(rfine), "k--", lw=1, label="analytic q₀(1+(r/a)²)")
    for m_n, lbl in [(1.0, "q=1"), (1.5, "q=3/2"), (2.0, "q=2")]:
        ax.axhline(m_n, color="0.6", ls=":", lw=0.9)
        ax.text(1.12, m_n, lbl, fontsize=7, va="center", color="0.4")
    ax.set_xlabel("minor radius r [m]"); ax.set_ylabel("q")
    ax.set_title("q-profile — rises outward; rational surfaces marked")
    ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.3)

    # ============ RIGHT: Solov'ev equilibrium, real shear =================
    sfield, saxis, (R, Z, psi) = solovev_equilibrium()
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")

    ax = fig.add_subplot(2, 2, 2)
    ax.contour(RR, ZZ, psi, levels=np.linspace(psi.min(), psi.max(), 14),
               colors="0.7", linewidths=0.5)
    for r in (0.08, 0.16, 0.26, 0.36):
        pc = toroidal_poincare(sfield, (saxis[0] + r, saxis[1]),
                               n_punctures=120, ds=0.02)
        ax.plot(pc[:, 0], pc[:, 1], ".", ms=2.0)
    ax.plot(*saxis, "rx", ms=10, mew=2)
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
    ax.set_title("Poincaré (φ=0) — Solov'ev surfaces\n(shaped, contours = ψ)")

    ax = fig.add_subplot(2, 2, 4)
    rs = np.linspace(0.05, 0.4, 9)
    qs = [abs(safety_factor(sfield, (saxis[0] + r, saxis[1]), saxis,
                            n_poloidal=8, ds=0.02))
          for r in rs]
    ax.plot(rs, qs, "s-", color="C3")
    ax.set_xlabel("minor radius r [m]"); ax.set_ylabel("q")
    ax.set_title("q-profile — rises outward (magnetic shear)")
    ax.grid(alpha=0.3)

    fig.suptitle("T1 — field-line topology: nested flux surfaces & the safety factor q(r)",
                 fontsize=14, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    print(f"sheared circular: q rises {qc[0]:.2f} -> {qc[-1]:.2f}; "
          f"analytic q(a)={q_of_r(1.0):.2f} (crosses q=1,3/2,2 — T3 island surfaces)")
    print(f"Solov'ev: axis R={saxis[0]:.3f} Z={saxis[1]:.3f}; "
          f"q rises {qs[0]:.2f} -> {qs[-1]:.2f} (shear)")

    if save:
        out = Path(__file__).resolve().parent / "outputs"
        out.mkdir(exist_ok=True)
        path = out / "tokamak_t1.png"
        fig.savefig(path, dpi=130)
        print(f"saved -> {path}")
    plt.show()


if __name__ == "__main__":
    main()

"""T4 — the linear resistive tearing mode (reduced MHD).

Rung T4 of `docs/3D_TOKAMAK_GUIDE.md` — the first *self-consistent* instability
of the ladder. A Harris current sheet B_y = tanh(x/a) reconnects at its neutral
line, growing a magnetic island. Four panels:

  1. equilibrium: the sheet field B_y(x) and its current −B_y''(x);
  2. Δ'(k): numeric Newcomb solve vs the exact 2(1/ka − ka)/a; unstable for ka<1;
  3. the γ ∝ S^(−3/5) Furth–Killeen–Rosenbluth law (log-log);
  4. the tearing eigenfunction ψ₁(x), and the reconnected magnetic island it makes.

Run:
    MPLBACKEND=Agg python tokamak_t4_viz.py    # headless -> outputs/tokamak_t4.png
    python tokamak_t4_viz.py                    # also show the window
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay.tearing import (
    delta_prime_analytic,
    delta_prime_slab,
    harris_By,
    harris_By_pp,
    tearing_growth_rate,
)


def main(save=True):
    fig = plt.figure(figsize=(13, 9))

    # --- panel 1: the Harris equilibrium ---------------------------------
    ax = fig.add_subplot(2, 2, 1)
    x = np.linspace(-5, 5, 400)
    ax.plot(x, harris_By(x), "C0", label="B_y = tanh(x/a)")
    ax.plot(x, -harris_By_pp(x), "C3", label="current  J = −B_y''")
    ax.axhline(0, color="0.7", lw=0.6); ax.axvline(0, color="0.7", lw=0.6)
    ax.set_xlabel("x / a"); ax.set_title("Harris current sheet (equilibrium)")
    ax.legend(fontsize=8)

    # --- panel 2: Δ'(k) numeric vs analytic ------------------------------
    ax = fig.add_subplot(2, 2, 2)
    kas = np.linspace(0.3, 1.5, 13)
    dn = [delta_prime_slab(ka) for ka in kas]
    ax.plot(kas, delta_prime_analytic(kas), "k-", lw=1.5, label="analytic 2(1/ka − ka)/a")
    ax.plot(kas, dn, "C2o", ms=4, label="numeric Newcomb")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.axvline(1.0, color="C3", ls=":", lw=1)
    ax.text(1.02, ax.get_ylim()[1] * 0.6, "ka=1\nthreshold", fontsize=7, color="C3")
    ax.set_xlabel("k a"); ax.set_ylabel("Δ'·a")
    ax.set_title("Tearing index Δ' — unstable where Δ' > 0 (ka < 1)")
    ax.legend(fontsize=8)

    # --- panel 3: γ ∝ S^(−3/5) -------------------------------------------
    ax = fig.add_subplot(2, 2, 3)
    Ss = np.array([1e4, 3e4, 1e5, 3e5, 1e6])
    gam = np.array([tearing_growth_rate(0.8, S, N=3000) for S in Ss])
    slope = np.polyfit(np.log(Ss), np.log(gam), 1)[0]
    ax.loglog(Ss, gam, "C0o", ms=6, label=f"computed γ  (fit slope {slope:.2f})")
    ref = gam[0] * (Ss / Ss[0]) ** (-0.6)
    ax.loglog(Ss, ref, "k--", lw=1, label="∝ S^(−3/5)  (FKR)")
    ax.set_xlabel("Lundquist number S"); ax.set_ylabel("growth rate γ τ_A")
    ax.set_title("Resistive tearing growth rate (ka=0.8)")
    ax.legend(fontsize=8)

    # --- panel 4: eigenfunction + the reconnected island -----------------
    ax = fig.add_subplot(2, 2, 4)
    g, xg, psi1, phih = tearing_growth_rate(0.8, 1e5, N=2000, return_mode=True)
    psi1 = np.convolve(psi1, np.ones(9) / 9, mode="same")   # de-noise grid-scale ripple
    psi1 = psi1 / np.max(np.abs(psi1)) * np.sign(psi1[np.argmax(np.abs(psi1))])
    a = 1.0
    k = 0.8
    # total helical flux ψ0(x) + ε ψ1(x) cos(k y); ψ0 = −a ln cosh(x/a)
    yy = np.linspace(-np.pi / k, np.pi / k, 200)
    X, Y = np.meshgrid(xg, yy, indexing="ij")
    psi0 = -a * np.log(np.cosh(xg / a))
    eps = 0.25
    PSI = psi0[:, None] + eps * (psi1[:, None] * np.cos(k * Y))
    ax.contour(Y, X, PSI, levels=24, colors="C0", linewidths=0.6)
    ax.plot(yy, 0 * yy, "C3:", lw=0.8)
    ax.set_xlabel("y"); ax.set_ylabel("x / a")
    ax.set_ylim(-2.5, 2.5)
    ax.set_title("Reconnected island (ψ₀ + ε ψ₁ cos ky)")
    axin = ax.inset_axes([0.06, 0.66, 0.32, 0.3])
    axin.plot(xg, psi1, "C2"); axin.set_xlim(-4, 4)
    axin.set_title("ψ₁(x)", fontsize=7); axin.tick_params(labelsize=5)

    fig.suptitle("T4 — the linear resistive tearing mode: Δ', the S^(−3/5) law, and the island",
                 fontsize=13, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    print(f"Δ'(numeric) at ka=0.5: {delta_prime_slab(0.5):.4f}  (analytic {delta_prime_analytic(0.5):.4f})")
    print(f"γ(S) fit slope: {slope:.3f}  (FKR −0.60)")

    if save:
        out = Path(__file__).resolve().parent / "outputs"
        out.mkdir(exist_ok=True)
        path = out / "tokamak_t4.png"
        fig.savefig(path, dpi=130)
        print(f"saved -> {path}")
    plt.show()


if __name__ == "__main__":
    main()

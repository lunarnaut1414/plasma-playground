"""T2 — particle orbits in the tokamak field: banana orbits & trapping.

Rung T2 of `docs/3D_TOKAMAK_GUIDE.md`. With the T0 field in hand, launching
particles is pure reuse of the experiment-01 Boris pusher and the experiment-02
guiding-center pusher. The physics on show:

  * |B| ∝ 1/R is larger on the **inboard** side (small R), so the flux surface is
    a magnetic mirror in miniature. A particle with small parallel velocity
    v∥/v can't climb the |B| hill — it **mirror-traps** and bounces, tracing a
    crescent **"banana"** in the poloidal (R, Z) plane.
  * The trapped/passing boundary sits at pitch λ = v∥/v = √(2ε/(1+ε)) on the
    outboard midplane (ε = r/R₀); for an isotropic distribution that fraction
    √(2ε/(1+ε)) ≈ √(2ε) of particles is trapped.
  * The magnetic moment μ = m v⊥²/(2|B|) is an adiabatic invariant — conserved
    over many gyro-orbits along the full Boris trajectory.

Run:
    MPLBACKEND=Agg python tokamak_t2_viz.py    # headless -> outputs/tokamak_t2.png
    python tokamak_t2_viz.py                    # also show the window
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay.constants import e, m_p
from plasmaplay.guiding_center import gc_push, magnetic_moment
from plasmaplay.pushers import boris_push
from plasmaplay.tokamak import equilibrium_field, to_cartesian, vacuum_F

Rc, B0c, a, Bth0 = 10.0, 1.0, 1.0, 0.05      # large-aspect circular equilibrium
ZERO = lambda x: np.zeros(3)                  # noqa: E731  (no electric field)


def build_field(n=201):
    """ψ = C r² circular equilibrium; B_φ = Rc B0c / R is an exact 1/R mirror."""
    C = Rc * Bth0 / (2.0 * a)
    R = np.linspace(Rc - 1.2, Rc + 1.2, n)
    Z = np.linspace(-1.2, 1.2, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    psi = C * ((RR - Rc) ** 2 + ZZ ** 2)
    return equilibrium_field(R, Z, psi, vacuum_F(Rc, B0c))


def launch_gc(field, r, lam, energy_eV=1000.0, T=2.5e-3, nst=5000):
    """Guiding-center orbit from the outboard midplane (Rc+r, 0) at pitch λ=v∥/v."""
    v = np.sqrt(2.0 * energy_eV * e / m_p)
    x0 = to_cartesian(Rc + r, 0.0, 0.0)
    B0 = np.linalg.norm(field(x0))
    vpar0 = lam * v
    vperp = np.sqrt(max(v * v - vpar0 * vpar0, 0.0))
    mu = magnetic_moment(vperp, m_p, B0)
    t, pos, vpar = gc_push(x0, vpar0, mu, e, m_p, ZERO, field, T / nst, nst)
    return t, pos, vpar


def lambda_c_measured(field, r, iters=9, nst=2200):
    """Bisect the trapped/passing boundary pitch at minor radius r."""
    lo, hi = 0.0, 1.0
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        _, _, vpar = launch_gc(field, r, mid, T=2e-3, nst=nst)
        trapped = np.sign(vpar).min() != np.sign(vpar).max()
        lo, hi = (mid, hi) if trapped else (lo, mid)
    return 0.5 * (lo + hi)


def main(save=True):
    field = build_field()
    fig = plt.figure(figsize=(13, 9))

    # ---- panel A: banana (trapped) + passing orbit in the poloidal plane ----
    ax = fig.add_subplot(2, 2, 1)
    r = 0.5
    th = np.linspace(0, 2 * np.pi, 200)
    ax.plot(Rc + r * np.cos(th), r * np.sin(th), "0.8", lw=1, label="flux surface")
    _, pos_t, _ = launch_gc(field, r, 0.12)          # trapped -> banana
    _, pos_p, _ = launch_gc(field, r, 0.55)          # passing -> follows surface
    Rt, Zt = np.hypot(pos_t[:, 0], pos_t[:, 1]), pos_t[:, 2]
    Rp, Zp = np.hypot(pos_p[:, 0], pos_p[:, 1]), pos_p[:, 2]
    ax.plot(Rp, Zp, "C0", lw=0.6, alpha=0.5, label="passing (λ=0.55)")
    ax.plot(Rt, Zt, "C3", lw=1.3, label="trapped — banana (λ=0.12)")
    ax.plot(Rc + r, 0, "k.", ms=8)
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
    ax.set_title("Poloidal projection — the banana orbit")
    ax.legend(fontsize=7, loc="upper right")

    # ---- panel B: bounce — v∥ reverses where |B| peaks (inboard turning pts) --
    ax = fig.add_subplot(2, 2, 2)
    t, pos_t, vpar = launch_gc(field, r, 0.12)
    Bmag = np.array([np.linalg.norm(field(p)) for p in pos_t])
    v = np.sqrt(2.0 * 1000.0 * e / m_p)
    ax.plot(t * 1e3, vpar / v, "C3", label="v∥/v")
    ax.axhline(0, color="0.6", lw=0.8)
    ax.set_xlabel("t [ms]"); ax.set_ylabel("v∥ / v", color="C3")
    ax.tick_params(axis="y", colors="C3")
    ax2 = ax.twinx()
    ax2.plot(t * 1e3, Bmag, "C0", lw=0.8, alpha=0.7)
    ax2.set_ylabel("|B| [T]", color="C0"); ax2.tick_params(axis="y", colors="C0")
    ax.set_title("Bounce: v∥ → 0 at the |B| turning points")

    # ---- panel C: trapped/passing boundary vs ε ---------------------------
    ax = fig.add_subplot(2, 2, 3)
    rs = np.array([0.3, 0.45, 0.6, 0.75, 0.9])
    eps = rs / Rc
    pred = np.sqrt(2 * eps / (1 + eps))
    meas = np.array([lambda_c_measured(field, rr) for rr in rs])
    ax.plot(eps, pred, "k--", label="√(2ε/(1+ε))")
    ax.plot(eps, meas, "C3o", label="measured λ_c (bisection)")
    ax.set_xlabel("ε = r / R₀"); ax.set_ylabel("trapped/passing pitch λ_c")
    ax.set_title("Trapping boundary = trapped fraction (isotropic)")
    ax.legend(fontsize=8); ax.grid(alpha=0.3)

    # ---- panel D: μ adiabatic invariance along the Boris orbit ------------
    ax = fig.add_subplot(2, 2, 4)
    x0 = to_cartesian(Rc + r, 0.0, 0.0)
    Bv = field(x0); B0 = np.linalg.norm(Bv); bhat = Bv / B0
    vpar0 = 0.12 * v; vperp = np.sqrt(v * v - vpar0 * vpar0)
    e1 = np.cross(bhat, [0, 0, 1.0]); e1 /= np.linalg.norm(e1)
    vel = vpar0 * bhat + vperp * e1
    wc = e * B0 / m_p; dt = 2 * np.pi / wc / 30
    tb, posb, velb = boris_push(x0, vel, e, m_p, ZERO, field, dt, 15000)
    s = slice(0, len(posb), 50)
    Bm = np.array([np.linalg.norm(field(p)) for p in posb[s]])
    bh = np.array([field(p) / np.linalg.norm(field(p)) for p in posb[s]])
    vp = velb[s]
    vparr = np.einsum("ij,ij->i", vp, bh)
    mu = m_p * (np.sum(vp**2, axis=1) - vparr**2) / (2 * Bm)
    ax.plot(tb[s] * 1e6, mu / mu[0], "C2")
    ax.set_xlabel("t [µs]"); ax.set_ylabel("μ / μ₀")
    ax.set_title(f"μ conserved to {100*(mu.max()-mu.min())/mu.mean():.2f}% "
                 f"over {tb[-1]*wc/(2*np.pi):.0f} gyro-orbits")
    ax.grid(alpha=0.3)

    fig.suptitle("T2 — banana orbits, the trapping boundary, and the μ adiabatic invariant",
                 fontsize=14, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    print(f"trapping boundary λ_c: measured {np.round(meas,3)} vs √(2ε/(1+ε)) {np.round(pred,3)}")
    print(f"μ spread over Boris orbit: {100*(mu.max()-mu.min())/mu.mean():.3f}%")

    if save:
        out = Path(__file__).resolve().parent / "outputs"
        out.mkdir(exist_ok=True)
        path = out / "tokamak_t2.png"
        fig.savefig(path, dpi=130)
        print(f"saved -> {path}")
    plt.show()


if __name__ == "__main__":
    main()

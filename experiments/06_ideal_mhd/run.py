"""Experiment 06 — Ideal MHD: waves and shocks.

Forget individual particles — treat the plasma as one conducting fluid threaded
by magnetic field, and solve the ideal MHD conservation laws. Two faces of the
same model:

  alfven   a circularly-polarized Alfvén wave (an exact solution) translates at
           v_A = Bx/√ρ without changing shape — magnetic tension as a restoring
           force. (F1)
  briowu   the Brio-Wu shock tube — the MHD analogue of Sod's problem — develops
           the full nonlinear zoo: fast rarefactions, a slow compound wave, a
           contact, and a slow shock. (F2)

Run:
    python run.py --case briowu [--save]
    python run.py --case alfven
    python run.py --case all --save
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import fvm, plotting

GAMMA = 2.0


def mhd_speeds(rho, p, Bx, By, Bz, gamma=GAMMA):
    """Sound, Alfven (along x), and fast/slow magnetosonic speeds."""
    a = np.sqrt(gamma * p / rho)
    cA = abs(Bx) / np.sqrt(rho)
    a2, b2, bx2 = a**2, (Bx**2 + By**2 + Bz**2) / rho, Bx**2 / rho
    disc = np.sqrt((a2 + b2) ** 2 - 4 * a2 * bx2)
    cf = np.sqrt(0.5 * (a2 + b2 + disc))
    cs = np.sqrt(0.5 * (a2 + b2 - disc))
    return a, cA, cf, cs


def run_alfven(save=False, n=256, t=0.25):
    print("\n--- circularly-polarized Alfven wave (F1) ---")
    L = 1.0
    dx = L / n
    x = (np.arange(n) + 0.5) * dx
    A, k, rho, Bx, p = 0.1, 2 * np.pi / L, 1.0, 1.0, 0.5
    By0, Bz0 = A * np.cos(k * x), A * np.sin(k * x)
    W0 = np.zeros((n, 7))
    W0[:, 0], W0[:, 4] = rho, p
    W0[:, 5], W0[:, 6] = By0, Bz0
    W0[:, 2], W0[:, 3] = -By0 / np.sqrt(rho), -Bz0 / np.sqrt(rho)

    v_A = Bx / np.sqrt(rho)
    Wf = fvm.solve_mhd_1d(W0, dx, Bx=Bx, gamma=GAMMA, t_end=t, bc="periodic")
    print(f"  v_A = {v_A:.3f};  wave should shift by v_A·t = {v_A * t:.3f}")

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(x, By0, label="By at t=0", lw=1.2)
    ax.plot(x, Wf[:, 5], "--", label=f"By at t={t} (shifted by v_A·t)", lw=1.2)
    ax.axvline((v_A * t) % L, color="gray", ls=":", lw=1)
    ax.set_xlabel("x"); ax.set_ylabel("By")
    ax.set_title("Alfvén wave translates at v_A (magnetic tension)")
    ax.legend()
    fig.tight_layout()
    _finish(fig, "alfven_wave.png", save)


def run_briowu(save=False, n=800, t=0.1):
    print("\n--- Brio-Wu shock tube (F2) ---")
    dx = 1.0 / n
    x = (np.arange(n) + 0.5) * dx
    W0 = np.zeros((n, 7))
    left = x < 0.5
    W0[left] = [1.0, 0, 0, 0, 1.0, 1.0, 0]
    W0[~left] = [0.125, 0, 0, 0, 0.1, -1.0, 0]

    a, cA, cf, cs = mhd_speeds(1.0, 1.0, 0.75, 1.0, 0.0)
    print(f"  left-state speeds: sound={a:.3f}  Alfven={cA:.3f}  "
          f"fast={cf:.3f}  slow={cs:.3f}")

    Wf = fvm.solve_mhd_1d(W0, dx, Bx=0.75, gamma=GAMMA, t_end=t)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    panels = [(0, "density ρ"), (1, "x-velocity u"), (5, "By"), (4, "pressure p")]
    for ax, (col, name) in zip(axes.flat, panels):
        ax.plot(x, Wf[:, col], lw=1.0)
        ax.set_xlabel("x"); ax.set_ylabel(name); ax.set_title(name)
        ax.grid(alpha=0.3)
    fig.suptitle(f"Brio-Wu MHD shock tube (t={t}, n={n})", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    _finish(fig, "brio_wu.png", save)


def _finish(fig, name, save):
    if save:
        out = plotting.ensure_outputs_dir(__file__)
        path = out / name
        fig.savefig(path, dpi=150)
        print(f"  saved -> {path}")
    plt.show()


def main(case="briowu", save=False):
    print("=" * 64)
    print("Experiment 06 — ideal MHD (waves and shocks)")
    print("=" * 64)
    if case in ("alfven", "all"):
        run_alfven(save)
    if case in ("briowu", "all"):
        run_briowu(save)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", default="briowu", choices=["alfven", "briowu", "all"])
    p.add_argument("--save", action="store_true")
    args = p.parse_args()
    main(case=args.case, save=args.save)

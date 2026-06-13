"""Experiment 08 — Plasma waves & dispersion (the unifier).

Every plasma wave has a dispersion relation ω(k). This experiment computes the
Langmuir-wave dispersion three ways and shows they agree — closing the loop on
the whole playground:

  kinetic   solve the kinetic dispersion relation (plasma dispersion function Z)
            for the real frequency AND the Landau damping rate. (F2)
  pic       run a warm PIC plasma (experiment 03's solver), FFT the field history
            E(x, t) into the ω–k plane, and watch the simulated dispersion ridge
            land on the analytic curve. (F3 — kinetic theory meets simulation)

Run:
    python run.py --case kinetic [--save]
    python run.py --case pic
    python run.py --case all --save
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import pic, plotting
from plasmaplay.diagnostics import omega_k_spectrum
from plasmaplay.dispersion import bohm_gross, langmuir_dispersion


def run_kinetic(save=False):
    print("\n--- kinetic Langmuir dispersion (F2) ---")
    kld = np.linspace(0.05, 0.6, 40)
    wr = np.array([langmuir_dispersion(k)[0] for k in kld])
    gamma = np.array([langmuir_dispersion(k)[1] for k in kld])
    for k in (0.3, 0.5):
        w, g = langmuir_dispersion(k)
        print(f"  kλ_D={k}:  ω_r={w:.4f} ω_pe   γ={g:.4f} ω_pe")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    a1.plot(kld, wr, label="kinetic ω_r")
    a1.plot(kld, bohm_gross(kld), "--", label="Bohm-Gross (fluid)")
    a1.axhline(1.0, color="gray", ls=":", lw=1, label="cold (ω_pe)")
    a1.set_xlabel("k λ_D"); a1.set_ylabel("ω_r / ω_pe")
    a1.set_title("Langmuir frequency: cold → fluid → kinetic"); a1.legend()

    a2.semilogy(kld, -gamma)
    a2.set_xlabel("k λ_D"); a2.set_ylabel("−γ / ω_pe  (Landau damping)")
    a2.set_title("Landau damping grows steeply with k λ_D")
    a2.grid(which="both", alpha=0.3)
    fig.tight_layout()
    _finish(fig, "kinetic_dispersion.png", save)


def run_pic(save=False, n_particles=200000, n_steps=400):
    print("\n--- dispersion measured from a PIC run (F3) ---")
    rng = np.random.default_rng(0)
    L, NG = 2 * np.pi, 128
    v_th = 0.14                       # λ_D = v_th/√2 ≈ 0.1
    lambda_D = v_th / np.sqrt(2.0)
    # random (noisy) load so thermal fluctuations excite the Langmuir modes
    x, v = pic.load_maxwellian(n_particles, L, v_th, rng=rng, quiet=False)
    q, m = -L / n_particles, L / n_particles
    sim = pic.ElectrostaticPIC1D(L, NG, x, v, q, m, eps0=1.0)

    dt = 0.1
    E_xt = np.empty((n_steps, NG))
    for n in range(n_steps):
        E_xt[n] = sim.step(dt)

    k, omega, power = omega_k_spectrum(E_xt, dx=L / NG, dt=dt)
    print(f"  recorded E(x,t): {E_xt.shape};  λ_D={lambda_D:.3f}")

    # restrict to the band where the Langmuir branch lives (k λ_D up to ~1)
    k_max = 1.2 / lambda_D
    kpos = (k >= 0) & (k <= k_max)
    wpos = (omega >= 0) & (omega <= 3)
    sub = np.log10(power[np.ix_(wpos, kpos)] + 1e-12)
    vmin, vmax = np.percentile(sub, [50, 99.5])

    fig, ax = plt.subplots(figsize=(8, 6))
    extent = [k[kpos].min(), k[kpos].max(), omega[wpos].min(), omega[wpos].max()]
    ax.imshow(sub, origin="lower", aspect="auto", extent=extent,
              cmap="inferno", vmin=vmin, vmax=vmax)
    kk = np.linspace(0.01, k[kpos].max(), 200)
    ax.plot(kk, bohm_gross(kk * lambda_D), "c--", lw=2,
            label="kinetic theory: Bohm-Gross ω(k)")
    ax.set_xlabel("k"); ax.set_ylabel("ω / ω_pe")
    ax.set_title("ω–k spectrum from PIC fields → the Langmuir dispersion ridge")
    ax.legend(loc="upper left")
    fig.tight_layout()
    _finish(fig, "pic_omega_k.png", save)


def _finish(fig, name, save):
    if save:
        out = plotting.ensure_outputs_dir(__file__)
        path = out / name
        fig.savefig(path, dpi=150)
        print(f"  saved -> {path}")
    plt.show()


def main(case="kinetic", save=False):
    print("=" * 64)
    print("Experiment 08 — plasma waves & dispersion")
    print("=" * 64)
    if case in ("kinetic", "all"):
        run_kinetic(save)
    if case in ("pic", "all"):
        run_pic(save)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", default="kinetic", choices=["kinetic", "pic", "all"])
    p.add_argument("--save", action="store_true")
    args = p.parse_args()
    main(case=args.case, save=args.save)

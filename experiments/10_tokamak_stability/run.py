"""Experiment 10 — Tokamak MHD stability: the straight-tokamak (cylinder) limit.

The transport experiment (09) is the *slow* half of a discharge — energy and
particle balance over seconds. This experiment is the *fast* half: the ideal/
resistive **MHD instabilities** that move the plasma fluid on the microsecond
Alfven timescale. We start where it is cleanest — a periodic cylinder (a torus cut
and straightened), the "straight tokamak" — where a current profile gives a safety
factor q(r) and the two organising instabilities appear:

  * the **m=1/n=1 internal kink** — the rigid sideways shift of the core inside the
    q=1 surface, which exists exactly when q(0) < 1. This is the **sawtooth trigger**.

  * the **m/n tearing mode** — reconnection at the q=m/n rational surface, governed
    by the outer Newcomb equation and its stability index Delta' (sign = stability),
    with the resistive growth rate following the FKR gamma ~ S^(-3/5) layer law.

This is rung B1 of the MHD track (NIGHT.md): it lifts the slab tearing of T4
(plasmaplay/tearing.py) onto a real q(r). Later rungs add the nonlinear island
saturation (B2) and the sawtooth cycle (B3), and Track C couples a sawtooth/tearing
event into the transport burn of experiment 09.

Run:
    python run.py [--save]
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import cylinder_mhd as cm, plotting


def main(save=False):
    print("=" * 64)
    print("Experiment 10 — tokamak MHD stability (straight-tokamak / cylinder)")
    print("=" * 64)

    # --- the m=1 internal kink: unstable iff q(0) < 1 (the sawtooth trigger) ---
    print("\n--- m=1/n=1 internal kink (the sawtooth trigger) ---")
    for q0 in (0.7, 0.85, 0.95, 1.05, 1.2):
        r1 = cm.rational_surface(1, 1, q0, nu=1.0)
        flag = "UNSTABLE" if cm.internal_kink_unstable(q0) else "stable  "
        loc = f"q=1 surface at r1 = {r1:.3f}" if r1 is not None else "no q=1 surface"
        print(f"  q(0) = {q0:.2f}: {flag}  ({loc})")

    # --- tearing Delta': sign predicts stability, falls with mode number m ---
    print("\n--- m/n tearing modes: Delta' sign predicts stability ---")
    for (m, n, q0, nu) in [(2, 1, 1.2, 2.0), (3, 1, 1.2, 2.0), (3, 1, 0.9, 3.0)]:
        r_s = cm.rational_surface(m, n, q0, nu)
        dp = cm.delta_prime_cylinder(m, n, q0, nu)
        verdict = "tearing-UNSTABLE" if dp > 0 else "tearing-stable"
        print(f"  m={m}/n={n}, q0={q0}, nu={nu}: r_s={r_s:.3f}, "
              f"Delta'={dp:+.2f} -> {verdict}")
    d2 = cm.delta_prime_cylinder(2, 1, 1.2, 2.0)
    d3 = cm.delta_prime_cylinder(3, 1, 1.2, 2.0)
    print(f"  higher-m is more wall-stabilized: Delta'(m=2)={d2:.2f} > "
          f"Delta'(m=3)={d3:.2f}  ({d2 > d3})")

    # --- resistive growth rate: the FKR S^(-3/5) law (carried over from the slab) ---
    print("\n--- resistive growth rate (FKR layer law) ---")
    for S in (1e4, 1e5, 1e6, 1e7):
        g = cm.fkr_growth_rate(d2, S)
        print(f"  S = {S:.0e}: gamma*tau_A = {g:.3e}")
    ratio = cm.fkr_growth_rate(d2, 1e6) / cm.fkr_growth_rate(d2, 1e5)
    print(f"  gamma(10S)/gamma(S) = {ratio:.4f}  (FKR S^-3/5 -> {10 ** -0.6:.4f})")

    _plot(save)


def _plot(save):
    q0, nu = 0.85, 1.0
    r1 = cm.rational_surface(1, 1, q0, nu)
    rg = np.linspace(1e-3, 1.0, 200)
    xi = cm.internal_kink_xi(rg, r1)
    q = cm.screw_pinch_q(rg, q0, nu)

    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    ax[0].plot(rg, q, color="navy", label="q(r)")
    ax[0].plot(rg, xi, color="crimson", label=r"$\xi_r(r)$ (m=1 kink)")
    ax[0].axhline(1.0, color="0.6", ls=":")
    ax[0].axvline(r1, color="k", ls="--", lw=0.9)
    ax[0].text(r1, 1.7, "q=1", rotation=90, ha="right", va="bottom", fontsize=8)
    ax[0].set(xlabel="r/a", ylim=(0, 2), title=f"Internal kink (q0={q0}): rigid core shift")
    ax[0].legend()

    # Delta' vs q0 for the m=2/n=1 mode (where the q=2 surface exists)
    q0s = np.linspace(0.6, 1.9, 40)
    dps = [cm.delta_prime_cylinder(2, 1, x, 2.0) for x in q0s]
    dps = np.array([d if np.isfinite(d) else np.nan for d in dps])
    ax[1].plot(q0s, dps, color="seagreen")
    ax[1].axhline(0.0, color="k", lw=0.8)
    ax[1].set(xlabel="q(0)", ylabel=r"$\Delta'$ (m=2/n=1)",
              title="Tearing index: sign predicts stability")
    fig.tight_layout()
    if save:
        out = plotting.ensure_outputs_dir(__file__) / "kink_eigenmode.png"
        fig.savefig(out, dpi=130, bbox_inches="tight")
        print(f"\n  saved {out}")
    else:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--save", action="store_true", help="write figures to ./outputs/")
    args = p.parse_args()
    main(save=args.save)

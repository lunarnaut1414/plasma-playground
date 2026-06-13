"""Experiment 02 — Guiding-center drifts.

A magnetized particle is "a fast little circle riding on a slowly drifting
center." This experiment throws away the circle and follows only the center —
the way plasma physicists actually think about confinement.

  drifts    the guiding center in a grad-B field drifts sideways; overlay it on
            the full Boris orbit to see it IS the gyro-averaged motion. (F1/F2)
  adiabatic sweep the gyroradius/scale-length ratio r_L/L and watch the
            guiding-center approximation degrade as it leaves the adiabatic
            regime. (F2)

Run:
    python run.py --case drifts [--save]
    python run.py --case adiabatic
    python run.py --case all --save
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import constants as k
from plasmaplay import fields, guiding_center as gc, plotting, pushers

CHARGE, MASS = k.e, k.m_p


def _boris_and_gc(B0, L, v_perp, n_orbits=30, steps_per_orbit=80):
    omega_c = k.gyrofrequency(CHARGE, B0, MASS)
    dt = (2 * np.pi / omega_c) / steps_per_orbit
    n = n_orbits * steps_per_orbit
    _, pos_full, _ = pushers.boris_push(
        [0, 0, 0], [v_perp, 0, 0], CHARGE, MASS,
        fields.zero_field(), fields.gradient_B(B0, L), dt, n)
    mu = gc.magnetic_moment(v_perp, MASS, B0)
    _, pos_gc, _ = gc.gc_push(
        [0, 0, 0], 0.0, mu, CHARGE, MASS,
        fields.zero_field(), fields.gradient_B(B0, L), dt, n)
    return pos_full, pos_gc


def run_drifts(save=False):
    print("\n--- guiding center vs full orbit in a grad-B field (F1/F2) ---")
    B0, L, v_perp = 1.0, 0.05, 1.0e5
    r_L = k.gyroradius(v_perp, CHARGE, B0, MASS)
    print(f"  r_L = {r_L*1e3:.3f} mm,  L = {L*1e3:.0f} mm,  r_L/L = {r_L/L:.3f}")
    pos_full, pos_gc = _boris_and_gc(B0, L, v_perp)

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.plot(pos_full[:, 0] * 1e3, pos_full[:, 1] * 1e3, lw=0.5,
            color="C0", alpha=0.7, label="full Boris orbit")
    ax.plot(pos_gc[:, 0] * 1e3, pos_gc[:, 1] * 1e3, lw=2.0,
            color="C3", label="guiding center")
    ax.set_xlabel("x [mm]"); ax.set_ylabel("y [mm]")
    ax.set_title("Grad-B drift: the guiding center is the gyro-averaged orbit")
    ax.legend()
    fig.tight_layout()
    _finish(fig, "guiding_center_drift.png", save)


def run_adiabatic(save=False):
    print("\n--- adiabaticity sweep: GC error vs r_L/L (F2) ---")
    B0, L = 1.0, 0.05
    ratios, errors = [], []
    for v_perp in np.logspace(np.log10(2e4), np.log10(3e6), 12):
        r_L = k.gyroradius(v_perp, CHARGE, B0, MASS)
        pos_full, pos_gc = _boris_and_gc(B0, L, v_perp, n_orbits=20)
        per = 80
        drift_full = pos_full[-per:, 1].mean() - pos_full[:per, 1].mean()
        drift_gc = pos_gc[-per:, 1].mean() - pos_gc[:per, 1].mean()
        rel_err = abs(drift_full - drift_gc) / abs(drift_gc)
        ratios.append(r_L / L); errors.append(rel_err)
        print(f"  r_L/L={r_L/L:6.3f}  relative drift error={rel_err:.3e}")

    fig, ax = plt.subplots(figsize=(7, 5))
    ax.loglog(ratios, errors, "o-")
    ax.set_xlabel("adiabaticity parameter  r_L / L")
    ax.set_ylabel("relative GC drift error")
    ax.set_title("Guiding-center approximation degrades as r_L/L → 1")
    ax.grid(which="both", alpha=0.3)
    fig.tight_layout()
    _finish(fig, "adiabaticity_sweep.png", save)


def _finish(fig, name, save):
    if save:
        out = plotting.ensure_outputs_dir(__file__)
        path = out / name
        fig.savefig(path, dpi=150)
        print(f"  saved -> {path}")
    plt.show()


def main(case="drifts", save=False):
    print("=" * 64)
    print("Experiment 02 — guiding-center drifts")
    print("=" * 64)
    if case in ("drifts", "all"):
        run_drifts(save)
    if case in ("adiabatic", "all"):
        run_adiabatic(save)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", default="drifts", choices=["drifts", "adiabatic", "all"])
    p.add_argument("--save", action="store_true")
    args = p.parse_args()
    main(case=args.case, save=args.save)

"""Experiment 07 — An MHD drive for space propulsion.

The repo's namesake goal: accelerate plasma with the J×B (Lorentz) force to make
thrust. Two complementary views of the same physics:

  scaling   self-field MPD thruster — the discharge current crossed with its OWN
            azimuthal field gives the Maecker law T ∝ I² (F0).
  channel   a 1-D applied-field accelerator — a prescribed Lorentz force density
            accelerates the flow along a channel; thrust = B₀ I L, linear in I,
            and the exit velocity sets the specific impulse (F1).

The contrast (T ∝ I² self-field vs T ∝ I applied-field) is the headline.

Run:
    python run.py --case scaling [--save]
    python run.py --case channel
    python run.py --case all --save
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import plotting
from plasmaplay import propulsion as prop

# A representative high-power MPD operating point
R_ANODE, R_CATHODE = 0.05, 0.01     # m
MDOT = 1.0e-4                       # kg/s (~0.1 g/s)


def run_scaling(save=False):
    print("\n--- self-field MPD thrust scaling (F0) ---")
    I = np.linspace(500, 10000, 200)
    T = prop.maecker_thrust(I, R_ANODE, R_CATHODE)
    Isp = prop.specific_impulse(T, MDOT)
    P = prop.jet_power(T, MDOT)

    for Ik in (1e3, 3e3, 1e4):
        Tk = prop.maecker_thrust(Ik, R_ANODE, R_CATHODE)
        print(f"  I={Ik/1e3:4.0f} kA -> T={Tk:6.2f} N  "
              f"Isp={prop.specific_impulse(Tk, MDOT):7.0f} s  "
              f"P_jet={prop.jet_power(Tk, MDOT)/1e3:7.1f} kW")

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    a1.loglog(I / 1e3, T, lw=1.5)
    a1.set_xlabel("current I [kA]"); a1.set_ylabel("thrust T [N]")
    a1.set_title("Maecker self-field thrust  (slope = 2 → T ∝ I²)")
    a1.grid(which="both", alpha=0.3)

    a2.plot(I / 1e3, Isp, label="Isp [s]")
    a2b = a2.twinx()
    a2b.plot(I / 1e3, P / 1e3, "r--", label="jet power [kW]")
    a2.set_xlabel("current I [kA]"); a2.set_ylabel("Isp [s]")
    a2b.set_ylabel("jet power [kW]", color="r")
    a2.set_title(f"Isp and power vs current (ṁ = {MDOT*1e3:.1f} g/s)")
    fig.tight_layout()
    _finish(fig, "mpd_scaling.png", save)


def run_channel(save=False):
    print("\n--- applied-field accelerator channel (F1) ---")
    B0, L, area, mdot = 0.1, 0.2, 1.0e-3, 1.0e-3      # T, m, m², kg/s
    n = 401
    x = np.linspace(0, L, n)
    dx = x[1] - x[0]
    G = mdot / area

    I0 = 2000.0
    f = np.full(n, (I0 / area) * B0)                  # f = j × B₀ (uniform)
    u = prop.channel_velocity(f, dx, G, u_inlet=0.0)
    thrust = mdot * (u[-1] - u[0])
    # energy balance check (ideal cold model -> Lorentz work = KE gain exactly)
    lorentz_work = area * np.trapezoid(f * u, dx=dx)
    ke_gain = 0.5 * mdot * (u[-1] ** 2 - u[0] ** 2)
    print(f"  I={I0/1e3:.1f} kA, B0={B0} T, L={L} m ->")
    print(f"    exit u = {u[-1]/1e3:.1f} km/s   thrust = {thrust:.1f} N   "
          f"Isp = {prop.specific_impulse(thrust, mdot):.0f} s")
    print(f"    thrust check  B0·I·L = {B0*I0*L:.1f} N")
    print(f"    energy balance: Lorentz work {lorentz_work:.3e} W vs "
          f"KE gain {ke_gain:.3e} W  (ratio {lorentz_work/ke_gain:.4f})")

    # thrust vs current: applied-field (∝ I) vs self-field Maecker (∝ I²)
    Ivals = np.linspace(500, 5000, 100)
    T_applied = B0 * Ivals * L
    T_self = prop.maecker_thrust(Ivals, R_ANODE, R_CATHODE)

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))
    a1.plot(x * 100, u / 1e3, lw=1.5)
    a1.set_xlabel("position along channel [cm]"); a1.set_ylabel("flow speed u [km/s]")
    a1.set_title(f"Channel acceleration (I={I0/1e3:.0f} kA, B₀={B0} T)")
    a1.grid(alpha=0.3)

    a2.plot(Ivals / 1e3, T_applied, label="applied-field  T = B₀ I L  (∝ I)")
    a2.plot(Ivals / 1e3, T_self, "--", label="self-field Maecker  (∝ I²)")
    a2.set_xlabel("current I [kA]"); a2.set_ylabel("thrust T [N]")
    a2.set_title("Two thrust scalings"); a2.legend(); a2.grid(alpha=0.3)
    fig.tight_layout()
    _finish(fig, "channel_accelerator.png", save)


def _finish(fig, name, save):
    if save:
        out = plotting.ensure_outputs_dir(__file__)
        path = out / name
        fig.savefig(path, dpi=150)
        print(f"  saved -> {path}")
    plt.show()


def main(case="scaling", save=False):
    print("=" * 64)
    print("Experiment 07 — MHD drive for space propulsion")
    print("=" * 64)
    if case in ("scaling", "all"):
        run_scaling(save)
    if case in ("channel", "all"):
        run_channel(save)


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--case", default="scaling", choices=["scaling", "channel", "all"])
    p.add_argument("--save", action="store_true")
    args = p.parse_args()
    main(case=args.case, save=args.save)

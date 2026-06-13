"""Experiment 01 — Single-particle motion in E and B fields.

The "hello world" of plasma simulation. A single charged particle is the
building block of everything else: a plasma is ~10^20 of these, and the
collective behaviour (waves, instabilities, confinement) emerges from how they
move in fields that they themselves help create.

Here we look at three canonical cases using the Boris pusher:

  1. Uniform B only          -> circular gyro-orbit (Larmor motion)
  2. Crossed E and B         -> gyration + steady E x B drift
  3. Magnetic mirror         -> trapping / reflection of the particle

Run:
    python run.py            # shows plots
    python run.py --save     # also writes figures to ./outputs/

Everything below is plain NumPy + matplotlib via the shared `plasmaplay`
package, so it runs instantly on CPU — no GPU needed.
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import constants as k
from plasmaplay import fields, plotting, pushers


def case_gyro_orbit():
    """Proton in a uniform 1 T field along z. Pure circular Larmor motion."""
    B0 = 1.0  # tesla
    charge, mass = k.e, k.m_p

    v_perp = 1.0e5  # m/s
    omega_c = k.gyrofrequency(charge, B0, mass)
    r_L = k.gyroradius(v_perp, charge, B0, mass)
    T_c = 2 * np.pi / omega_c

    dt = T_c / 60          # resolve the gyro-period well
    n_steps = 60 * 3       # three full orbits

    t, pos, vel = pushers.boris_push(
        position=[0, 0, 0],
        velocity=[v_perp, 0, 0],
        charge=charge, mass=mass,
        E_func=fields.zero_field(),
        B_func=fields.uniform_B(B0),
        dt=dt, n_steps=n_steps,
    )

    # Energy conservation check — Boris should hold this flat to ~machine eps.
    speed2 = np.sum(vel**2, axis=1)
    drift = (speed2.max() - speed2.min()) / speed2.mean()
    print(f"[gyro]   gyrofreq = {omega_c:.3e} rad/s   r_L = {r_L*1e3:.3f} mm")
    print(f"[gyro]   relative energy drift over 3 orbits = {drift:.2e}")

    return t, pos, vel


def case_ExB_drift():
    """Crossed E and B: the guiding centre drifts at v = (E x B)/B^2."""
    B0 = 1.0
    Ey = 1.0e4  # V/m
    charge, mass = k.e, k.m_p

    omega_c = k.gyrofrequency(charge, B0, mass)
    T_c = 2 * np.pi / omega_c
    dt = T_c / 60
    n_steps = 60 * 6

    t, pos, vel = pushers.boris_push(
        position=[0, 0, 0],
        velocity=[1.0e5, 0, 0],
        charge=charge, mass=mass,
        E_func=fields.uniform_E(Ey=Ey),
        B_func=fields.uniform_B(B0),
        dt=dt, n_steps=n_steps,
    )

    v_drift_theory = Ey / B0  # |E x B| / B^2, directed along +x for E_y, B_z
    v_drift_measured = (pos[-1, 0] - pos[0, 0]) / t[-1]
    print(f"[ExB]    drift theory = {v_drift_theory:.3e} m/s   "
          f"measured = {v_drift_measured:.3e} m/s")

    return t, pos, vel


def case_magnetic_mirror():
    """Particle launched into a mirror field; large pitch angle -> reflected."""
    B0 = 1.0
    charge, mass = k.e, k.m_p

    omega_c = k.gyrofrequency(charge, B0, mass)
    T_c = 2 * np.pi / omega_c
    dt = T_c / 60
    n_steps = 60 * 400

    # Mostly perpendicular velocity (large pitch angle) -> well trapped.
    v_par, v_perp = 3.0e4, 1.2e5
    t, pos, vel = pushers.boris_push(
        position=[0.0, 0.0, 0.0],
        velocity=[v_perp, 0.0, v_par],
        charge=charge, mass=mass,
        E_func=fields.zero_field(),
        B_func=fields.magnetic_mirror(B0=B0, mirror_ratio=4.0, length=2.0),
        dt=dt, n_steps=n_steps,
    )
    z = pos[:, 2]
    print(f"[mirror] z range = [{z.min():.3f}, {z.max():.3f}] m  "
          f"(stays bounded -> trapped)")

    return t, pos, vel


def main(save: bool = False):
    print("=" * 64)
    print("Experiment 01 — single-particle motion (Boris pusher)")
    print("=" * 64)

    t1, pos1, _ = case_gyro_orbit()
    t2, pos2, _ = case_ExB_drift()
    t3, pos3, _ = case_magnetic_mirror()

    fig = plt.figure(figsize=(15, 5))

    ax1 = fig.add_subplot(131)
    ax1.plot(pos1[:, 0] * 1e3, pos1[:, 1] * 1e3)
    ax1.set_aspect("equal")
    ax1.set_xlabel("x [mm]"); ax1.set_ylabel("y [mm]")
    ax1.set_title("1. Gyro-orbit (uniform B)")

    ax2 = fig.add_subplot(132)
    ax2.plot(pos2[:, 0] * 1e3, pos2[:, 1] * 1e3)
    ax2.set_xlabel("x [mm]"); ax2.set_ylabel("y [mm]")
    ax2.set_title("2. E×B drift")

    ax3 = fig.add_subplot(133)
    ax3.plot(t3 * 1e6, pos3[:, 2])
    ax3.set_xlabel("t [µs]"); ax3.set_ylabel("z [m]")
    ax3.set_title("3. Magnetic mirror (z bounces)")

    fig.tight_layout()

    if save:
        out = plotting.ensure_outputs_dir(__file__)
        path = out / "single_particle_motion.png"
        fig.savefig(path, dpi=150)
        print(f"\nSaved figure -> {path}")

    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true",
                        help="write figures to ./outputs/")
    args = parser.parse_args()
    main(save=args.save)

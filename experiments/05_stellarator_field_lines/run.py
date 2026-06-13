"""Experiment 05 (F1) — Stellarator field lines in a screw-pinch field.

The defining trick of a stellarator: get the magnetic field lines to *twist*
(rotational transform ι) purely from the shape of the field, with no driven
plasma current. Here we use the simplest field that does this — a straight
"screw pinch", B = Bz ẑ + B_θ(r) θ̂ — and trace its field lines.

Two diagnostics, both the bread-and-butter of 3-D confinement:

  1. Poincaré section — puncture the field line through a series of equally
     spaced planes and plot the dots. Nested closed curves = good flux surfaces.
  2. Rotational transform ι(r) — the average field-line twist per period. With a
     sheared B_θ(r) profile ι varies with radius, exactly as in a real device.

This is F1 on the experiment's fidelity ladder (prescribed analytic field). F2
adds real coil filaments (Biot-Savart) and magnetic islands. See PLAN.md.

Run:
    python run.py            # show plots
    python run.py --save     # also write a figure to ./outputs/
"""

from __future__ import annotations

import argparse

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay import fields, plotting
from plasmaplay.diagnostics import poincare_section, rotational_transform

# --- field definition ----------------------------------------------------
BZ = 1.0
L = 2 * np.pi          # axial period
TWIST0 = 0.3
R0 = 0.5               # profile scale length


def b_theta(r):
    """Sheared azimuthal profile -> ι decreases with radius (magnetic shear)."""
    return TWIST0 * r / (1.0 + (r / R0) ** 2)


def iota_analytic(r):
    """ι(r) = (B_θ/r)/Bz * L/(2π) for the screw pinch."""
    return (b_theta(r) / r) / BZ * L / (2 * np.pi)


def main(save: bool = False):
    print("=" * 64)
    print("Experiment 05 (F1) — stellarator field lines (screw pinch)")
    print("=" * 64)

    B = fields.screw_pinch(Bz=BZ, b_theta=b_theta)

    # --- Poincaré section: several nested surfaces -----------------------
    surface_radii = [0.15, 0.30, 0.45, 0.60, 0.75]
    fig, (axp, axi) = plt.subplots(1, 2, figsize=(12, 5.5))

    for r0 in surface_radii:
        pts = poincare_section(B, x0=[r0, 0.0, 0.0], period=L,
                               n_crossings=40, ds=0.015)
        axp.scatter(pts[:, 0], pts[:, 1], s=6, label=f"r₀={r0:.2f}")
        iota = rotational_transform(pts)
        print(f"  surface r₀={r0:.2f}:  ι_measured={iota:+.4f}  "
              f"ι_analytic={iota_analytic(r0):+.4f}")

    axp.set_aspect("equal")
    axp.set_xlabel("x [m]"); axp.set_ylabel("y [m]")
    axp.set_title("Poincaré section (nested flux surfaces)")
    axp.legend(fontsize=8, loc="upper right")

    # --- ι(r) profile: measured vs analytic ------------------------------
    radii = np.linspace(0.1, 0.85, 9)
    iota_meas = []
    for r0 in radii:
        pts = poincare_section(B, x0=[r0, 0.0, 0.0], period=L,
                               n_crossings=20, ds=0.015)
        iota_meas.append(rotational_transform(pts))

    rr = np.linspace(0.05, 0.9, 200)
    axi.plot(rr, iota_analytic(rr), "k-", lw=1.5, label="analytic")
    axi.plot(radii, iota_meas, "o", ms=6, label="traced")
    axi.set_xlabel("radius r [m]"); axi.set_ylabel("rotational transform ι")
    axi.set_title("ι(r): magnetic shear from the B_θ profile")
    axi.legend()

    fig.tight_layout()
    if save:
        out = plotting.ensure_outputs_dir(__file__)
        path = out / "screw_pinch_field_lines.png"
        fig.savefig(path, dpi=150)
        print(f"\nSaved figure -> {path}")
    plt.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--save", action="store_true",
                        help="write figure to ./outputs/")
    args = parser.parse_args()
    main(save=args.save)

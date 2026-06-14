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
from plasmaplay.diagnostics import (
    poincare_section, rotational_transform, trace_field_line,
)

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


def main(save: bool = False, n_crossings=40, profile_n=9):
    print("=" * 64)
    print("Experiment 05 (F1) — stellarator field lines (screw pinch)")
    print("=" * 64)

    B = fields.screw_pinch(Bz=BZ, b_theta=b_theta)

    # --- Poincaré section: several nested surfaces -----------------------
    surface_radii = [0.15, 0.30, 0.45, 0.60, 0.75]
    fig, (axp, axi) = plt.subplots(1, 2, figsize=(12, 5.5))

    for r0 in surface_radii:
        pts = poincare_section(B, x0=[r0, 0.0, 0.0], period=L,
                               n_crossings=n_crossings, ds=0.015)
        axp.scatter(pts[:, 0], pts[:, 1], s=6, label=f"r₀={r0:.2f}")
        iota = rotational_transform(pts)
        print(f"  surface r₀={r0:.2f}:  ι_measured={iota:+.4f}  "
              f"ι_analytic={iota_analytic(r0):+.4f}")

    axp.set_aspect("equal")
    axp.set_xlabel("x [m]"); axp.set_ylabel("y [m]")
    axp.set_title("Poincaré section (nested flux surfaces)")
    axp.legend(fontsize=8, loc="upper right")

    # --- ι(r) profile: measured vs analytic ------------------------------
    radii = np.linspace(0.1, 0.85, profile_n)
    iota_meas = []
    for r0 in radii:
        pts = poincare_section(B, x0=[r0, 0.0, 0.0], period=L,
                               n_crossings=max(5, n_crossings // 2), ds=0.015)
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


def run_stellarator(save=False):
    """A GENUINE stellarator (F2): rotational transform from a current-free 3-D helical
    vacuum field, not from plasma current (the screw-pinch above twists via current —
    the tokamak mechanism). `fields.helical_stellarator` is curl-free, so iota comes
    purely from geometry and appears at 2nd order, iota ~ eps^2."""
    print("\n--- Genuine stellarator: iota from external coils, no plasma current ---")
    L = 2 * np.pi
    eps, ll = 0.5, 2
    B = fields.helical_stellarator(eps=eps, l=ll, h=1.0)

    # (1) no net plasma current: the loop integral of B around the axis vanishes
    th = np.linspace(0, 2 * np.pi, 360, endpoint=False)
    r = 0.4
    circ = sum(np.dot(B([r * np.cos(t), r * np.sin(t), 0.0])[:2],
                      np.array([-r * np.sin(t), r * np.cos(t)]) * (2 * np.pi / 360))
               for t in th)
    print(f"  loop integral B.dl around axis = {circ:.1e}  (~0: NO net plasma current)")

    # (2) the transform grows with helical shaping (iota ~ eps^2 from geometry)
    print("  rotational transform from geometry (current-free):")
    for e in (0.3, 0.5, 0.7):
        pts = poincare_section(fields.helical_stellarator(eps=e, l=ll, h=1.0),
                               x0=[0.5, 0.0, 0.0], period=L, n_crossings=40, ds=0.01)
        print(f"    eps={e}: iota = {abs(rotational_transform(pts)):.4f}")

    # (3) E2 — the steady-state contrast: q from coils (q>1) -> no kink, no sawteeth
    from plasmaplay import sawtooth as st
    rho = np.linspace(0.0, 1.0, 65)
    q_ext = st.external_q_profile(rho)
    print(f"  E2: stellarator q from coils stays q>1 (min {q_ext.min():.2f}) -> NO "
          "q=1 surface -> no sawteeth/disruptions (inherently steady-state;")
    print("      contrast: the same burning core sawteeth ~179x in a current-driven "
          "tokamak). See gif stellarator_burn.")

    # (4) nested flux surfaces in the Poincare section
    fig = plt.figure(figsize=(11, 5.0))
    axp = fig.add_subplot(1, 2, 1)
    for r0 in (0.15, 0.30, 0.45, 0.60):
        pts = poincare_section(B, x0=[r0, 0.0, 0.0], period=L, n_crossings=80, ds=0.01)
        axp.scatter(pts[:, 0], pts[:, 1], s=5)
    axp.set_aspect("equal")
    axp.set(xlabel="x", ylabel="y", title="Poincaré: nested flux surfaces (l=2 vacuum)")
    ax3 = fig.add_subplot(1, 2, 2, projection="3d")
    for r0, c in ((0.25, "#ef476f"), (0.5, "#118ab2")):
        for j in range(4):
            a0 = np.pi * j / 2
            p = trace_field_line(B, [r0 * np.cos(a0), r0 * np.sin(a0), 0.0], ds=0.02,
                                 n_steps=900)
            phi = p[:, 2] / 3.0
            ax3.plot((3 + p[:, 0]) * np.cos(phi), (3 + p[:, 0]) * np.sin(phi), p[:, 1],
                     color=c, lw=0.6)
    ax3.set_axis_off(); ax3.set_title("twisty flux surfaces (mapped to a torus)")
    fig.tight_layout()
    if save:
        out = plotting.ensure_outputs_dir(__file__) / "stellarator_flux_surfaces.png"
        fig.savefig(out, dpi=140, bbox_inches="tight")
        print(f"  saved {out}")
    else:
        plt.show()
    plt.close(fig)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=["screwpinch", "stellarator"],
                        default="screwpinch",
                        help="screwpinch (F1, iota from current) or stellarator "
                             "(F2, iota from current-free 3-D geometry)")
    parser.add_argument("--save", action="store_true",
                        help="write figure to ./outputs/")
    args = parser.parse_args()
    if args.mode == "stellarator":
        run_stellarator(save=args.save)
    else:
        main(save=args.save)

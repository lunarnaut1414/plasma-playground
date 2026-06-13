"""Visual check of the T0 tokamak field (plasmaplay/tokamak.py).

Builds the experiment-04 Solov'ev equilibrium, wraps it as the 3-D field, and
renders four views that together *show* the three T0 validations:

  1. Poloidal cross-section  — flux surfaces (ψ contours) + poloidal B streamlines
  2. |B| over (R, Z)         — the 1/R falloff (strong on the inboard side)
  3. A 3-D field line        — winding helically around the torus on a flux surface
  4. log |∇·B|               — divergence-free to ~machine precision

Run:
    MPLBACKEND=Agg python tokamak_t0_viz.py          # headless -> outputs/tokamak_t0.png
    python tokamak_t0_viz.py                          # also shows the window
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay.diagnostics import trace_field_line
from plasmaplay.solvers import grad_shafranov_solve
from plasmaplay.tokamak import (
    divergence,
    equilibrium_field,
    to_cartesian,
    vacuum_F,
)

R0, B0 = 1.0, 2.0


def build_equilibrium(n=161):
    R = np.linspace(R0 - 0.55, R0 + 0.55, n)
    Z = np.linspace(-0.7, 0.7, n)
    RR, _ = np.meshgrid(R, Z, indexing="ij")
    psi = grad_shafranov_solve(R, Z, -(RR**2 + 1.0), boundary=0.0)
    field = equilibrium_field(R, Z, psi, vacuum_F(R0, B0))
    return R, Z, psi, field


def main(save=True):
    R, Z, psi, field = build_equilibrium()
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")

    # magnetic axis = ψ extremum
    i, j = np.unravel_index(np.argmax(psi), psi.shape)
    R_axis, Z_axis = R[i], Z[j]

    # sample the field on the φ=0 poloidal plane: at (R,0,Z), Bx=B_R, Bz=B_Z
    ng = 36
    Rs = np.linspace(R0 - 0.5, R0 + 0.5, ng)
    Zs = np.linspace(-0.62, 0.62, ng)
    BR = np.empty((ng, ng)); BZ = np.empty((ng, ng)); Bmag = np.empty((ng, ng))
    for a, r in enumerate(Rs):
        for b, z in enumerate(Zs):
            B = field(to_cartesian(r, 0.0, z))
            BR[a, b], BZ[a, b] = B[0], B[2]
            Bmag[a, b] = np.linalg.norm(B)

    # divergence map on the same plane (h spans a grid cell)
    h = 2.0 * (R[1] - R[0])
    logdiv = np.empty((ng, ng))
    for a, r in enumerate(Rs):
        for b, z in enumerate(Zs):
            logdiv[a, b] = np.log10(abs(divergence(field, to_cartesian(r, 0.0, z), h=h)) + 1e-18)

    # a 3-D field line off the axis (the "tokamak in 3-D space" money shot)
    x0 = to_cartesian(R_axis + 0.18, 0.0, Z_axis)
    line = trace_field_line(field, x0, ds=0.01, n_steps=90000)

    fig = plt.figure(figsize=(13, 10))

    # --- panel 1: flux surfaces + poloidal B streamlines -----------------
    ax1 = fig.add_subplot(2, 2, 1)
    lv = np.linspace(psi.min(), psi.max(), 16)
    ax1.contour(RR, ZZ, psi, levels=lv, colors="0.6", linewidths=0.6)
    ax1.contour(RR, ZZ, psi, levels=[0], colors="k", linewidths=1.2)
    # streamplot needs ascending grids with shape (nz, nr) -> transpose
    speed = np.hypot(BR, BZ).T
    ax1.streamplot(Rs, Zs, BR.T, BZ.T, color=speed, cmap="viridis", density=1.2,
                   linewidth=0.8, arrowsize=0.8)
    ax1.plot(R_axis, Z_axis, "rx", ms=12, mew=2.5)
    ax1.set_aspect("equal")
    ax1.set_xlabel("R [m]"); ax1.set_ylabel("Z [m]")
    ax1.set_title("Poloidal B on the flux surfaces\n(streamlines = B_pol, contours = ψ)")

    # --- panel 2: |B| showing the 1/R falloff ----------------------------
    ax2 = fig.add_subplot(2, 2, 2)
    BmagGrid = np.empty_like(RR)
    for a in range(R.size):
        for b in range(Z.size):
            BmagGrid[a, b] = np.linalg.norm(field(to_cartesian(R[a], 0.0, Z[b])))
    pc = ax2.pcolormesh(RR, ZZ, BmagGrid, shading="auto", cmap="inferno")
    ax2.contour(RR, ZZ, psi, levels=[0], colors="w", linewidths=1.0)
    ax2.plot(R_axis, Z_axis, "cx", ms=10, mew=2)
    ax2.set_aspect("equal")
    ax2.set_xlabel("R [m]"); ax2.set_ylabel("Z [m]")
    ax2.set_title("|B| — strong on the inboard side (∝ 1/R)")
    fig.colorbar(pc, ax=ax2, label="|B| [T]", shrink=0.85)

    # --- panel 3: 3-D field line winding around the torus ----------------
    ax3 = fig.add_subplot(2, 2, 3, projection="3d")
    s = np.linspace(0, 1, len(line))
    ax3.scatter(line[:, 0], line[:, 1], line[:, 2], c=s, cmap="plasma", s=0.4)
    ax3.scatter(*line[0], color="lime", s=40)
    ax3.set_xlabel("x [m]"); ax3.set_ylabel("y [m]"); ax3.set_zlabel("z [m]")
    ax3.set_title("One field line — helical winding on a flux surface")
    ax3.set_box_aspect((1, 1, 0.5))

    # --- panel 4: divergence map -----------------------------------------
    ax4 = fig.add_subplot(2, 2, 4)
    RsG, ZsG = np.meshgrid(Rs, Zs, indexing="ij")
    pc4 = ax4.pcolormesh(RsG, ZsG, logdiv, shading="auto", cmap="cividis")
    ax4.plot(R_axis, Z_axis, "rx", ms=10, mew=2)
    ax4.set_aspect("equal")
    ax4.set_xlabel("R [m]"); ax4.set_ylabel("Z [m]")
    ax4.set_title("log₁₀ |∇·B| — ≈ 0 (~1e-4, the interp residual; ≪ |B|≈1.8)")
    fig.colorbar(pc4, ax=ax4, label="log₁₀ |∇·B|", shrink=0.85)

    fig.suptitle("T0 — the 3-D tokamak field built from the experiment-04 equilibrium",
                 fontsize=14, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))

    print(f"magnetic axis: R={R_axis:.3f} Z={Z_axis:.3f}  (Shafranov shift +{R_axis-R0:.3f})")
    print(f"|B| on axis  : {np.linalg.norm(field(to_cartesian(R_axis,0,Z_axis))):.3f} T")
    print(f"max log10|divB| (interior): {logdiv[2:-2,2:-2].max():.1f}")

    if save:
        out = Path(__file__).resolve().parent / "outputs"
        out.mkdir(exist_ok=True)
        path = out / "tokamak_t0.png"
        fig.savefig(path, dpi=130)
        print(f"saved -> {path}")
    plt.show()


if __name__ == "__main__":
    main()

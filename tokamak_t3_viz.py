"""T3 — breaking axisymmetry: magnetic islands and the road to stochasticity.

Rung T3 of `docs/3D_TOKAMAK_GUIDE.md` — the first genuinely **3-D** field
structure. A tokamak equilibrium is axisymmetric (2-D); add a non-axisymmetric
helical perturbation δψ = δ·env(r)·cos(mθ − nφ) (built through the flux, so the
total field stays exactly divergence-free) and the nested flux surfaces **tear**:

  * at a rational surface **q = m/n** the perturbation is resonant and opens an
    **m-island chain**; its width scales like **√(amplitude)**;
  * away from rationals the KAM surfaces survive (just rippled);
  * when two island chains grow until their separatrices **overlap** (Chirikov),
    the field lines between them become **chaotic** — a stochastic sea, the loss
    of good surfaces that limits tokamak confinement.

The equilibrium is a sheared circular one (q₀=1.2) chosen so the q=3/2 and q=2
surfaces sit at mid-radius; we drive the 3/2 and 2/1 modes and watch them go from
isolated islands to a stochastic layer as the amplitude is raised.

Run:
    MPLBACKEND=Agg python tokamak_t3_viz.py    # headless -> outputs/tokamak_t3.png
    python tokamak_t3_viz.py                    # also show the window
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from plasmaplay.tokamak import (
    equilibrium_field,
    helical_perturbation,
    safety_factor,
    superpose,
    toroidal_poincare,
    vacuum_F,
)

Rc, B0c, a, q0 = 5.0, 1.0, 1.0, 1.2
AXIS = (Rc, 0.0)


def base_field(n=261):
    R = np.linspace(Rc - 1.15, Rc + 1.15, n)
    Z = np.linspace(-1.15, 1.15, n)
    RR, ZZ = np.meshgrid(R, Z, indexing="ij")
    psi = (B0c * a**2 / (2.0 * q0)) * np.log1p(((RR - Rc) ** 2 + ZZ ** 2) / a**2)
    return equilibrium_field(R, Z, psi, vacuum_F(Rc, B0c))


def q_profile(field, rs):
    return np.array([abs(safety_factor(field, (Rc + r, 0.0), AXIS,
                                       n_poloidal=5, ds=0.07)) for r in rs])


def gaussian_env(r0, sigma=0.16):
    return lambda r: np.exp(-((r - r0) ** 2) / (2 * sigma**2))


def poincare_multi(field, starts, npunc=120, ds=0.07):
    """Trace a Poincaré section for each explicit (R, Z) start point."""
    return [toroidal_poincare(field, st, n_punctures=npunc, ds=ds) for st in starts]


# background surfaces: launched on the outboard midplane (θ = 0)
BG = [(Rc + s, 0.0) for s in np.linspace(0.2, 0.98, 9)]


def island_seeds(*r_res):
    """Seeds at the X-point (θ = π/2, i.e. R = Rc, Z = r) of each chain — these
    trace the separatrix and so reveal the full island lobes (a midplane seed
    would sit on an O-point and trace only a small contour)."""
    return [(Rc, r) for r in r_res]


def _plot_section(ax, sections, title, n_bg=len(BG)):
    for k, pc in enumerate(sections):
        if k < n_bg:
            ax.plot(pc[:, 0], pc[:, 1], ".", ms=0.8, color="0.6")
        else:
            ax.plot(pc[:, 0], pc[:, 1], ".", ms=1.4)     # island lines, coloured
    ax.plot(Rc, 0, "k+", ms=8)
    ax.set_aspect("equal")
    ax.set_xlabel("R [m]"); ax.set_ylabel("Z [m]")
    ax.set_title(title)


def main(save=True):
    field0 = base_field()

    # locate the resonant surfaces
    rs = np.linspace(0.15, 1.0, 12)
    qs = q_profile(field0, rs)
    r_32 = float(np.interp(1.5, qs, rs))
    r_21 = float(np.interp(2.0, qs, rs))
    print(f"resonances: q=3/2 at r≈{r_32:.3f}, q=2 at r≈{r_21:.3f}")

    env_32 = gaussian_env(r_32)
    env_21 = gaussian_env(r_21)

    fig = plt.figure(figsize=(13, 10))

    # --- panel A: unperturbed nested surfaces + q-profile ----------------
    ax = fig.add_subplot(2, 2, 1)
    _plot_section(ax, poincare_multi(field0, BG), "Unperturbed: nested flux surfaces")
    axq = ax.inset_axes([0.62, 0.08, 0.34, 0.32])
    rfine = np.linspace(0.15, 1.0, 30)
    axq.plot(q_profile(field0, rfine), rfine, "k-", lw=1)
    for qv in (1.5, 2.0):
        axq.axvline(qv, color="C3", ls=":", lw=0.8)
    axq.set_xlabel("q", fontsize=6); axq.set_ylabel("r", fontsize=6)
    axq.tick_params(labelsize=5)

    # --- panel B: single 2/1 mode -> one clean island chain --------------
    f_single = superpose(field0, helical_perturbation(6e-4, 2, 1, AXIS, envelope=env_21))
    ax = fig.add_subplot(2, 2, 2)
    _plot_section(ax, poincare_multi(f_single, BG + island_seeds(r_21), npunc=200),
                  f"2/1 mode (small δ): island chain at q=2 (r≈{r_21:.2f})")

    # --- panel C: 3/2 + 2/1, moderate -> two island chains ---------------
    f_two = superpose(field0,
                      helical_perturbation(6e-4, 3, 2, AXIS, envelope=env_32),
                      helical_perturbation(6e-4, 2, 1, AXIS, envelope=env_21))
    ax = fig.add_subplot(2, 2, 3)
    _plot_section(ax, poincare_multi(f_two, BG + island_seeds(r_32, r_21), npunc=200),
                  "3/2 + 2/1 (moderate δ): two island chains")

    # --- panel D: same modes, large amplitude -> stochastic sea ----------
    f_chaos = superpose(field0,
                        helical_perturbation(4e-3, 3, 2, AXIS, envelope=env_32),
                        helical_perturbation(4e-3, 2, 1, AXIS, envelope=env_21))
    ax = fig.add_subplot(2, 2, 4)
    _plot_section(ax, poincare_multi(f_chaos, BG + island_seeds(r_32, r_21), npunc=260),
                  "3/2 + 2/1 (large δ): overlap → stochastic sea")

    fig.suptitle("T3 — magnetic islands at q=m/n, and the Chirikov road to stochasticity",
                 fontsize=14, weight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.96))

    if save:
        out = Path(__file__).resolve().parent / "outputs"
        out.mkdir(exist_ok=True)
        path = out / "tokamak_t3.png"
        fig.savefig(path, dpi=130)
        print(f"saved -> {path}")
    plt.show()


if __name__ == "__main__":
    main()

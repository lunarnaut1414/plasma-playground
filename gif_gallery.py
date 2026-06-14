"""gif_gallery.py — regenerate every animated showcase .gif.

One function per animation, registered in GALLERY. Each prints its validation
number(s) and writes a .gif (and, for operation-mode showcases, a PNG still) to
`outputs/`. The night's autonomous build (see NIGHT.md) adds one entry per rung.
(Static-figure montage lives in the separate `gallery.py`.)

Usage:
    python gif_gallery.py                 # list available gifs
    python gif_gallery.py <name> [...]    # regenerate the named gif(s)
    python gif_gallery.py all             # regenerate everything
"""

from __future__ import annotations

import sys

import numpy as np

from plasmaplay import animate as anim, transport as tr

OUT = "outputs"


def smoke_diffusion():
    """G1 foundation smoke: a 1-D diffusing Gaussian. Validates the gif pipeline
    against the analytic diffusion solution (mass conservation + peak-decay law)."""
    x = np.linspace(-8, 8, 161)
    times = np.linspace(0, 6, 90)
    D, s0 = 0.7, 1.0

    def gauss(t):
        s2 = s0**2 + 2 * D * t
        return 1.0 / np.sqrt(2 * np.pi * s2) * np.exp(-x**2 / (2 * s2))

    frames = anim.make_frames(gauss, times)
    mass = np.trapezoid(frames, x, axis=1)
    peak = frames.max(axis=1)
    exp_peak = peak[0] / np.sqrt(1 + 2 * D * times / s0**2)
    print(f"  [smoke_diffusion] mass drift = {np.abs(mass / mass[0] - 1).max():.2e} "
          f"(conserved); peak-law err = {np.abs(peak / exp_peak - 1).max():.2e}")
    out = anim.animate_profiles(x, frames, times, path=f"{OUT}/_smoke_diffusion.gif",
                                ylabel="c(x,t)", title="G1 smoke: 1-D diffusion",
                                fps=20, dpi=90)
    print(f"  wrote {out}")


def burn_0d_ignition():
    """A1 (F1): 0-D burn with He ash + beta-limit. Phase-space (n,T) ignition track
    igniting onto the beta-limited operating point, the point colored by ash
    fraction as ash builds up. Validates the steady ash balance n_He = tau_he*R."""
    beta_lim = tr.troyon_limit(3.0, 7.0, 1.0, 5.3)   # ~3.96%

    def paux(t):
        return 5.0e5 if t < 5.0 else 0.0

    r = tr.burn_0d_ash(1.0e20, 5.0, tau_E=3.0, p_aux=paux, B=5.3, tau_p=6.0,
                       tau_he=10.0, fuel_rate=1.0e20 / 6.0, beta_limit=beta_lim,
                       t_end=40.0)
    R = tr.reaction_rate_dt(r["n_DT"][-1], r["T"][-1])
    print(f"  [burn_0d_ignition] steady T = {r['T'][-1]:.1f} keV, "
          f"beta = {r['beta'][-1]*100:.2f}% (limit {beta_lim*100:.2f}%), "
          f"f_He = {r['f_He'][-1]*100:.1f}%, ash balance = "
          f"{r['n_He'][-1]/(10.0*R):.3f}")
    # subsample to a sane frame count
    s = slice(0, None, max(1, r["t"].size // 100))
    out = anim.animate_phase_track(
        r["n_e"][s], r["T"][s], r["t"][s], path=f"{OUT}/burn_0d_ignition.gif",
        color=r["f_He"][s] * 100, xlabel=r"$n_e$ [m$^{-3}$]", ylabel="T [keV]",
        clabel="ash fraction [%]", title="0-D ignition onto the burning point",
        ylim=(0, 20), fps=20, dpi=90)
    print(f"  wrote {out}")


GALLERY = {
    "smoke_diffusion": smoke_diffusion,
    "burn_0d_ignition": burn_0d_ignition,
}


def main(argv):
    if not argv:
        print("Available gifs:")
        for name in GALLERY:
            print(f"  {name}")
        print("\nUsage: python gif_gallery.py <name>[ ...] | all")
        return
    names = list(GALLERY) if argv == ["all"] else argv
    for name in names:
        if name not in GALLERY:
            print(f"  ! unknown gif '{name}' (skipped)")
            continue
        print(f"[{name}]")
        GALLERY[name]()


if __name__ == "__main__":
    main(sys.argv[1:])

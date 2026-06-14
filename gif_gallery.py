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

from plasmaplay import animate as anim

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


GALLERY = {
    "smoke_diffusion": smoke_diffusion,
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

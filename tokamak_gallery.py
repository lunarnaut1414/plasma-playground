"""Stitch the five 3-D-tokamak-ladder figures (T0–T4) into one montage.

    python tokamak_gallery.py            # build outputs/tokamak_gallery.png and show
    python tokamak_gallery.py --no-show  # write file only (headless)

Reads the per-rung PNGs already written by the tokamak_t{0..4}_viz.py scripts
(it does *not* re-render them — those traces take minutes). Regenerate a missing
panel with e.g. `MPLBACKEND=Agg python tokamak_t2_viz.py`.
"""

from __future__ import annotations

import argparse
import pathlib

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent
OUT = ROOT / "outputs"

PANELS = [
    ("T0 · 3-D equilibrium field (∇·B≈0, |B|∝1/R)", "tokamak_t0.png"),
    ("T1 · field-line topology: Poincaré & q-profile", "tokamak_t1.png"),
    ("T2 · banana orbits, trapping & the μ invariant", "tokamak_t2.png"),
    ("T3 · magnetic islands & stochasticity (real 3-D)", "tokamak_t3.png"),
    ("T4 · linear resistive tearing mode (γ ∝ S^−3/5)", "tokamak_t4.png"),
]


def build(show=True):
    fig, axes = plt.subplots(3, 2, figsize=(16, 17))
    fig.suptitle("plasma-playground — the 3-D tokamak ladder (T0 → T4)", fontsize=17, y=0.997)
    missing = []
    for ax, (title, fname) in zip(axes.flat, PANELS):
        path = OUT / fname
        if path.exists():
            ax.imshow(mpimg.imread(path))
            ax.set_title(title, fontsize=11)
        else:
            ax.text(0.5, 0.5, f"missing {fname}\n(run its tokamak_t*_viz.py)",
                    ha="center", va="center", fontsize=10, color="C3")
            missing.append(fname)
        ax.axis("off")
    axes.flat[-1].axis("off")     # 5 panels in a 3×2 grid; last cell blank

    fig.tight_layout(rect=[0, 0, 1, 0.985])
    OUT.mkdir(exist_ok=True)
    path = OUT / "tokamak_gallery.png"
    fig.savefig(path, dpi=110)
    print(f"Tokamak gallery written -> {path}")
    if missing:
        print(f"  (placeholders for not-yet-rendered: {', '.join(missing)})")
    if show:
        plt.show()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--no-show", action="store_true", help="write file only")
    args = p.parse_args()
    build(show=not args.no_show)

"""Visual progress check — run every experiment and assemble one gallery image.

    python gallery.py            # build outputs/gallery.png and show it
    python gallery.py --no-show  # just write the file (headless)

Each experiment regenerates its own figure (into its experiments/*/outputs/), then
we stitch them into a single montage at outputs/gallery.png so you can see
everything we've built so far at a glance.
"""

from __future__ import annotations

import argparse
import importlib.util
import pathlib

import matplotlib.image as mpimg
import matplotlib.pyplot as plt

ROOT = pathlib.Path(__file__).resolve().parent
EXPERIMENTS = ROOT / "experiments"


def _load(exp_dir):
    path = EXPERIMENTS / exp_dir / "run.py"
    spec = importlib.util.spec_from_file_location(f"run_{exp_dir}", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def build(show=True):
    # Suppress each experiment's blocking plt.show(); we render to files instead.
    real_show = plt.show
    plt.show = lambda *a, **k: None

    print("Running experiments (this regenerates each figure)...")
    _load("01_single_particle_motion").main(save=True)
    _load("04_tokamak_equilibrium").main(save=True)
    _load("05_stellarator_field_lines").main(save=True)
    _load("06_ideal_mhd").run_briowu(save=True)
    exp03 = _load("03_pic_1d")
    exp03.run_cold(save=True)
    exp03.run_landau(save=True)
    exp03.run_twostream(save=True)
    plt.close("all")
    plt.show = real_show

    panels = [
        ("01 · single-particle motion",
         "experiments/01_single_particle_motion/outputs/single_particle_motion.png"),
        ("04 · tokamak equilibrium (Grad-Shafranov)",
         "experiments/04_tokamak_equilibrium/outputs/grad_shafranov_equilibrium.png"),
        ("05 · stellarator field lines",
         "experiments/05_stellarator_field_lines/outputs/screw_pinch_field_lines.png"),
        ("06 · ideal MHD (Brio-Wu shock tube)",
         "experiments/06_ideal_mhd/outputs/brio_wu.png"),
        ("03 · Landau damping",
         "experiments/03_pic_1d/outputs/landau_damping.png"),
        ("03 · two-stream instability",
         "experiments/03_pic_1d/outputs/two_stream.png"),
    ]

    fig, axes = plt.subplots(3, 2, figsize=(16, 14))
    fig.suptitle("plasma-playground — progress gallery", fontsize=16, y=0.995)
    for ax, (title, rel) in zip(axes.flat, panels):
        ax.imshow(mpimg.imread(ROOT / rel))
        ax.set_title(title, fontsize=11)
        ax.axis("off")
    fig.tight_layout(rect=[0, 0, 1, 0.98])

    out = ROOT / "outputs"
    out.mkdir(exist_ok=True)
    path = out / "gallery.png"
    fig.savefig(path, dpi=120)
    print(f"\nGallery written -> {path}")
    if show:
        plt.show()


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--no-show", action="store_true", help="write file only")
    args = p.parse_args()
    build(show=not args.no_show)

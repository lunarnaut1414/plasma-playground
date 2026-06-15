"""Build web-friendly 'money shot' gifs for the showcase (assets/README.md).

The full-resolution gifs in outputs/ are large (the 3-D burns are ~10 MB at 200
frames). This downsamples the curated showcase set into assets/ — fewer frames + an
adaptive palette + optional scaling — so they embed in a README and load fast, while
the originals stay in outputs/ (gitignored). Run from the repo root after the gifs
exist: `python scripts/build_assets.py`.
"""
import os

from PIL import Image

OUT, ASSETS = "outputs", "assets"

# (source, frame_step, scale, fps, colors) — one entry per showcase gif.
# No dithering (it destroys GIF run-length compression); fewer frames + a smaller
# adaptive palette + scaling are the size levers.
SHOTS = [
    ("tokamak_3d_discharge", 2, 0.7, 18, 144),
    ("stellarator_3d_burn", 2, 0.7, 18, 144),
    ("tokamak_discharge_full", 2, 0.78, 14, 180),
    ("stellarator_flux_surfaces", 1, 0.78, 16, 144),
    ("tearing_island_saturation", 2, 0.78, 14, 180),
    ("operating_modes", 2, 0.85, 18, 200),
]


def optimize(name, frame_step, scale, fps, colors):
    src = f"{OUT}/{name}.gif"
    dst = f"{ASSETS}/{name}.gif"
    im = Image.open(src)
    frames = []
    for i in range(0, im.n_frames, frame_step):
        im.seek(i)
        f = im.convert("RGB")
        if scale != 1.0:
            f = f.resize((round(f.width * scale), round(f.height * scale)), Image.LANCZOS)
        frames.append(f.quantize(colors=colors, method=Image.MEDIANCUT, dither=Image.NONE))
    frames[0].save(dst, save_all=True, append_images=frames[1:], loop=0,
                   duration=int(1000 / fps), optimize=True, disposal=2)
    mb = os.path.getsize(dst) / 1048576
    print(f"  {name:30s} {len(frames):3d} frames  {frames[0].size}  ->  {mb:5.2f} MB")


def main():
    os.makedirs(ASSETS, exist_ok=True)
    for name, step, scale, fps, colors in SHOTS:
        if os.path.exists(f"{OUT}/{name}.gif"):
            optimize(name, step, scale, fps, colors)
        else:
            print(f"  SKIP {name} (no outputs/{name}.gif)")


if __name__ == "__main__":
    main()

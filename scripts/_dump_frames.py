"""Sample representative frames from a gif to PNG so they can be eyeballed.

Usage: python scripts/_dump_frames.py outputs/<name>.gif
Writes outputs/_review/<name>_NNN.png for a spread of frames (weighted to the ends).
Part of the REVIEW.md gif-review loop. outputs/ is gitignored.
"""
import os
import sys

from PIL import Image


def main(path):
    im = Image.open(path)
    n = im.n_frames
    name = os.path.splitext(os.path.basename(path))[0]
    out = "outputs/_review"
    os.makedirs(out, exist_ok=True)
    picks = sorted({0, n // 8, n // 4, n // 2, (3 * n) // 4, (7 * n) // 8, n - 1})
    for f in picks:
        im.seek(f)                                   # seek+copy: do NOT list() the iterator
        im.convert("RGB").save(f"{out}/{name}_{f:03d}.png")
    print(f"{name}: n_frames={n} per_frame={im.info.get('duration', '?')}ms saved={picks}")


if __name__ == "__main__":
    main(sys.argv[1])

# G1 — Visualization foundation (`plasmaplay/animate.py`)

The gif pipeline every later rung depends on (NIGHT.md Track G). Built first so the
burn/MHD rungs have a validated way to turn a time-series of fields into a `.gif`.

## What was built

- `make_frames(field_fn, times)` — pure sampler: stacks `field_fn(t)` over a time
  grid into `(n_t, ...)`. Separating frame-data construction from rendering is the
  whole design point — the *physics* is unit-tested without ever writing gif bytes.
- `torus_surface(R0, a)` — parametric torus `(X,Y,Z)`; every point satisfies the
  implicit torus equation (the unit-test invariant).
- `animate_profiles` — radial line profiles over time (single or multi-series, e.g.
  Te & Ti) → gif.
- `animate_cross_section` — revolves a radial profile into a poloidal disk heatmap
  → gif (the "watch it burn" view).
- `animate_torus_3d` — rotating 3-D torus colored per-frame (showcase view; a scalar
  edge-color stand-in until a later rung colors by a poloidally-resolved field).

All write via matplotlib `PillowWriter` — pure Python, no external binary.

## Validation (falsifiable, in `tests/test_animate.py`)

The reference is the exact 1-D diffusion of a Gaussian, `σ²(t)=σ₀²+2Dt`:

- **Mass conservation:** `∫c dx` invariant across frames → drift `9.08e-03` (wide grid).
- **Peak-decay law:** peak falls as `1/√(1+2Dt/σ₀²)` → err `2.22e-16`.
- **Torus geometry:** `max|(√(X²+Y²)−R0)² + Z² − a²| < 1e-9`.
- Each `animate_*` writes a non-empty `.gif` (pipeline smoke).

7 tests, all green (full suite 147 passed, ruff clean).

## Deliverable

`outputs/_smoke_diffusion.gif` — the diffusing Gaussian, proof the pipeline works.
Regenerate with `python gif_gallery.py smoke_diffusion` (prints both validation
numbers). `gif_gallery.py` is the registry the night's later rungs append to.

## Gotcha

`matplotlib.use("Agg", force=False)` respects an already-set backend; run headless
work with `MPLBACKEND=Agg`. Keep gifs modest (~60–120 frames, dpi≈90) — the smoke
gif is 460K at 90 frames / dpi 90.

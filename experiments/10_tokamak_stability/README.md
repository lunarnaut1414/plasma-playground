# 10 — Tokamak MHD stability (the straight-tokamak / cylinder limit)

The **fast** half of a tokamak discharge: the ideal/resistive MHD instabilities that
move the plasma fluid on the Alfvén (µs) timescale — the counterpart to the slow
*transport* arc of [experiment 09](../09_burning_plasma/). We start in the cleanest
geometry that still has a safety factor: a **periodic cylinder** (a torus cut and
straightened — the "straight tokamak").

## What it shows

A current profile `J_z(r) = J0 (1 − (r/a)²)^ν` gives a poloidal field `B_θ(r)` and a
safety factor `q(r) = r B_z / (R B_θ)`, rising from `q(0)=q0` to `q(a)=(ν+1)q0`. On it:

- **m=1/n=1 internal kink** — the rigid sideways shift of the core inside the q=1
  surface, unstable **exactly when q(0) < 1**. The *sawtooth trigger*.
- **m/n tearing modes** — reconnection at the `q = m/n` rational surface, governed by
  the outer **Newcomb equation** and its stability index **Δ′** (sign = stability),
  with the resistive growth rate following the FKR `γ ∝ S^(−3/5)` layer law.

## Run

```bash
python run.py            # prints the kink criterion, Delta' signs, and the FKR scaling
python run.py --save     # also writes outputs/kink_eigenmode.png
MPLBACKEND=Agg python ../../gif_gallery.py kink_eigenmode   # the eigenmode gif
```

## Concepts

Safety factor & rational surfaces · Newcomb equation · tearing stability index Δ′ ·
the m=1 internal kink and the q(0)<1 sawtooth trigger · Furth–Killeen–Rosenbluth
resistive-layer growth `γ ∝ Δ′^{4/5} S^{−3/5}`.

Rungs **B1** (cylinder linear stability) and **B2** (nonlinear reduced-MHD island +
Rutherford saturation, `run.py --island`) are implemented, reusing the slab-tearing
layer physics of [`plasmaplay/tearing.py`](../../plasmaplay/tearing.py) (T4) on a real
`q(r)` and Harris sheet. Next: the sawtooth cycle (B3), then the staged coupling of a
sawtooth/tearing event into the experiment-09 burn (Track C). See `PLAN.md`.

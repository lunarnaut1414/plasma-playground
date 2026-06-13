# 08 — Plasma waves & dispersion

Every plasma wave has a fingerprint: its dispersion relation ω(k). This
experiment is the **unifier** — it computes the Langmuir dispersion from kinetic
theory and then *measures the same curve out of a PIC simulation*, tying the
kinetic (03) and analytic pictures together.

**Status:** F2 + F3 implemented (`run.py`). See [`PLAN.md`](PLAN.md).

## What it shows

| Case (`--case`) | Physics |
|------|---------|
| `kinetic` (F2) | solve the kinetic dispersion relation via the plasma dispersion function Z(ζ) → the real frequency **and** the Landau damping rate, vs the cold (ω_pe) and Bohm-Gross (fluid) limits |
| `pic` (F3) | run a warm PIC plasma, FFT the field history E(x, t) into the **ω–k plane**, and watch the simulated dispersion ridge land on the analytic curve |

The kinetic solver reproduces the textbook Landau benchmarks exactly
(k λ_D = 0.5 → ω_r = 1.4157 ω_pe, γ = −0.1534 ω_pe) — and that γ is the same
theory line experiment 03's PIC Landau damping (V6) was measured against.

## Run it

```bash
python run.py --case kinetic [--save]
python run.py --case pic
python run.py --case all --save
```

## Kernels exercised (validated in `tests/test_dispersion.py`)

- `plasmaplay.dispersion` — `plasma_dispersion_function` (Z, cross-checked vs
  PlasmaPy), `langmuir_dispersion` (complex ω) — **V9**
- `plasmaplay.diagnostics.omega_k_spectrum` — 2-D FFT of E(x, t) — **V15**

## Concepts / keywords

- Dispersion relation ω(k); cold → warm-fluid (Bohm-Gross) → kinetic
- The plasma dispersion function Z(ζ); Landau damping from complex ω
- ω–k spectral analysis of simulation data

## How it closes the loop

Experiment 03 generates the PIC field data; experiment 08 measures its dispersion
and matches it to the kinetic theory derived here. High-k modes wash into noise
precisely because they are strongly Landau-damped — the F2 and F3 results are two
views of the same physics. Next (F4): the full hot, magnetized dispersion
(Bernstein modes, cyclotron harmonics).

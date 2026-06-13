# 03 — 1-D electrostatic Particle-in-Cell

The first **self-consistent, many-particle** experiment: the particles create the
electric field that pushes them. That feedback loop is what makes a plasma a
plasma — and it produces collective behavior (oscillations, damping,
instabilities) with no single-particle analogue.

**Status:** F1 + F2 implemented (`run.py`). See [`PLAN.md`](PLAN.md).

## The PIC cycle

```
deposit charge to grid  ->  solve Poisson for E  ->  gather E to particles  ->  push
   (cic_deposit)            (solvers.solve_efield_1d)   (cic_interpolate)      (leapfrog)
```

All in normalized units (ε0 = 1) with the macroparticle charge/mass chosen so
**ω_pe = 1**, making every measured frequency/rate directly comparable to theory.

## What it shows

| Case (`--case`) | Physics | Result | Test |
|------|---------|--------|------|
| `cold` (F1) | perturbed cold plasma | oscillates at **ω_pe** (measured 1.005) | V5 |
| `landau` (F2) | warm Langmuir wave | **collisionless Landau damping**, γ≈0.15 ω_pe | V6 |
| `twostream` (F2) | two counter-streaming beams | **two-stream instability** + phase-space hole | V7 |

## Run it

```bash
python run.py --case cold       [--save]
python run.py --case landau
python run.py --case twostream
python run.py --case all --save
```

## Kernels exercised (all validated in `tests/`)

- `plasmaplay.solvers.solve_efield_1d` / `solve_poisson_1d` — V4
- `plasmaplay.pic` — `cic_deposit`, `cic_interpolate`, loaders, `ElectrostaticPIC1D`
- `plasmaplay.diagnostics.dominant_frequency` — spectral peak of a time series

## Things to play with (the physics knobs)

- **Resolve the Debye length.** A truly cold plasma (`v_th = 0`) has λ_D = 0, so
  the grid is always under-resolved and a numerical *grid-heating instability*
  grows. The cold case adds a small `v_th` so λ_D > dx/π — try setting it back to
  0 and watch the energy blow up. (This is why `dx < ~π λ_D` is a PIC rule.)
- **Landau damping** vanishes as k λ_D → 0 and strengthens as it grows; the
  benchmark here is k λ_D = 0.5 (γ ≈ 0.1533 ω_pe).
- **Two-stream** is unstable only for k·v_beam < ω_pe; push v_beam up past ω_pe
  and the instability switches off.

## Concepts / keywords

- Plasma frequency ω_pe, Debye length λ_D, Bohm–Gross dispersion
- Cloud-in-cell weighting, leapfrog, the CFL/Debye resolution conditions
- Landau damping (collisionless!), two-stream instability, BGK phase-space holes

## Next rung (F3 / F4)

Electromagnetic (1D3V) PIC, or 2-D, or benchmarking against a production code
(Smilei / WarpX). See [`PLAN.md`](PLAN.md). The field histories here also feed
experiment 08, where we *measure* the dispersion relation from the ω–k spectrum.

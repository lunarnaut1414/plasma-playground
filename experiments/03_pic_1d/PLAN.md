# 03 — 1-D Particle-in-Cell — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** not started. The first *self-consistent, many-particle* experiment.

## The question

What happens when the particles make the fields that push them? Drop the
"prescribed field" crutch of experiments 01–02 and let charge density drive a
self-consistent electric field. This is where collective behavior — oscillations,
damping, instabilities — emerges, and where it starts to feel like real plasma
simulation.

## Why it matters (real devices)

Particle-in-Cell (PIC) is *the* workhorse kinetic method — used for everything
from laser-plasma accelerators to thruster plumes to fusion edge physics.
Landau damping and the two-stream instability you'll reproduce here are textbook
phenomena with no fluid analogue; they only exist because plasma is made of
particles with a velocity distribution.

## Prerequisites

Experiment 01 (the leapfrog/particle push). Comfort with FFTs and 1-D grids.
Concept of a distribution function f(x, v).

## Fidelity ladder

### F0 — Analytic / sanity
- **Models:** plasma frequency ω_pe, Debye length λ_D, the Landau damping rate of a Langmuir wave, the two-stream growth rate, Bohm–Gross dispersion ω² = ω_pe² + 3k²v_th².
- **Assumes:** linear theory, Maxwellian (or two-beam) distribution.
- **Method & tools:** arithmetic; PlasmaPy `plasmapy.dispersion` for cross-checks.
- **You'll learn:** the numbers your simulation must hit (oscillation frequency, damping rate, growth rate).
- **Validation:** reference; the F1–F2 sims are judged against these.
- **Compute:** instant.

### F1 — Cold electrostatic PIC (the machine)
- **Models:** 1D1V electrostatic PIC loop — weight particles to grid → solve Poisson (FFT) → interpolate field → push (leapfrog). Cold plasma slab oscillating at ω_pe.
- **Assumes:** 1 spatial + 1 velocity dimension, electrostatic, cold (no thermal spread), immobile ions.
- **Method & tools:** NumPy + FFT Poisson solve; CIC (cloud-in-cell) weighting; this is the reusable PIC core.
- **You'll learn:** the full PIC cycle end-to-end; that a perturbed cold slab oscillates at exactly ω_pe.
- **Validation:** measured oscillation frequency = ω_pe to a few %; total energy bounded; require dx < λ_D to avoid grid heating.
- **Compute:** seconds; 10⁴–10⁵ particles, single core (numba optional).

### F2 — Warm plasma: Landau damping & two-stream
- **Models:** Maxwellian electrons → collisionless Landau damping of a Langmuir wave; two counter-streaming beams → two-stream instability that grows then saturates by particle trapping.
- **Assumes:** still 1D1V electrostatic; now a real velocity distribution.
- **Method & tools:** the F1 core + Maxwellian/two-beam loaders; diagnostics on mode energy vs. time and the v–x phase space.
- **You'll learn:** Landau damping with *no collisions* (the famous result); nonlinear saturation and phase-space "vortices."
- **Validation:** measured damping rate matches the F0 Landau rate; two-stream growth rate matches theory; phase-space holes appear at saturation.
- **Compute:** seconds–minutes; numba recommended for nicer-resolved runs.

### F3 — Electromagnetic / nonlinear extensions
- **Models:** go to 1D3V electromagnetic (add B, currents, a proper field solve) — or stay electrostatic and push deep into the nonlinear/BGK-mode regime. Careful energy and momentum diagnostics.
- **Assumes:** 1 spatial dimension but 3 velocity/field components; explicit time stepping (CFL-limited).
- **Method & tools:** NumPy + `numba`; a current-deposition + EM field update (e.g. FDTD-style); stability/heating diagnostics.
- **You'll learn:** numerical heating, the CFL condition, why EM PIC is harder; electromagnetic instabilities (e.g. Weibel) if you add the right setup.
- **Validation:** energy conservation within scheme tolerance; reproduce a known EM instability growth rate.
- **Compute:** minutes; numba parallel helps.

### F4 — 2-D, or run a real code
- **Models:** a 2-D electromagnetic PIC of your own (numba-parallel), **or** set up the same physics in an established code — Smilei or WarpX — and compare.
- **Assumes:** depends on path; this rung is about scale and validation against a reference implementation.
- **Method & tools:** `numba`/`jax` for your own 2-D; or Smilei/WarpX (both Python-scriptable) as the reference.
- **You'll learn:** what production PIC codes do that toy codes don't (load balancing, current conservation schemes, absorbing boundaries); how to read and trust a real code.
- **Validation:** your 2-D result matches the reference code on a shared benchmark (e.g. Weibel or two-stream growth rate).
- **Compute:** minutes–hours; the 96 GB unified memory is the enabler for larger 2-D grids on CPU.

## Diagnostics you'll reuse
Phase-space (x–v) scatter/heatmaps, mode-energy-vs-time (semilog for growth/damping rates), field energy partition, ω–k spectra (links to experiment 08).

## Key references
- Birdsall & Langdon, *Plasma Physics via Computer Simulation* — the PIC bible.
- Hockney & Eastwood, *Computer Simulation Using Particles*.
- Smilei (smileipic.github.io) and WarpX (ecp-warpx.github.io) docs.

## Stretch goals
- Animate phase-space hole formation in the two-stream run.
- Reproduce the Landau damping rate across a sweep of k and overlay on the F0 curve.
- Feed F2 field data into experiment 08 to *measure* the dispersion relation.

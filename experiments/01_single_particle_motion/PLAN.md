# 01 — Single-particle motion — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Current status:** F1 implemented (`run.py`). F0 helpers exist in `plasmaplay.constants`.

## The question

How does a single charged particle move in prescribed electric and magnetic
fields? This is the irreducible atom of plasma physics — every other experiment
is, at bottom, many copies of this with fields the particles create themselves.

## Why it matters (real devices)

Magnetic confinement *is* the statement "particles are stuck to field lines and
only leak slowly." Gyro-motion, drifts, and mirror trapping are why a tokamak or
stellarator can hold a 100-million-kelvin plasma off the wall. Get the intuition
here and the device experiments stop being mysterious.

## Prerequisites

None. This is the entry point. Helpful: vector cross products, Newton's law.

## Fidelity ladder

### F0 — Analytic / sanity
- **Models:** gyrofrequency ω_c = qB/m, Larmor radius r_L = mv⊥/qB, E×B drift speed E/B, magnetic moment μ = mv⊥²/2B.
- **Assumes:** uniform fields, non-relativistic, single particle.
- **Method & tools:** arithmetic; `plasmaplay.constants` helpers.
- **You'll learn:** the natural scales (how fast, how big) that every later rung must reproduce.
- **Validation:** this *is* the reference. Cross-check against PlasmaPy's `plasmapy.formulary`.
- **Compute:** instant.

### F1 — Minimal numerical  ✅ implemented
- **Models:** full Lorentz-force orbit in uniform B, crossed E×B, and a paraxial magnetic mirror.
- **Assumes:** prescribed (fixed) fields — the particle does not affect them. Non-relativistic.
- **Method & tools:** Boris pusher (`plasmaplay.pushers`), pure NumPy.
- **You'll learn:** the gyro-orbit, the charge-independent E×B drift, mirror reflection and the loss cone; *why* the Boris algorithm (energy-conserving, time-reversible) beats RK4 for this.
- **Validation:** relative energy drift ~1e-15 over many orbits; measured E×B drift matches E/B; mirror keeps z bounded. (All printed by `run.py`.)
- **Compute:** seconds, single core.

### F2 — Guiding-center & higher-order drifts
- **Models:** add grad-B and curvature drifts in *non-uniform* fields; compare full-orbit to the guiding-center approximation; relativistic Boris (Vay or Higuera–Cary) for energetic particles.
- **Assumes:** still prescribed fields, but now spatially varying; slow variation over a gyro-radius (adiabatic).
- **Method & tools:** NumPy; extend `pushers.py` with a guiding-center integrator; introduce the adiabaticity parameter r_L/L.
- **You'll learn:** that the messy gyrating orbit *is* a smooth guiding-center drift plus a circle — the central simplification of magnetized-plasma theory. (This is the bridge to experiment 02.)
- **Validation:** guiding-center drift velocity matches the analytic grad-B/curvature formulas; full-orbit minus guiding-center error scales as expected with r_L/L.
- **Compute:** seconds.

### F3 — Orbits in a real tokamak field
- **Models:** particle orbits in an analytic large-aspect-ratio tokamak field (toroidal + poloidal) → trapped "banana" orbits and passing orbits.
- **Assumes:** analytic concentric-circle flux surfaces (no real equilibrium yet — that's experiment 04).
- **Method & tools:** NumPy; a toroidal-field helper added to `plasmaplay.fields`; Poincaré section at a fixed toroidal angle.
- **You'll learn:** why a purely toroidal field fails (grad-B drift separates charges → vertical E-field → everything drifts out), and why you need the poloidal twist (rotational transform) to cancel it. Banana orbits explain neoclassical transport.
- **Validation:** banana width ≈ q·r_L/√ε from theory; trapped/passing boundary at the expected pitch angle.
- **Compute:** seconds to minutes.

### F4 — Particle ensembles with collisions
- **Models:** 10⁴–10⁵ particles with Monte-Carlo pitch-angle (Lorentz) collisions → a first taste of neoclassical transport and the loss-cone leakage of a mirror.
- **Assumes:** test particles (still no self-consistent fields); collisions as a stochastic operator.
- **Method & tools:** `numba` for the vectorized push; histogram diagnostics; optional comparison to PlasmaPy collision frequencies.
- **You'll learn:** how collisions fill the loss cone and set confinement time; the statistical view that leads into kinetic theory.
- **Validation:** measured scattering rate matches the Spitzer/Lorentz collision frequency; mirror confinement time scales with the loss-cone solid angle.
- **Compute:** seconds–minutes with numba; 10⁵ particles fits easily in RAM.

## Diagnostics you'll reuse
Trajectory plots, energy-conservation trace, Poincaré sections, velocity-space (pitch-angle) histograms.

## Key references
- Chen, *Introduction to Plasma Physics and Controlled Fusion* — ch. 2 (single-particle motion).
- Birdsall & Langdon, *Plasma Physics via Computer Simulation* — Boris pusher.
- PlasmaPy `plasmapy.formulary` and `plasmapy.simulation.particle_integrators`.

## Stretch goals
- Animate the mirror bounce and the banana orbit.
- Reproduce the loss cone as a function of mirror ratio.
- Swap in the F3 field from your *own* experiment-04 equilibrium once it exists.

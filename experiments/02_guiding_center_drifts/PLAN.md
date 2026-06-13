# 02 — Guiding-center drifts — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** not started. Builds directly on experiment 01.

## The question

If a magnetized particle is really just "a fast little circle riding on a slowly
drifting center," can we throw away the circle and follow only the center? What
do those drifts look like, when is the approximation valid, and what breaks it?

## Why it matters (real devices)

The guiding-center picture is how plasma physicists actually *think*. Confinement,
neoclassical transport, and the need for rotational transform are all
guiding-center arguments. Drifts also cause real problems: the grad-B/curvature
drift in a tokamak separates charge and would destroy confinement without the
poloidal field. This experiment makes "drift" a thing you can compute, not just a
word.

## Prerequisites

Experiment 01 (at least F1; ideally F2). Concept of an adiabatic invariant helps.

## Fidelity ladder

### F0 — Analytic / sanity
- **Models:** the drift catalogue — E×B = (E×B)/B², grad-B = (mv⊥²/2qB³)(B×∇B), curvature = (mv∥²/qB²)... , polarization drift, and the magnetic moment μ as an adiabatic invariant.
- **Assumes:** uniform or slowly-varying fields; r_L ≪ L.
- **Method & tools:** arithmetic; PlasmaPy for cross-checks.
- **You'll learn:** which drifts depend on charge sign (→ currents) and which don't (→ bulk motion).
- **Validation:** reference values for a chosen field configuration.
- **Compute:** instant.

### F1 — Integrate the guiding-center equations
- **Models:** the guiding-center equations of motion (drift velocity + parallel acceleration along B) in prescribed analytic fields.
- **Assumes:** adiabatic (μ conserved); prescribed fields.
- **Method & tools:** NumPy ODE integration (RK4 is fine here — no fast gyration to conserve); a new `guiding_center.py` integrator.
- **You'll learn:** drift trajectories directly, with no gyro-circle to average over.
- **Validation:** computed drift velocity matches the F0 formulas in a controlled field gradient.
- **Compute:** seconds.

### F2 — Guiding-center vs. full orbit (the payoff)
- **Models:** run the full Boris orbit from experiment 01 and the guiding-center orbit in the *same* non-uniform field; overlay them.
- **Assumes:** adiabatic regime, but now you measure the error as you leave it.
- **Method & tools:** reuse `pushers.boris_push` + the F1 integrator; sweep the adiabaticity parameter r_L/L.
- **You'll learn:** the guiding center *is* the gyro-averaged orbit; the approximation degrades smoothly as r_L/L → 1. This is the single most important idea in the experiment.
- **Validation:** position error between the two scales like (r_L/L)ⁿ as predicted; μ conservation breaks where theory says it should.
- **Compute:** seconds.

### F3 — Drifts in tokamak geometry
- **Models:** grad-B + curvature drift in a toroidal field → the vertical drift, charge separation, and how the poloidal field's rotational transform short-circuits it. Trapped-particle banana orbits as a guiding-center phenomenon.
- **Assumes:** analytic tokamak field (concentric flux surfaces).
- **Method & tools:** NumPy; shares the toroidal-field helper with experiment 01 F3.
- **You'll learn:** *why tokamaks are shaped the way they are* — the single clearest payoff of guiding-center theory.
- **Validation:** banana width and bounce frequency match neoclassical formulas.
- **Compute:** seconds–minutes.

### F4 — Neoclassical transport flavor
- **Models:** ensemble of guiding centers + pitch-angle collisions → diffusion across flux surfaces; recover the banana / plateau / Pfirsch–Schlüter transport regimes vs. collisionality.
- **Assumes:** test particles; simplified collision operator.
- **Method & tools:** `numba`; many guiding centers; measure radial diffusivity vs. collision frequency.
- **You'll learn:** why toroidal geometry makes transport much larger than naive (classical) collisional diffusion — the neoclassical result.
- **Validation:** the D-vs-collisionality curve shows the three classic regimes with the right scalings.
- **Compute:** minutes with numba.

## Diagnostics you'll reuse
Drift-velocity vector fields, GC-vs-orbit overlays, μ-conservation traces, radial-diffusion (mean-square displacement) plots.

## Key references
- Chen, ch. 2–3. Hazeltine & Meiss, *Plasma Confinement* (guiding-center & neoclassical).
- Helander & Sigmar, *Collisional Transport in Magnetized Plasmas* (for F4).

## Stretch goals
- Animate a banana orbit projected onto a poloidal cross-section.
- Map the trapped-particle fraction vs. inverse aspect ratio.

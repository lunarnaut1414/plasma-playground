# 06 — Ideal MHD basics — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** not started. The fluid view of plasma — prerequisite for the space drive (07).

## The question

Forget individual particles: treat the plasma as a single electrically
conducting fluid threaded by magnetic field. What waves does it support, and what
happens when those flows go nonlinear and form shocks? This is
magnetohydrodynamics — the model behind most large-scale plasma behavior.

## Why it matters (real devices)

MHD governs the macroscopic stability of *every* confinement device (kink,
ballooning, tearing modes are all MHD), the dynamics of the solar wind and
magnetospheres, and the plasma acceleration in MHD thrusters (experiment 07).
The J×B force you'll meet here is literally the thrust mechanism of a magneto-
plasma drive. It's also a gateway into computational fluid dynamics done right.

## Prerequisites

Conservation laws (mass/momentum/energy), basic CFD ideas (finite volume, fluxes,
CFL condition). Vector calculus. No particle background needed.

## Fidelity ladder

### F0 — Analytic / sanity
- **Models:** the MHD wave family — Alfvén speed v_A = B/√(μ₀ρ), sound speed, and the fast/slow magnetosonic speeds; the J×B force decomposed into magnetic pressure + tension.
- **Assumes:** linear perturbations of a uniform state.
- **Method & tools:** arithmetic; Friedrichs diagram sketch.
- **You'll learn:** that magnetized fluid has *three* wave speeds, not one; magnetic tension as a restoring force.
- **Validation:** reference speeds for the F1/F2 solvers to reproduce.
- **Compute:** instant.

### F1 — Linear waves in 1-D
- **Models:** linearize the ideal MHD equations about a uniform background; launch a small perturbation and watch an Alfvén / magnetosonic wave propagate.
- **Assumes:** small-amplitude (linear), 1-D.
- **Method & tools:** NumPy; a simple linear-advection / wave-equation integrator.
- **You'll learn:** that the F0 speeds are real and measurable; how perturbations propagate along vs. across B.
- **Validation:** measured propagation speed = F0 Alfvén/magnetosonic speed.
- **Compute:** seconds.

### F2 — Nonlinear 1-D MHD: the Brio–Wu shock tube
- **Models:** full nonlinear ideal MHD conservation laws in 1-D; the Brio–Wu shock tube (the MHD analogue of Sod's tube) producing shocks, rarefactions, and a compound wave.
- **Assumes:** 1-D (so ∇·B = 0 is automatic); ideal (no resistivity); finite-volume discretization.
- **Method & tools:** NumPy finite-volume solver with an approximate Riemann solver (HLL/HLLD) or Lax–Friedrichs; `numba` for speed.
- **You'll learn:** how to write a conservative hyperbolic solver; what MHD shocks look like; the famous Brio–Wu structure.
- **Validation:** your solution matches the published Brio–Wu reference profiles at t = 0.1.
- **Compute:** seconds–minutes.

### F3 — 2-D MHD: the Orszag–Tang vortex
- **Models:** 2-D ideal MHD; the Orszag–Tang vortex — the standard test that develops MHD turbulence and current sheets.
- **Assumes:** 2-D ideal MHD; now ∇·B = 0 must be *enforced* numerically.
- **Method & tools:** extend the F2 solver to 2-D with constrained transport or divergence cleaning; `numba`; the 96 GB RAM helps for fine grids.
- **You'll learn:** the ∇·B problem (the one thing that makes MHD codes special) and how constrained transport solves it; the onset of MHD turbulence.
- **Validation:** reproduce the canonical Orszag–Tang density/pressure structure at the standard output time; ∇·B stays at round-off.
- **Compute:** minutes; numba/parallel valuable here.

### F4 — Resistive MHD / a real framework
- **Models:** add resistivity → magnetic reconnection (tearing modes, current-sheet dynamics), **or** reproduce these tests in an established framework (Dedalus, Athena++, or PLUTO) and compare.
- **Assumes:** resistive MHD or a chosen production code.
- **Method & tools:** your F3 solver + a resistive term, or an external code with Python scripting (Dedalus is pure-Python spectral and pleasant on a Mac).
- **You'll learn:** reconnection — how field topology changes despite "frozen-in" ideal MHD; what production MHD frameworks provide.
- **Validation:** tearing-mode growth rate matches the resistive-MHD scaling (∝ S^-3/5); or your result matches the reference code.
- **Compute:** minutes–hours.

## Diagnostics you'll reuse
1-D profile plots vs. reference, 2-D pseudocolor of density/pressure/|B|, ∇·B monitor, energy partition (kinetic/magnetic/thermal), wave-speed measurements.

## Key references
- Freidberg, *Ideal MHD*. Goedbloed & Poedts, *Principles of Magnetohydrodynamics*.
- Brio & Wu (1988) and Orszag & Tang (1979) — the canonical test problems.
- Dedalus (dedalus-project.org), Athena++, PLUTO docs.

## Stretch goals
- Animate the Orszag–Tang vortex forming current sheets.
- Measure all three wave speeds from a single 2-D simulation via the ω–k spectrum (links to experiment 08).

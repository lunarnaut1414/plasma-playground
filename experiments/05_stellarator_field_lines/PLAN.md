# 05 — Stellarator field lines — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** not started. The 3-D cousin of experiment 04.

## The question

A tokamak gets its twist (rotational transform) by driving a current in the
plasma. A stellarator gets it purely from the *shape of external coils* — no
plasma current needed. Can we trace magnetic field lines through such a 3-D field
and see the nested flux surfaces (and the islands) appear?

## Why it matters (real devices)

Stellarators (W7-X, NCSX, HSX) are the leading steady-state fusion concept
precisely because they need no driven current — but the price is fully 3-D
geometry and fiendish coil design. Field-line tracing + Poincaré plots are the
basic diagnostic of whether a stellarator has good flux surfaces at all. This
experiment is your entry into 3-D magnetic confinement.

## Prerequisites

Experiment 01 (field-line/orbit integration). Biot–Savart law. The idea of a
Poincaré section. Helpful: experiment 04 for the flux-surface concept.

## Fidelity ladder

### F0 — Analytic / sanity
- **Models:** rotational transform ι (iota), the field-line winding number; a toy field = uniform toroidal field + a helical perturbation.
- **Assumes:** simple analytic model field; large aspect ratio.
- **Method & tools:** arithmetic + a sketch; concept of rational vs. irrational surfaces.
- **You'll learn:** what ι means, why rational surfaces (ι = m/n) are where islands form.
- **Validation:** reference for the traced ι in F1.
- **Compute:** instant.

### F1 — Trace field lines in a model field
- **Models:** integrate the field-line ODE dx/ds = B/|B| through an analytic stellarator-like field; make a Poincaré plot at one toroidal angle.
- **Assumes:** prescribed analytic field (no real coils yet).
- **Method & tools:** NumPy/SciPy ODE integrator; a 3-D field helper in `plasmaplay.fields`.
- **You'll learn:** that field lines trace out nested surfaces; the Poincaré "puncture plot" as the fundamental 3-D diagnostic.
- **Validation:** measured ι (twist per toroidal turn) matches the F0 model value; surfaces close up on themselves.
- **Compute:** seconds.

### F2 — Coil filaments (Biot–Savart) + islands
- **Models:** build the field from discrete current-carrying coil filaments via Biot–Savart; trace lines; identify magnetic islands at rational surfaces and the onset of stochasticity.
- **Assumes:** infinitely thin coils; vacuum field (no plasma response).
- **Method & tools:** NumPy Biot–Savart (vectorized; numba if slow); Poincaré analysis to find island chains.
- **You'll learn:** how coil geometry *creates* the rotational transform; that imperfect coils break good surfaces into islands and chaos — the central stellarator design problem.
- **Validation:** island width and location match the resonant-perturbation estimate; ι profile is smooth where surfaces are intact.
- **Compute:** seconds–minutes (Biot–Savart over many points benefits from numba).

### F3 — A real stellarator (simsopt)
- **Models:** load a real coil set — W7-X or the NCSX/`simsopt` example coils — and trace/Poincaré the actual vacuum field; compute a quasi-symmetry metric.
- **Assumes:** vacuum or a provided equilibrium field; published/example coil data.
- **Method & tools:** **simsopt** (`pip install simsopt`) — field-line tracing, Poincaré, and quasi-symmetry tools are built in.
- **You'll learn:** what real stellarator flux surfaces look like (W7-X's bean-shaped cross-sections); what "quasi-symmetry" means and why it matters for confinement.
- **Validation:** reproduce a simsopt example's Poincaré plot; quasi-symmetry residual in the expected range.
- **Compute:** minutes.

### F4 — 3-D equilibrium or coil optimization
- **Models:** either (a) compute a 3-D MHD equilibrium with **DESC** (JAX-based) or VMEC, or (b) run a small **simsopt** optimization that reshapes coils to improve quasi-symmetry.
- **Assumes:** depends on path; this is the actual research frontier of stellarator design.
- **Method & tools:** DESC (JAX — CPU fine on M2, GPU optional) or simsopt's optimization drivers.
- **You'll learn:** that stellarator design *is* an optimization problem over coil/boundary shapes; how equilibria and coils are co-designed.
- **Validation:** optimized configuration shows measurably lower quasi-symmetry residual than the start; DESC equilibrium matches a reference case.
- **Compute:** minutes–hours; DESC on JAX-CPU is reasonable, JAX-Metal experimental.

## Diagnostics you'll reuse
Poincaré puncture plots, ι(radius) profiles, island-width measurements, quasi-symmetry residual maps, 3-D coil/field-line renderings.

## Key references
- Boozer, "Physics of magnetically confined plasmas" (Rev. Mod. Phys.) — stellarator theory.
- Imbert-Gérard, Paul & Wright, *An Introduction to Stellarators* (arXiv:1908.05360).
- simsopt docs (simsopt.readthedocs.io) and DESC docs (desc-docs.readthedocs.io).

## Stretch goals
- Side-by-side Poincaré of a tokamak (experiment 04) vs. a stellarator — same diagnostic, very different geometry.
- Deliberately misalign a coil and watch good surfaces break into islands.

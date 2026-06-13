# Fundamentals & common functions — prove these out first

Reading across all eight experiment plans, the same handful of **numerical
kernels** show up again and again. Building and *validating* those shared pieces
first — bottom-up, each with a falsifiable test — means every later experiment
stands on already-trusted code instead of re-deriving (and re-debugging) the same
machinery.

This doc is the build order. It lists (1) the common functions that belong in the
`plasmaplay` package, and (2) a prioritized validation suite — **1-D tests first,
then 2-D** — where each test pins one kernel against an analytic result or a
published reference.

> Rule: a kernel is not "done" until its validation test passes. See the
> validation philosophy in [`FIDELITY.md`](FIDELITY.md).

---

## 1. The shared kernels (what goes in `plasmaplay/`)

Ordered by how many experiments depend on them — build the high-leverage ones first.

| Kernel | Lives in | Used by | Status |
|--------|----------|---------|--------|
| **Constants & formulary** (ω_c, r_L, ω_pe, λ_D, v_A) | `constants.py` | all | ✅ done |
| **Analytic field models** (uniform, mirror, toroidal, coil) | `fields.py` | 01, 02, 05, 07 | ✅ partial |
| **Boris pusher** (particle orbit) | `pushers.py` | 01, 02 | ✅ done |
| **RK4 / ODE integrator** (field lines, guiding center) | `integrators.py` | 02, 05 | ✅ done |
| **Grid ↔ particle weighting** (CIC deposit/interpolate) | `pic.py` | 03 | ✅ done |
| **FFT Poisson solver** (1-D & 2-D periodic) | `solvers.py` | 03 | ✅ done |
| **Finite-difference elliptic solver** (Laplacian / Grad–Shafranov Δ*) | `solvers.py` | 04 | ✅ done |
| **Finite-volume hyperbolic solver** (1-D → 2-D, Riemann) | `fvm.py` *(new)* | 06, 07 | ☐ |
| **∇·B control** (constrained transport / cleaning) | `fvm.py` *(new)* | 06, 07 | ☐ |
| **Biot–Savart** (field from coil filaments) | `fields.py` | 05 | ✅ done |
| **Distribution loaders** (Maxwellian, two-beam) | `pic.py` | 03 | ✅ done |
| **Collision operator** (Monte-Carlo pitch-angle) | `collisions.py` *(new)* | 01, 02 | ☐ |
| **Plasma dispersion function** Z(ζ) | `dispersion.py` *(new)* | 08 | ☐ |
| **Spectral diagnostics** (ω–k FFT, Poincaré section) | `diagnostics.py` | 03, 05, 06, 08 | ◐ Poincaré, ι, dominant-frequency done; ω–k pending |
| **Conservation monitors** (energy, μ, ∇·B) | `diagnostics.py` *(new)* | all | ☐ |

Design conventions these must follow (so they compose):
- **SI units everywhere.** Use `astropy.units`/`plasmapy` only for cross-checks.
- **Fields are callables** `f(position) -> (3,) ndarray`, so any field drops into any
  pusher/tracer (already the pattern in `fields.py`).
- **Solvers take arrays in, arrays out** — no global state, easy to test.
- **Every kernel ships with its validation test** (below) before any experiment uses it.

---

## 2. Validation suite — 1-D first, then 2-D

Each test has: the kernel it proves, the setup, and the **pass criterion** (a number,
not a vibe). Suggested home: a top-level `tests/` dir runnable with `pytest`, plus a
plot-producing version in the relevant experiment.

### Tier 0 — analytic anchors (no simulation)

- **V0 · Formulary** *(✅ in `tests/test_constants.py`)* — `constants.py` outputs match a
  first-principles SI computation, canonical reference numbers, and (optionally)
  `plasmapy.formulary` for ω_c, r_L, ω_pe, λ_D, v_A. *Pass: < 1% vs reference.*

### Tier 1 — 1-D / single-particle kernels

- **V1 · Boris energy conservation** *(✅ in exp 01)* — proton in uniform B over many
  orbits. *Pass: relative energy drift < 1e-10; radius = r_L.*
- **V2 · E×B drift** *(✅ in exp 01)* — crossed E, B. *Pass: measured guiding-center
  drift = E/B within 1%.*
- **V3 · ODE integrator accuracy** *(✅ in `tests/test_integrators.py`)* — RK4 on the
  harmonic oscillator. *Pass: log-error vs log-dt slope in (3.8, 4.2) — 4th order.*
- **V4 · FFT Poisson 1-D** *(✅ in `tests/test_solvers.py`)* — periodic source vs analytic
  φ. *Pass: round-off for any grid-resolved source (spectral solver is exact per mode,
  stronger than the 2nd-order a finite-difference solve would give).*
- **V5 · Cold-plasma oscillation** *(✅ in `tests/test_pic_physics.py`)* — perturbed cold
  plasma in 1-D PIC. *Pass: measured frequency = ω_pe within 5% (measured 1.005),
  independent of mode.*
- **V6 · Landau damping** *(✅ in `tests/test_pic_landau.py`)* — Maxwellian Langmuir wave,
  k λ_D = 0.5. *Pass: measured γ matches the Landau rate 0.1533 ω_pe within ~35%
  (PIC noise); measured ≈ 0.17.*
- **V7 · Two-stream growth rate** *(✅ in `tests/test_pic_landau.py`)* — counter-streaming
  beams. *Pass: linear growth rate matches the cold-dispersion root within ~30%
  (theory 0.353, measured ≈ 0.32); plus marginal-stability checks.*
- **V8 · 1-D shock tubes** — finite-volume solver on Sod (hydro sanity) then **Brio–Wu**
  (MHD). *Pass: profiles overlay the published reference at the standard output time.*
- **V9 · Kinetic dispersion** — Z(ζ) root-find gives Bohm–Gross frequency **and**
  Landau damping. *Pass: matches `plasmapy.dispersion` and V6.*

### Tier 2 — 2-D / field kernels

- **V10 · FFT Poisson 2-D** *(✅ in `tests/test_solvers.py`)* — ρ = sin(k_x x)·sin(k_y y)
  vs analytic φ. *Pass: round-off (spectral).*
- **V11 · Biot–Savart loop** *(✅ in `tests/test_fields.py`)* — field of a single circular
  current loop vs the on-axis analytic formula B_z(0,0,z). *Pass: < 0.1% on axis.*
- **V12 · Grad–Shafranov** *(✅ in `tests/test_solvers.py`)* — fixed-boundary Δ* solver
  vs a manufactured exact solution (the Solov'ev idea). *Pass: 2nd-order convergence
  (slope in (1.8, 2.2)); < 0.1% on a fine grid.*
- **V13 · Field-line ι** *(✅ in `tests/test_diagnostics.py`)* — trace the screw-pinch
  field; Poincaré → rotational transform. *Pass: measured ι = analytic within 1%;
  surfaces close (punctures at constant radius).*
- **V14 · Orszag–Tang** — 2-D MHD vortex. *Pass: density/pressure structure matches the
  canonical reference; **∇·B stays at round-off** (the real test of the CT scheme).*
- **V15 · ω–k recovery** — FFT field data from V5/V6/V8 runs into the ω–k plane.
  *Pass: power ridge falls on the analytic dispersion curve (closes the loop to exp 08).*

---

## 3. Recommended build order

A path that makes each step usable immediately and keeps you on already-tested ground:

1. ✅ **Finish Tier 0 + V1–V2** — formulary cross-check; Boris pusher locked in. *Foundation for 01/02.*
2. ✅ **`integrators.py` + V3**, then **Biot–Savart + V11** — unlocked experiment 05 field-line tracing (V13). (Guiding-center / experiment 02 still to do.)
3. ✅ **`solvers.py` Poisson (V4, V10)** + **`pic.py` weighting/loaders** — unlocked the whole PIC experiment 03 (V5 → V6 → V7).
4. ✅ **`solvers.py` elliptic + V12** — unlocked tokamak equilibrium (experiment 04).
5. **`fvm.py` 1-D (V8)** → **2-D + ∇·B (V14)** — unlocks MHD (06) and the space drive (07). ← next
6. **`dispersion.py` + Z (V9)** and **ω–k diagnostic (V15)** — unlocks experiment 08 and ties the kinetic/fluid pictures together.

Kernels 1–3 are pure NumPy/SciPy, **seconds on one core**. Steps 5–6 are where
`numba` and the 96 GB RAM start to earn their keep (see per-experiment compute notes).

---

## 4. What "common functions" buys you

- **One bug fix, everywhere.** A validated Poisson solver serves PIC *and* (in elliptic
  form) the tokamak equilibrium. The ω–k FFT serves PIC, MHD, and dispersion.
- **Cross-experiment reuse is already designed in.** Exp 01's banana orbits want exp
  04's equilibrium field; exp 08 measures dispersion from exp 03's PIC output; exp 07
  reuses exp 06's MHD solver. Shared kernels make those couplings trivial instead of
  copy-paste.
- **Trust.** Because each kernel has a number it must hit, when an experiment shows
  something surprising you know it's physics, not a bug in machinery you never tested.

# 10 — Tokamak MHD stability — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** B1 implemented (cylinder linear stability), kernel
> `plasmaplay/cylinder_mhd.py`, tests `tests/test_cylinder_mhd.py` (14 passing).
> B2 (nonlinear island saturation) and B3 (sawtooth cycle) next.

## The question

A tokamak does not just sit there and burn — its current-carrying plasma is prone to
**MHD instabilities** (kinks, tearing modes, sawteeth) that move the fluid on the
microsecond Alfvén timescale. Which ones go unstable, and what sets the threshold?

## Why a cylinder first

The fully toroidal, shaped, nonlinear problem is JOREK/NIMROD territory. The
**periodic cylinder** (straight tokamak) keeps the essential ingredients — a safety
factor `q(r)`, rational surfaces, the kink and tearing branches — while staying
analytically anchored, so every number can be checked. (Honest boundary: no
toroidal coupling of poloidal harmonics, no shaping, no ballooning — those are the
later, toroidal rungs.)

## Fidelity ladder

### B1 — Cylindrical linear stability  ✅ implemented (`run.py`)
- **Models:** `plasmaplay/cylinder_mhd.py` — the screw-pinch q-profile
  `screw_pinch_q`, the rational surface `rational_surface`, the outer Newcomb
  equation and tearing index `delta_prime_cylinder`, the m=1 internal-kink criterion
  `internal_kink_unstable` + displacement `internal_kink_xi`, and the FKR growth rate
  `fkr_growth_rate` (reusing the slab-layer law from `tearing.py` / T4).
- **Validation:** q(0)=q0, q(a)=(ν+1)q0; the m=1 kink is unstable **iff q(0)<1**; the
  **sign of Δ′** predicts tearing stability and Δ′ falls with m; `γ ∝ S^(−3/5)`.
  *(14 tests)* Honest caveat: the *absolute* Δ′ is resolution-dependent near the
  singular layer (a known feature); the **sign and ordering** are the robust outputs.
- **Deliverable:** `outputs/kink_eigenmode.gif` — the m=1 internal kink: the core
  shifting into the characteristic crescent; `kink_eigenmode.png` still.

### B2 — Nonlinear reduced-MHD island saturation  ◻ not yet
- **Models:** evolve ψ and vorticity in (r,θ) for one n with resistivity; watch a
  tearing island grow and **saturate** (Rutherford `dW/dt ∝ Δ′(W)`).
- **Validation:** linear phase matches B1 growth; island width saturates at Δ′(W)=0.
- **Deliverable:** `tearing_island_saturation.gif`.

### B3 — The sawtooth cycle (Kadomtsev)  ◻ not yet
- **Models:** when q(0)<1, an m=1 reconnection flattens the core (helical-flux
  conserving); q(0) relaxes >1, resistive diffusion re-peaks it, repeat.
- **Validation:** helical-flux conservation; period scales with the resistive time.
- **Deliverable:** `sawtooth_cycle.gif`.

## Toward the coupling (Track C)

B1–B3 are the MHD-event library. Track C of the overnight charter couples them into
the experiment-09 transport burn: monitor q(0) and Δ′ from the evolving profiles and
fire a sawtooth crash / tearing island as a profile-redistribution *event* during the
discharge — the staged two-timescale "flight simulator".

## Key references
- Newcomb, *Ann. Phys.* 10, 232 (1960) — the cylindrical stability equation.
- Furth, Killeen & Rosenbluth, *Phys. Fluids* 6, 459 (1963) — tearing & Δ′.
- Wesson, *Tokamaks* — kink/tearing/sawtooth; the q(0)<1 internal-kink trigger.

# 10 — Tokamak MHD stability — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** B1 (cylinder linear stability, `cylinder_mhd.py`, 14 tests), B2
> (nonlinear reduced-MHD tearing island + Rutherford saturation, `reduced_mhd.py`,
> 4 tests), and B3a (Kadomtsev reconnection, `sawtooth.py`, 4 tests) implemented. The
> tuned periodic sawtooth cycle (B3b) and Track C (MHD↔transport coupling) are next.

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

### B2 — Nonlinear reduced-MHD island + Rutherford saturation  ✅ implemented (`run.py --island`)
- **Models:** `plasmaplay/reduced_mhd.py` — the Strauss reduced-MHD equations (ψ and
  vorticity U=∇²φ) on a 2-D slab (x finite-difference, y spectral, vectorized
  FFT+tridiagonal elliptic solve for φ, SSP-RK2), on the Harris sheet of T4. A tearing
  mode grows exponentially, reconnects an island, and **saturates** (Rutherford).
- **Validation:** the elliptic inversion is exact; a seeded mode grows for ka<1 and
  decays for ka>1; the **linear growth obeys γ ∝ S^(−3/5)** (measured exponent −0.58,
  FKR); the **island saturates** — dW/dt rises, **peaks, then declines** (W → W_sat ~ 2
  sheet widths) rather than growing exponentially. *(4 tests)*
- **Deliverable:** `outputs/tearing_island_saturation.gif` — W(t) bending over beside
  the flux contours of the reconnecting/saturating island; `..._saturation.png` still.
- **Scope (honest):** the absolute linear growth rate is ~0.54× the T4 eigenvalue (an
  O(1) convention difference — the scaling/threshold are asserted, not the value); the
  asymptotic plateau is approached on the slow resistive timescale and the exact W_sat
  is influenced by the wall. The Δ′(W)→0 Rutherford form is shown qualitatively via the
  dW/dt turnover, not fit to the analytic Rutherford coefficient.

### B3 — The Kadomtsev sawtooth  ◧ B3a reconnection done; periodic cycle = B3b
- **Models:** `plasmaplay/sawtooth.py` — the helical flux ψ*(r), the Kadomtsev mixing
  radius, the energy-conserving `kadomtsev_flatten`, and a `SawtoothCycle` (resistive
  induction re-peaking + the crash). When q(0)<1 the m=1 reconnection flattens the core
  (q→1) and the temperature, conserving thermal energy.
- **Validation (B3a):** ψ* peaks at the q=1 surface; r_mix sits outside it; the flatten
  conserves ∫T·r dr exactly; a single crash flattens T (energy conserved to 2e-16),
  reconnects the helical flux core to ~0, and resets q(0)→1. *(4 tests)*
- **Deliverable (B3a):** `outputs/sawtooth_crash.png` — before/after a single crash
  (`run.py --sawtooth`).
- **B3b (next):** the tuned **periodic cycle** — the crashes are presently
  over-frequent and the period scales only weakly with τ_R (~τ_R^0.6); a clean
  period ∝ τ_R sawtooth and `sawtooth_cycle.gif` need a re-peaking timescale set by
  the global resistive diffusion rather than the fast near-axis dynamics.

## Toward the coupling (Track C)

B1–B3 are the MHD-event library. Track C of the overnight charter couples them into
the experiment-09 transport burn: monitor q(0) and Δ′ from the evolving profiles and
fire a sawtooth crash / tearing island as a profile-redistribution *event* during the
discharge — the staged two-timescale "flight simulator".

## Key references
- Newcomb, *Ann. Phys.* 10, 232 (1960) — the cylindrical stability equation.
- Furth, Killeen & Rosenbluth, *Phys. Fluids* 6, 459 (1963) — tearing & Δ′.
- Wesson, *Tokamaks* — kink/tearing/sawtooth; the q(0)<1 internal-kink trigger.

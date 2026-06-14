# B3 — The Kadomtsev sawtooth crash (reconnection of the core)

The third rung of the MHD-instability track (NIGHT.md Track B). B1 found *when* the
m=1/n=1 internal kink goes unstable (q(0) < 1); B3 is the **reconnection** that fires
once it does — the Kadomtsev crash that flattens the core and resets q(0) toward 1,
the elementary event of the sawtooth oscillation.

**Status: B3a — the reconnection model is built and validated; the tuned periodic
cycle (and the period ∝ τ_R scaling + `sawtooth_cycle.gif`) is the follow-on B3b.**

## What was built (`plasmaplay/sawtooth.py`)

- `helical_flux(r, q)` — the (1,1) helical flux `ψ*(r) = ∫₀ʳ B_θ(1−q) dr'` (B_z=R=1):
  rises from 0 on axis, **peaks exactly at the q=1 surface** (the integrand changes
  sign there), then falls. The non-monotonic profile that reconnects.
- `mixing_radius(r, q)` — the Kadomtsev mixing radius r_mix where ψ* returns to its
  axis value (0); reconnection mixes the core out to here. `None` if q(0) ≥ 1.
- `kadomtsev_flatten(r, field, r_mix)` — flattens a field inside r_mix to its
  area-weighted mean, **conserving ∫field·r dr exactly** (so flattening a temperature
  conserves its thermal energy).
- `SawtoothCycle` — a 1-D model: the poloidal field evolves by the **resistive
  induction equation** with a core-peaked conductivity (current peaks on axis, driving
  q(0) below 1), and a `_crash()` applies the Kadomtsev flatten when q(0) crosses the
  trigger, resetting q→~1 and flattening T. Reuses the B1 screw-pinch q-profile.

## Validation (falsifiable — `tests/test_sawtooth.py`, 4 tests)

- **ψ* peaks at the q=1 surface** (to grid resolution).
- **r_mix lies outside the q=1 surface** and inside the wall; it is absent (no crash)
  when q(0) ≥ 1.
- **The flatten conserves energy exactly** (∫field·r dr to 1e-12) and leaves a flat core.
- **A single crash reconnects the core:** q(0) resets to ~1, the core temperature
  flattens (std → 0), the **helical flux in the core goes to ~0** (ψ*_max 0.0198 → 0),
  and **thermal energy is conserved to 2×10⁻¹⁶**.

## Deliverable

`outputs/sawtooth_crash.png` — before/after a single Kadomtsev crash: T flattening, q
reconnecting to 1, and the helical flux ψ* collapsing in the core
(`run.py --sawtooth`).

## Scope boundary (stated honestly)

- **B3a vs B3b.** The reconnection event is validated. The full **periodic cycle**
  runs (q(0) and the core T do oscillate), but the crashes are over-frequent and the
  **period scales only weakly with τ_R** (~τ_R^0.6, not linearly) in the present
  driver — the central current re-peaks faster than the global resistive time. A
  clean, well-separated sawtooth with period ∝ τ_R, and the `sawtooth_cycle.gif`, are
  the follow-on rung **B3b** (needs a better-separated re-peaking timescale, e.g. a
  current-diffusion model where q(0) relaxes on τ_R rather than locally).
- **Reduced model.** 1-D cylinder, full Kadomtsev "complete reconnection" idealisation
  (q→1 flat core with the residual current sheet at r_mix not resolved); the exactly
  conserved, tested invariant is **thermal energy**, and the helical flux is shown to
  reconnect (it is the resistive process, not a conserved quantity).

## References
- Kadomtsev, *Sov. J. Plasma Phys.* 1, 389 (1975) — the sawtooth reconnection model.
- von Goeler, Stodiek & Sauthoff, *Phys. Rev. Lett.* 33, 1201 (1974) — sawtooth discovery.
- Wesson, *Tokamaks* — the sawtooth cycle and the q(0) < 1 trigger.

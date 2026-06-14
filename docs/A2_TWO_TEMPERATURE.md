# A2 / F2.5 — Two temperatures (Te, Ti) + Spitzer equipartition

The burning-plasma transport rung that stops pretending the electrons and ions
share a temperature (NIGHT.md Track A, rung A2). Single-temperature transport (F2)
is a useful first cut, but it hides three facts that matter once you ask *which*
species is hot:

1. **Fusion is an ion reaction.** The Bosch–Hale `<σv>` and therefore the alpha
   self-heating follow the **ion** temperature `T_i`, not a blended `T`.
2. **Bremsstrahlung is an electron loss.** It follows `T_e`.
3. **The two species exchange energy only collisionally**, on the Spitzer
   equipartition time — which at fusion conditions is *hundreds of milliseconds*,
   comparable to `τ_E`. So if the heating is lopsided, the two temperatures stay
   apart long enough to matter.

Heating *is* lopsided in real machines: neutral beams and ion-cyclotron RF deposit
on the ions, while fusion alphas (3.5 MeV) slow down predominantly on the
electrons (`f_αe ≈ 0.85`). The result is the regime beam-heated tokamaks actually
live in — **`T_i > T_e`** (TFTR/JET "supershots" reached `T_i/T_e ~ 2–3`).

## What was built (`plasmaplay/transport.py`)

- `coulomb_logarithm(n_e, T_e)` — NRL Plasma Formulary `lnΛ = 24 − ln(√n_e[cm⁻³]/T_e[eV])`.
- `collision_frequency_ei(n_e, T_e, z_eff)` — Spitzer electron collision rate,
  `ν_e = 2.91e-6·Z·n_e[cm⁻³]·lnΛ·T_e[eV]^(−3/2)` (NRL). The clock for energy exchange.
- `equipartition_power(n_e, T_e, T_i)` — Braginskii electron→ion exchange power
  density `Q_Δ = 3(m_e/m_i)n_e ν_ei k_B(T_e − T_i)` [W/m³]; positive when electrons
  are hotter, zero at `T_e = T_i`. `m_i = μ_i·m_p`, `μ_i ≈ 2.5` amu for 50:50 D-T.
- `equipartition_time(n_e, T_e)` — `τ_eq = 1/(4(m_e/m_i)ν_ei)`, the e-folding time
  of the temperature *difference* at `n_e = n_i`. Scales as `T_e^{3/2}/n_e`.
- `two_temperature_relax_0d(...)` — a 0-D relaxation (only `Q_Δ`, no heating/transport);
  the clean validation anchor.
- `TwoTempTransport1D(Transport1D)` — evolves `T_e(ρ,t)`, `T_i(ρ,t)`, `n(ρ,t)` with
  separate `χ_e`, `χ_i`, equipartition coupling, and a heating mix
  (`p_aux_e`, `p_aux_i`, `frac_alpha_e`). Reuses the parent's implicit (Thomas)
  diffusion per channel and its `_vol_avg` diagnostics.

  ```
  (3/2)∂(n T_e)/∂t = ∇·(χ_e ∇T_e) + p_aux_e + f_αe p_α − p_brem − Q_Δ
  (3/2)∂(n T_i)/∂t = ∇·(χ_i ∇T_i) + p_aux_i + (1−f_αe) p_α        + Q_Δ
  ```

## Validation (falsifiable — `tests/test_transport.py`, 7 tests)

- **Equipartition-time magnitude:** `τ_eq(n=1e20, T_e=10 keV, μ=2.5) = 231 ms`,
  matching a hand computation from the NRL `ν_e` formula (~230 ms).
- **Scaling:** `τ_eq ∝ T_e^{3/2}/n_e` — doubling `n_e` ≈ halves it (the residual ~2%
  is `lnΛ`'s own `√n` dependence, kept honest in the test), doubling `T_e` raises it ~2.7×.
- **Coupling sign:** `Q_Δ > 0` when `T_e > T_i`, `< 0` when colder, exactly 0 when equal.
- **0-D relaxation:** an unheated `(12, 4) keV` split relaxes to the energy-conserving
  mean (8 keV); the difference decays at exactly `1/τ_eq` measured from the first step.
- **1-D coupling & separation:** the solver drags a split together when unheated, and
  sustains `T_i0 = 24 keV > T_e0 = 13 keV` (`T_i/T_e = 1.89`) under NBI ion heating.

## Deliverables

- `outputs/burn_1d_two_temperature.gif` — `T_e(ρ,t)` and `T_i(ρ,t)` profiles
  separating as the beam ramps onto the ions (regen: `python gif_gallery.py
  burn_1d_two_temperature`).
- `outputs/burn_1d_two_temperature.png` still + a 0-D relaxation panel
  (`python experiments/09_burning_plasma/run.py --mode twotemp --save`).

## Scope boundary (what this is *not*)

The 1-D transport model still has **no β-limit** (that is rung A4), so with good
confinement it would run away hot; the showcase scenario is deliberately a
**beam-heated, sub-ignition** discharge (higher χ) that sits at a realistic
`T_i ~ 24 keV` instead. Transport coefficients `χ_e`, `χ_i` are still *prescribed*
inputs, not computed from a turbulence model (that is F4). The equipartition term is
the standard Braginskii/Spitzer form with a constant-`μ_i` D-T average mass; it does
not resolve separate D and T channels or fast-ion slowing-down spectra.

## References
- NRL Plasma Formulary — collision rates, `lnΛ`, equipartition.
- Braginskii, *Reviews of Plasma Physics* 1 (1965) — the transport closure.
- Wesson, *Tokamaks* — `T_i ≠ T_e` in beam-heated discharges.

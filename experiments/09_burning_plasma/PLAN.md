# 09 — Burning plasma (transport) — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** F0 + F1 + F2 + F2.5 implemented (`run.py`), kernel
> `plasmaplay/transport.py`, tests `tests/test_transport.py` (29 passing). F3
> (couple to the experiment-04 equilibrium) next.

## The question

A tokamak discharge has an *arc*: you heat the plasma, fusion alpha particles
start self-heating it, it ignites and settles into a burning steady state, and
you feed it fuel to keep it going. What does that arc look like, and what sets
whether it ignites at all?

## Why this is a transport problem, not MHD/CFD

"Ignition → steady state → fuel injection" plays out over the **energy- and
particle-confinement timescale (~seconds)**. That is the domain of *transport*
modelling — evolving the radial temperature and density **profiles** under
sources (heating, fusion α, fuelling) and sinks (transport, radiation). It is a
different simulation from MHD/CFD, which resolves the plasma *fluid motion* on the
**Alfvén timescale (microseconds–milliseconds)** — flows, waves, instabilities.
You cannot march an MHD/CFD code across a whole discharge; the right tool for the
arc is a transport code (TRANSP / ASTRA / RAPTOR). This experiment is their toy
cousin. (The MHD/CFD layer is experiment 06 and the T4 rung of the 3-D guide.)

## Prerequisites

0-D power balance (Lawson criterion). The D-T reaction and its 17.6 MeV split
(3.5 MeV α stays and heats; 14.1 MeV neutron escapes). 1-D diffusion / implicit
finite differences. A flux-surface-averaged radial coordinate ρ = r/a.

## Fidelity ladder

### F0 — 0-D burn dynamics / Lawson  ✅ implemented (`--mode zerod`)
- **Models:** two coupled ODEs for the volume-averaged plasma energy W = 3nT and
  density n: `dW/dt = P_aux + P_α − P_brem − W/τ_E`, `dn/dt = S_fuel − n/τ_p`.
- **Assumes:** zero-dimensional (one number per quantity); τ_E lumps all transport
  into one confinement time; 50:50 D-T, single temperature, quasineutral.
- **Method & tools:** RK4; Bosch-Hale `<σv>`; bremsstrahlung loss. Pure NumPy.
- **You'll learn:** ignition as a *threshold* — α self-heating overtaking losses —
  and the thermally-stable burning point it runs away to; the Lawson triple product.
- **Validation:** ignites above the Lawson triple product (n T τ_E ≳ 3×10²¹
  keV·s·m⁻³) and dies below it; at steady state P_α = P_loss + P_brem. *(tests)*
- **Compute:** instant.

### F1 — 0-D with He ash, dilution & β-limit  ✅ implemented (`--mode ash`)
- **Models:** `transport.burn_0d_ash` — three coupled ODEs for fuel-ion density
  n_DT, helium-ash density n_He, and energy W. Ash is born one-per-reaction and
  pumped on τ_He*; fuel is burned (−2·R_fus) and refuelled; quasineutrality
  n_e = n_DT + 2n_He gives fuel **dilution**; ash raises **Z_eff** (more brem); a
  **soft one-sided β-limit** degrades confinement only above β_limit (Troyon),
  pinning the operating point in the real burning band.
- **Method & tools:** RK4 on (n_DT, n_He, W); Bosch-Hale `reaction_rate_dt`;
  `beta_thermal` / `troyon_limit` helpers. Pure NumPy.
- **You'll learn:** why machines must pump ash and refuel continuously, why Z_eff
  matters, and why β-limited burns sit at ~15 keV not ~80 keV.
- **Validation:** steady ash balance n_He = τ_He*·R_fus (0.99); dilution lowers
  P_fusion; the β-limit pins β at the Troyon limit and lands T in the 10–25 keV
  band (14.3 keV vs 45 keV unlimited). *(6 tests)*
- **Deliverable:** `outputs/burn_0d_ignition.gif` — the (n,T) ignition track
  colored by ash fraction; `burn_0d_ash.png` still.
- **Compute:** instant.

### F2 — 1-D radial transport  ✅ implemented (default mode)
- **Models:** evolve T(ρ,t) and n(ρ,t) with diffusion equations,
  `(3/2)∂(nT)/∂t = ∇·(nχ∇T) + p_aux + p_α − p_brem`, `∂n/∂t = ∇·(D∇n) + S_fuel`,
  in circular (ρ-weighted) geometry. Scripted three-phase discharge:
  **ignition** (heating ramp) → **steady burn** (heating off, self-sustained) →
  **fuel injection** (a deep pellet pulse; density and fusion power respond).
- **Assumes:** transport coefficients χ, D are *prescribed* (the anomalous
  transport is an input, not computed); fixed-edge pedestal; single temperature;
  geometry is a large-aspect-ratio circular column (no real flux surfaces yet);
  **no β-limit**, so the stable burning point sits hotter (~80 keV) than a real
  β-limited machine (~15–25 keV).
- **Method & tools:** backward-Euler (implicit, tridiagonal/Thomas) diffusion so
  conduction doesn't force a tiny timestep; explicit reaction/heating sources.
- **You'll learn:** how profiles (not just averages) ignite and peak; how pellet
  (deep) vs gas-puff (edge) fuelling deposit differently; that fusion power ∝ n²
  so density control is a fusion-power knob.
- **Validation:** χ tuned to τ_E ≈ a²/(5.78χ); core heats and the profile is
  monotonically peaked; fuelling raises ⟨n⟩; reduces to the F0 power balance when
  volume-averaged. *(tests)*
- **Compute:** seconds; pure NumPy, single core.

### F2.5 — Two-temperature (Te, Ti) + heating mix  ✅ implemented (`--mode twotemp`)
- **Models:** `transport.TwoTempTransport1D` evolves T_e(ρ,t), T_i(ρ,t), n(ρ,t) with
  separate heat diffusivities χ_e, χ_i and a collisional energy-exchange term
  `equipartition_power` Q_Δ = 3(m_e/m_i)n_e ν_ei k_B(T_e−T_i) (Braginskii/Spitzer).
  The split single-T transport hides: **fusion follows T_i** (it is an ion
  reaction), **bremsstrahlung follows T_e**, fusion alphas heat mostly electrons
  (f_αe≈0.85) while **neutral-beam / RF-ion heating heats the ions** — so a
  beam-heated plasma runs T_i > T_e until equipartition closes the gap.
- **Method & tools:** reuses the parent's implicit (Thomas) diffusion per channel;
  explicit sources; the Spitzer ν_e is the NRL Plasma-Formulary collision rate.
  `equipartition_time` = 1/(4(m_e/m_i)ν_ei) and `two_temperature_relax_0d` (a 0-D
  relaxation) are the validation anchors.
- **You'll learn:** why real plasmas run T_i ≠ T_e; the equipartition time that sets
  how fast the channels couple (∝ T_e^{3/2}/n_e); why beam-heated discharges
  (TFTR/JET supershots) reach T_i/T_e ~ 2.
- **Validation:** τ_eq @ n=1e20, Te=10 keV = 231 ms (matches the Spitzer formula,
  ~230 ms); a 0-D split relaxes to the energy-conserving mean at the formula rate;
  the 1-D beam-heated discharge settles at T_i0=24 keV > T_e0=13 keV (T_i/T_e=1.89).
  *(7 tests)*
- **Deliverable:** `outputs/burn_1d_two_temperature.gif` — T_e(ρ,t) and T_i(ρ,t)
  profiles separating as the NBI heats the ions; `burn_1d_two_temperature.png` still.
- **Compute:** seconds; pure NumPy.

### F3 — On the real equilibrium  ◻ not yet (the next rung)
- **Models:** run the F2 transport on flux-surface-averaged coordinates of the
  experiment-04 Grad–Shafranov equilibrium — real ψ(R,Z), V'(ψ), ⟨|∇ρ|²⟩ metrics,
  a D-shaped cross-section. Recompute the equilibrium as the pressure evolves.
- **Reuses:** `solvers.grad_shafranov_solve` (exp 04), the `tokamak.py` field bridge.
- **Validation:** flux-surface geometry matches the equilibrium; energy confinement
  consistent with an L/H scaling for the device parameters.
- **Compute:** minutes.

### F4 — Predictive transport / real code  ◻ not yet
- **Models:** replace prescribed χ, D with a transport model (a critical-gradient /
  TGLF-like closure), or drive an external integrated-modelling code (ASTRA/RAPTOR).
- **Validation:** reproduce a published discharge's profile evolution or Q.
- **Compute:** minutes–hours.

## Toward MHD/CFD (the "both, staged" plan)

This experiment is the **transport** half. The **MHD/CFD** half — actually
watching the plasma fluid move and go unstable (sawtooth, tearing, ELM) — lives in
experiment 06 (`fvm.py`) and the T4 rung of [`docs/3D_TOKAMAK_GUIDE.md`]. The
staged goal is to let an MHD instability fire as an *event* during the burn (e.g.
a sawtooth crash that flattens the core profile), coupling the two timescales.

## Diagnostics you'll reuse
Triple-product / Lawson check, Q = P_fusion/P_aux, profile snapshots, the poloidal
cross-section render, volume-integrated power balance.

## Key references
- Lawson, *Proc. Phys. Soc. B* 70, 6 (1957) — the ignition criterion.
- Bosch & Hale, *Nucl. Fusion* 32, 611 (1992) — the D-T reactivity used here.
- Wesson, *Tokamaks* — confinement scalings, burn physics.
- ASTRA / RAPTOR / TRANSP — the real integrated-transport codes this mimics.

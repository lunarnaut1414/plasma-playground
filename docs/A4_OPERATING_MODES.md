# A4 / F3.5 — Operational limits & confinement modes

The rung that turns the single happy burn into an **operating window** (NIGHT.md
Track A, rung A4). A real tokamak does not run wherever you like: hard limits bound
the (density, power) plane and a confinement bifurcation splits it into L-mode and
H-mode. This rung adds the three pieces of physics that decide *whether and how* a
discharge runs — and fixes the "no β-limit" caveat carried since A2/A3.

## What was built

- `plasmaplay/operating_limits.py`
  - `greenwald_density(Ip, a)` → n_G = Ip/(πa²)·1e20 — the empirical density limit.
  - `lh_power_threshold(n20, B, S)` → P_LH (Martin 2008 / ITPA scaling) — the heating
    power needed to access H-mode.
  - `confinement_factor_lh(P, P_LH, h_factor)` — a smooth tanh **bifurcation**: τ_E
    multiplier ≈ 1 in L-mode, ≈ h_factor (~2) once heating exceeds threshold.
  - `confinement_factor_greenwald(n, n_G)` — a **collapse**: τ_E multiplier ≈ 1 below
    the limit, falling toward a floor as n → n_G (the edge-cooling density limit).
- `transport.burn_0d_ash` gained a `tau_factor(t, n_e, T, p_heat_density)` hook — a
  state-dependent multiplier on τ_E, applied on top of the β-limit. The two factors
  above plug straight in.
- `transport.Transport1D` gained a **soft β-limit** (`B`, `beta_limit`,
  `beta_stiffness`): above the Troyon β the heat diffusivity χ is raised so the
  volume-averaged pressure is pinned at the limit — the 1-D analogue of the 0-D cap
  in `burn_0d_ash`. This finally caps the otherwise-runaway 1-D burn (A2/A3 caveat).
- `animate.animate_operating_space` — animates several trajectories sweeping one (n,T)
  operating diagram together, with limit verticals and a shaded burning band.

## The three modes (the deliverable)

A small ITER-like toy device (R₀=3, a=1, B=5.3 T, Ip=7 MA → n_G = 2.2e20 m⁻³):

- **L-mode** — heating (≈7 MW) below the L→H threshold (≈16 MW): stays in low
  confinement, sits cool at T₀ ≈ 5 keV.
- **H-mode** — heating above threshold: confinement bifurcates up, the plasma ignites
  into the β-limited burning band, T₀ ≈ 22 keV, β ≈ 4.6 %.
- **Disruption** — over-fuel past the Greenwald limit: confinement collapses, the burn
  dies (T₀ → <1 keV at n/n_G ≈ 1.1). Reversible — back the fuelling off and it
  recovers to ≈22 keV.

## Validation (falsifiable — `tests/test_operating_limits.py` + `test_transport.py`)

- **n_G** = 1.19e20 m⁻³ for the ITER baseline (15 MA, a=2 m); scales ∝ Ip, ∝ 1/a².
- **P_LH** ≈ 52 MW at the ITER operating point (n20≈0.5, B=5.3, S≈680 m²; published
  ~50 MW); rises with density.
- The L→H factor bifurcates 1 → h_factor across the threshold; the Greenwald factor
  collapses 1 → floor across n_G.
- **H-mode runs >2× hotter than L-mode** across the L→H threshold (a real jump).
- The **over-fuel collapse past n_G is reversible** (back off the fuelling → recovers).
- The **1-D soft β-limit pins ⟨β⟩ ≈ 4 %** and lands a cooler core than the runaway.
  *(9 tests total: 8 operating-limits + 1 for the 1-D β-limit.)*

## Scope boundary (stated honestly)

- The L→H transition is modelled as a **prescribed confinement multiplier**, not a
  first-principles edge bifurcation; there is no real pedestal structure, ELM cycle,
  or L→H hysteresis (the back-transition uses the same threshold).
- The density-limit collapse is a **confinement-degradation model** of the Greenwald
  limit (edge cooling), not MARFE / detachment / impurity-radiation thermal-instability
  physics. The reversibility is a property of that model, chosen to match the
  qualitative experimental picture.
- These are **0-D** scenarios for the operating-window showcase; wiring the same
  τ_factor + 1-D β-limit into the F3 D-shaped burn is a natural follow-up.

## References
- Greenwald et al., *Nucl. Fusion* 28, 2199 (1988) — the density limit.
- Martin et al., *J. Phys. Conf. Ser.* 123, 012033 (2008) — the L→H threshold scaling.
- Troyon et al., *Plasma Phys. Control. Fusion* 26, 209 (1984) — the β-limit.
- Wesson, *Tokamaks* — H-mode, density limit, operating boundaries.

# 09 — Burning plasma: ignition → steady state → fuel injection

The whole life of a fusion discharge as a **transport** simulation: heat the
plasma, watch fusion α-particles take over the heating, ignite, settle into a
self-sustained burn, then inject fuel to drive it.

> **Why transport, not MHD/CFD?** This arc is an *energy/particle-balance* story on
> the confinement timescale (~seconds). MHD/CFD is the plasma *fluid moving* on the
> Alfvén timescale (~microseconds) — a different simulation you can't march across a
> whole discharge. The right tool for "ignition → steady → fuelling" is a transport
> code (TRANSP / ASTRA / RAPTOR); this is their laptop-scale toy cousin. The MHD/CFD
> side lives in experiment 06 and the T4 rung of the 3-D tokamak guide.

## What it shows

| Mode (`--mode`) | Fidelity | Result |
|------|---------|--------|
| `zerod` | F0 | 0-D Lawson/POPCON burn — two ODEs. The clearest view of **ignition**: α-heating overtakes losses, temperature runs away and settles at a thermally-stable burning point. The validation anchor. |
| `burn` (default) | F2 | 1-D radial transport. Three scripted phases — **ignition** (heating ramp), **steady burn** (heating off, Q→∞), **fuel injection** (a deep pellet pulse boosts density and fusion power). |

The physics terms (`plasmaplay.transport`): Bosch–Hale D-T reactivity `<σv>(T)`,
fusion α self-heating (3.5 MeV of the 17.6 MeV stays in the plasma), bremsstrahlung
radiation, prescribed-diffusivity transport, and pellet/gas-puff fuelling sources.

## Run it

```bash
python run.py                 # 1-D burn arc: ignition -> steady -> pellet (3 figures)
python run.py --mode zerod    # 0-D ignition / Lawson demo
python run.py --save          # write figures to ./outputs/
```

Outputs (`--save`): `burn_1d_traces.png` (T, n, power vs time with the three
phases shaded), `burn_1d_profiles.png` (radial profile snapshots),
`burn_1d_cross_sections.png` (the poloidal cross-section heating up and burning),
and `burn_0d_ignition.png` for the 0-D mode.

## What you should see

- **Ignition:** the core temperature climbs as auxiliary heating ramps; once hot
  enough, α-self-heating (∝ n²⟨σv⟩) overtakes the losses and the plasma stays lit
  when the heating is switched off — Q → ∞.
- **Steady burn:** profiles hold flat; volume-integrated α-power balances transport
  + radiation losses (the 0-D check passes to ~1%).
- **Fuel injection:** a pellet deposited inside the plasma raises the density, and
  because fusion power ∝ n² the fusion output jumps (≈ +30% here) — density control
  *is* a fusion-power knob.

## Honest caveats (what this model is not)

- Transport coefficients are **prescribed**, not predicted — χ is tuned so
  τ_E ≈ a²/(5.78 χ) ≈ 1 s. F3/F4 relax this.
- **No β-limit**, so the stable burning point sits hot (~80 keV); real machines are
  pressure-limited to ~15–25 keV. Single temperature, 50:50 D-T, no He-ash dilution
  yet (that's F1), circular geometry (real flux surfaces are F3).
- It is *not* MHD/CFD — no flows, no instabilities. That's the staged next half.

## Kernel exercised (validated in `tests/test_transport.py`, 16 tests)

`plasmaplay.transport`: `reactivity_dt` (Bosch-Hale, matches published `<σv>` to
~3%), `fusion_power_density`, `bremsstrahlung_density`, `burn_0d` (F0 ODE burn),
`Transport1D` (F2 implicit transport), `gaussian_deposition`.

## Concepts

Lawson criterion & the triple product nTτ_E; ignition vs. driven (Q) operation;
α-particle self-heating; thermal stability of the burning point; pellet vs.
gas-puff fuelling and deposition depth; why fusion power scales as n².

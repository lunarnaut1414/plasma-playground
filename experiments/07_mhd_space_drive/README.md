# 07 — An MHD drive for space propulsion

The repo's namesake goal: use the **J×B (Lorentz) force** — the same force that
drives the MHD dynamics of experiment 06 — to accelerate plasma and make thrust.
This is electric propulsion at its most powerful: magnetoplasmadynamic (MPD)
thrusters and applied-field accelerators reach thrust densities far beyond ion
engines, at specific impulses chemical rockets can't touch.

**Status:** F0 + F1 implemented (`run.py`). See [`PLAN.md`](PLAN.md).

## What it shows

| Case (`--case`) | Physics | Result |
|------|---------|--------|
| `scaling` (F0) | self-field MPD (Maecker) | thrust from the current's own field: **T ∝ I²** — why MPD wants megawatts |
| `channel` (F1) | applied-field accelerator | a prescribed Lorentz force accelerates the flow; **thrust = B₀ I L** (∝ I), exit speed sets Isp |

The headline is the **contrast in scaling**: self-field thrust grows as I² while
applied-field thrust grows as I — the two thrusters live in different regimes.

Sample numbers (F1, I = 2 kA, B₀ = 0.1 T, L = 0.2 m): exit speed 40 km/s, thrust
40 N, Isp ≈ 4080 s — and the ideal cold model conserves energy exactly (Lorentz
work = kinetic-energy gain, ratio 1.0000).

## Run it

```bash
python run.py --case scaling [--save]
python run.py --case channel
python run.py --case all --save
```

## Kernel exercised (validated in `tests/test_propulsion.py`)

`plasmaplay.propulsion` — `maecker_thrust` (T ∝ I²), `channel_velocity`
(impulse-momentum + energy balance hold exactly), `specific_impulse`, `jet_power`.

## Concepts / keywords

- J×B body force, magnetic pressure & tension (from experiment 06)
- Self-field vs applied-field MPD; Maecker's formula
- Thrust, exhaust velocity, specific impulse Isp, the thrust–power–Isp trade
- The rocket equation; electric propulsion vs chemical

## Caveats (what this ideal model leaves out)

The F1 channel is a *cold, lossless* model — hence 100% efficiency. Real MPD
thrusters run ~30–50% efficient: resistive dissipation, frozen-flow losses,
electrode falls, and onset instabilities all matter. Those enter at F2 (self-
consistent resistive MHD, reusing experiment 06's solver) and F4 (real propellant
+ comparison to measured thruster data). See [`PLAN.md`](PLAN.md).

## Next rung (F2 / F3)

F2: couple the flow to the field self-consistently (resistive MHD accelerator
channel, built on `plasmaplay.fvm`). F3: 2-D axisymmetric MHD of an applied-field
MPD or a **magnetic nozzle**, including the plasma-detachment problem.

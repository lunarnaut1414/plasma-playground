# 01 — Single-particle motion

The foundation of plasma physics: how one charged particle moves in electric
and magnetic fields. Master this and the rest of the playground (confinement,
drives, instabilities) is "just" many particles plus self-consistent fields.

## What it shows

| Case | Field setup | What you see | Why it matters |
|------|-------------|--------------|----------------|
| Gyro-orbit | uniform **B** | circle of radius `r_L` | basic magnetization — particles are tied to field lines |
| E×B drift | crossed **E**, **B** | gyration + sideways drift | drifts move plasma across field lines; charge-independent |
| Magnetic mirror | converging **B** | bouncing along z | the trapping principle behind mirror machines & tokamak banana orbits |

## Run it

```bash
python run.py          # show plots
python run.py --save   # also save a PNG to ./outputs/
```

## Concepts / keywords to look up

- Larmor (cyclotron) motion, gyrofrequency `ω_c = qB/m`, gyroradius `r_L`
- Guiding-center approximation and **drifts** (E×B, grad-B, curvature)
- Magnetic moment `μ = mv⊥²/2B` as an adiabatic invariant → mirror force
- **Boris algorithm** — why PIC codes use it instead of RK4 (energy conservation)

## The physics knob to play with

In `run.py`, try changing the launch pitch angle in `case_magnetic_mirror`
(ratio of `v_par` to `v_perp`). Particles inside the **loss cone**
(too much parallel velocity) escape instead of reflecting — that loss cone is
exactly why simple mirror machines leak and why tokamaks/stellarators were
invented.

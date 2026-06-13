# 05 — Stellarator field lines

How a stellarator confines plasma: make the magnetic field lines twist
(rotational transform ι) from the field's *shape* alone — no driven plasma
current. This experiment traces field lines and visualizes the flux surfaces and
the twist.

**Status:** F1 implemented (`run.py`). See [`PLAN.md`](PLAN.md) for the full ladder.

## What F1 shows

A straight **screw-pinch** field, `B = Bz ẑ + B_θ(r) θ̂` — the simplest field
with rotational transform — using the two core 3-D-confinement diagnostics:

| Diagnostic | What you see | Meaning |
|------------|--------------|---------|
| Poincaré section | nested closed curves (here, circles) | good flux surfaces — the plasma is confined on them |
| ι(r) profile | traced points on the analytic curve | magnetic **shear** (ι varies with radius), set by the B_θ(r) profile |

The traced ι matches the closed form `ι = (B_θ/r)/Bz · L/2π` to ~4 significant
figures at every radius.

## Run it

```bash
python run.py          # show plots
python run.py --save   # save a PNG to ./outputs/
```

## Kernels exercised (all validated in `tests/`)

- `plasmaplay.integrators.integrate_ode` / `rk4_step` — V3
- `plasmaplay.diagnostics.trace_field_line`, `poincare_section`, `rotational_transform` — V13
- `plasmaplay.fields.screw_pinch` — analytic field with a known ι

## Concepts / keywords

- Rotational transform ι (and safety factor q = 1/ι), magnetic shear
- Flux surfaces, Poincaré (puncture) plot, rational vs. irrational surfaces
- Screw pinch / "straight stellarator" as a teaching model

## Next rung (F2)

Replace the analytic field with real **coil filaments** via Biot–Savart
(`plasmaplay.fields.biot_savart`, already built and tested — V11) and look for
**magnetic islands** at rational surfaces. Then F3 brings in real stellarator
coils via `simsopt`. See [`PLAN.md`](PLAN.md).

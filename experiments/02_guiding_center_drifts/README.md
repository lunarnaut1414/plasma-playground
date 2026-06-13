# 02 — Guiding-center drifts

If a magnetized particle is just "a fast circle riding on a slowly drifting
center," can we drop the circle and follow only the center? Yes — and that
guiding-center picture is how confinement, transport, and the need for
rotational transform are actually reasoned about.

**Status:** F1 + F2 implemented (`run.py`). See [`PLAN.md`](PLAN.md).

## What it shows

| Case (`--case`) | Physics |
|------|---------|
| `drifts` (F1/F2) | the guiding center in a grad-B field, overlaid on the full Boris orbit — it is the gyro-averaged motion |
| `adiabatic` (F2) | sweep r_L/L and watch the GC error grow from ~3×10⁻⁴ (deeply adiabatic) to ~0.7 (r_L/L → 1) |

The drift catalogue (F0) lives in the kernel and its tests:
v_E = (E×B)/B² (E×B), v_∇B = (μ/q)(b̂×∇B)/B (grad-B), v_curv (curvature).

## Run it

```bash
python run.py --case drifts [--save]
python run.py --case adiabatic
python run.py --case all --save
```

## Kernel exercised (validated in `tests/test_guiding_center.py`)

`plasmaplay.guiding_center` — `gc_push` integrates the guiding-center equations
of motion (drifts + parallel mirror force) with numerical field gradients.
Validated: E×B drift exact, grad-B drift vs analytic, and GC trajectory matches
the gyro-averaged Boris orbit.

## Concepts / keywords

- Guiding-center approximation, the adiabaticity parameter r_L/L
- E×B / grad-B / curvature drifts; the magnetic moment μ as an adiabatic invariant
- Connection to neoclassical transport and banana orbits

## Next rung (F3)

Drifts in a real tokamak field → trapped-particle **banana orbits** and why the
poloidal field's rotational transform is needed to cancel the vertical grad-B
drift. Then F4: an ensemble + collisions → neoclassical transport regimes. See
[`PLAN.md`](PLAN.md). The guiding center can also run in your *own* experiment-04
equilibrium field.

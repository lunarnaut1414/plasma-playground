# 06 — Ideal MHD: waves and shocks

Treat the plasma as a single conducting fluid threaded by magnetic field and
solve the ideal magnetohydrodynamics (MHD) conservation laws. MHD governs the
*macroscopic* life of a plasma — stability of fusion devices, the solar wind,
and (experiment 07) the J×B thrust of an MHD drive.

**Status:** F1 + F2 implemented (`run.py`). See [`PLAN.md`](PLAN.md).

## What it shows

| Case (`--case`) | Physics | Result |
|------|---------|--------|
| `alfven` (F1) | circularly-polarized Alfvén wave | translates at **v_A = Bx/√ρ** with shape preserved — magnetic tension as a restoring force |
| `briowu` (F2) | Brio–Wu MHD shock tube | the full nonlinear structure: fast rarefactions, the **compound wave**, a contact, and a slow shock |

The solver (`plasmaplay.fvm`) is a finite-volume Godunov scheme: MUSCL (minmod)
reconstruction → HLL Riemann flux → SSP-RK2 in time. MHD has *three* wave speeds
(slow, Alfvén, fast); the printout reports them for the Brio–Wu left state.

## Run it

```bash
python run.py --case briowu [--save]
python run.py --case alfven
python run.py --case all --save
```

## Kernel exercised (validated in `tests/test_fvm.py`, V8)

`plasmaplay.fvm.solve_mhd_1d` — checked against the exact Sod hydro star values
(p* = 0.30313), an exact Alfvén wave (speed v_A, amplitude preserved), and
Brio–Wu (positivity, conservation, By sign reversal).

## Concepts / keywords

- Ideal MHD conservation laws, the J×B force = magnetic pressure + tension
- Alfvén / fast / slow magnetosonic waves; the Friedrichs diagram
- Godunov / finite-volume methods, HLL Riemann solver, MUSCL reconstruction, CFL
- Brio–Wu shock tube, compound waves

## Next rung (F3)

2-D ideal MHD: the **Orszag–Tang vortex**, which forces the ∇·B = 0 problem and
constrained-transport (the one thing that makes MHD codes special). Then F4 adds
resistivity → magnetic reconnection. See [`PLAN.md`](PLAN.md). This same solver
is the engine of experiment 07 (the MHD space drive).

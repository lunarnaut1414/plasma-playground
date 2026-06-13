# 04 — Tokamak equilibrium (Grad–Shafranov)

What shape does a magnetically confined plasma settle into? In an axisymmetric
torus, the force balance ∇p = J×B reduces to one elliptic PDE for the poloidal
flux ψ — the **Grad–Shafranov equation** — whose level sets are the magnetic
flux surfaces. This is *the* equation of tokamak design.

**Status:** F1 implemented (`run.py`). See [`PLAN.md`](PLAN.md).

## What F1 shows

A **Solov'ev** equilibrium (linear pressure/current profiles, so the solve is a
single linear elliptic problem) on a fixed rectangular boundary:

- nested closed **flux surfaces** (contours of ψ),
- the **magnetic axis** (the ψ extremum), and
- the **Shafranov shift** — the axis sits *outboard* of the geometric center
  (measured here: +0.118 m, ~26% of the minor radius), pushed there by the
  toroidal 1/R term in Δ*. Real tokamaks show exactly this.

## Run it

```bash
python run.py [--save]
```

## Kernels exercised (validated in `tests/`)

- `plasmaplay.solvers.grad_shafranov_solve` — the Δ* fixed-boundary solver,
  validated to 2nd-order against a manufactured exact solution (**V12**).

## Concepts / keywords

- Grad–Shafranov equation, poloidal flux ψ, the Δ* (toroidal elliptic) operator
- Flux surfaces, magnetic axis, separatrix
- Shafranov shift, Solov'ev equilibrium

## The physics knobs

In `run.py`, the source is `-(c_p R² + c_0)` (Solov'ev linear profiles):
`c_p` is the pressure-gradient term, `c_0` the FF′ (poloidal-current) term.
Change their ratio and watch the axis and surface shapes respond — more pressure
→ larger Shafranov shift.

## Next rung (F2)

Free-boundary equilibria with real poloidal-field **coils** via
[FreeGS](https://github.com/freegs-plasma/freegs): the plasma boundary, X-point,
and divertor *emerge* from the coil currents instead of being prescribed. Then
F3 matches a real machine (MAST-U / DIII-D / ITER) and extracts the q-profile and
β. See [`PLAN.md`](PLAN.md). The equilibrium B-field here can also become the F3
field for experiments 01/02 (banana orbits in a *real* equilibrium).

# E1 — The stellarator: rotational transform without plasma current

The stretch rung (NIGHT.md Track E). Everything up to here has been a **tokamak** —
its field-line twist (rotational transform ι) comes from a driven toroidal plasma
current, which is also what makes it prone to disruptions and sawteeth. A
**stellarator** gets ι from the 3-D shape of its *external* coils instead, with
essentially **zero net plasma current** — so it is inherently steady-state and has no
current-driven kinks, no sawteeth, no disruptions. E1 builds a genuine stellarator
vacuum field and shows exactly that contrast.

## Why experiment 05's existing field was not (yet) a stellarator

Experiment 05's F1 used `fields.screw_pinch`: B = B_z ẑ + B_θ(r) θ̂. Its twist comes
from B_θ — i.e. from an **axial current**. That is the *tokamak* mechanism dressed up;
it is current-driven. A true stellarator must produce ι with **no net current**.

## The field (`fields.helical_stellarator`)

The straight-stellarator limit: a strong axial guide field plus a single helical
harmonic, written as the gradient of a **harmonic scalar potential** so the field is
**curl-free (current-free) in the plasma region by construction**:

    Phi = B0 z + eps I_l(h r) cos(l θ − h z),     ∇²Phi = 0   (modified Bessel I_l)

    B_r     =  eps h I_l'(h r) cos(l θ − h z)
    B_θ     = −eps (l/r) I_l(h r) sin(l θ − h z)
    B_z     =  B0 + eps h I_l(h r) sin(l θ − h z)

`l` is the field periodicity (l = 2 is the classic), `h = 2π/L` the axial wavenumber.
A single helicity averages to zero twist at *first* order; ι appears at **second
order, ι ∝ eps²** — the current-free geometric transform that defines a stellarator.

## Validation (falsifiable — `tests/test_fields.py`, 2 tests)

- **Current-free:** the numerical curl ∇×B ≈ 0 (< 1e-4) everywhere in the plasma
  region, and the loop integral ∮B·dl around the magnetic axis vanishes (< 1e-9) — so
  there is **no net plasma current**, the defining stellarator property.
- **Transform from geometry:** a current-free field still twists its lines (ι ≠ 0,
  traced from the actual field), and ι **grows with the helical shaping amplitude**
  (ι(0.6) > ι(0.3)) — the 2nd-order ι ∝ eps² origin. (In `run.py --mode stellarator`,
  at higher tracing resolution: ι = 0.015, 0.021, 0.039 for eps = 0.3, 0.5, 0.7.)

The Poincaré section (run.py) shows the **nested flux surfaces**.

## Deliverable

`outputs/stellarator_flux_surfaces.gif` — field lines traced on several nested flux
surfaces of the l=2 vacuum field, with the straight-stellarator cylinder mapped onto a
torus (its large-aspect limit) so the elliptical cross-section rotates around it — the
iconic twisty stellarator, rotating. Regen: `python gif_gallery.py
stellarator_flux_surfaces`; `python experiments/05_stellarator_field_lines/run.py
--mode stellarator --save` for the Poincaré + 3-D still.

## Scope boundary (stated honestly)

- **Straight-stellarator / single-helicity** model: a vacuum field from one helical
  harmonic, mapped onto a torus only for *rendering* (the physics is the straight
  large-aspect limit). Not a real optimized 3-D equilibrium (W7-X / NCSX use many
  harmonics + finite β + a Grad–Shafranov-like 3-D solve — VMEC territory).
- ι is small (~0.02–0.04) and the absolute value from field-line tracing carries a
  finite-step floor (so the tests assert curl-free + ι≠0 + monotone-in-eps, not a
  precise ι∝eps² coefficient).
- The field is a **vacuum** field (no pressure, no self-consistent currents); E2 would
  run the Track-A transport on these flux surfaces (inherently steady-state — no
  sawteeth/disruptions, the explicit tokamak contrast).

## E2 — Transport on the stellarator: the steady-state contrast

The payoff of building a genuine current-free field: the **operating modes collapse**.
A tokamak needs a driven plasma current for both confinement and ι — and that current
is what makes q(0) fall below 1 as the core heats (the sawtooth, Track C) and what the
Greenwald limit caps (the disruption, F3.5). A stellarator has **no net current**, so:

- its q-profile is fixed by the coils (`sawtooth.external_q_profile`, q > 1 everywhere);
- the burning core **cannot drive q(0) < 1** → no m=1 kink → **no sawteeth, no
  current-driven disruptions**;
- the discharge is just **startup → steady burn → fuelling**, inherently steady-state.

Validation (`tests/test_sawtooth.py::test_stellarator_burn_is_sawtooth_free`): the same
hot, peaked core that fires a sawtooth in a current-driven tokamak (q_from_temperature
crosses 1) is **sawtooth-free** on the stellarator — its `external_q_profile` has no
q = 1 surface, so `mixing_radius` returns None (no reconnection possible, ever).

Deliverable: `outputs/stellarator_burn.gif` — the same exp-09 transport burn run on the
stellarator (`gif_gallery stellarator_burn`): **0 sawtooth crashes vs the tokamak's
179**, the smooth steady stellarator core-T0 shown beside the tokamak's sawtoothing
trace. The honest reduction: same transport physics, the device topology removes the
MHD-event coupling — it does *not* add stellarator-specific transport (neoclassical /
turbulent χ from 3-D geometry, the real stellarator confinement story, is beyond F2).

## References
- Spitzer, *Phys. Fluids* 1, 253 (1958) — the stellarator concept.
- Freidberg, *Ideal MHD* — the straight-stellarator helical vacuum field.
- Boozer, *Rev. Mod. Phys.* 76, 1071 (2004) — stellarator transform & symmetry.
- Dinklage et al., *Nat. Phys.* 14, 855 (2018) — W7-X: steady-state, current-free operation.

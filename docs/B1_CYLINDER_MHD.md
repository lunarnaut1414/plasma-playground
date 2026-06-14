# B1 — Cylindrical linear MHD: kink & tearing on a real q(r)

The first rung of the MHD-instability track (NIGHT.md Track B) and the start of
experiment 10. It lifts the slab tearing mode of T4 (`plasmaplay/tearing.py`) onto a
**periodic cylinder** — the "straight tokamak" — where a current profile gives a
safety factor `q(r)` and the two organising tokamak instabilities finally appear:
the **m=1 internal kink** (the sawtooth trigger) and the **m/n tearing mode**.

## The straight tokamak

A torus cut and straightened: radius `a`, length `2πR`, strong axial guide field
`B_z`. A peaked current `J_z(r) = J0 (1 − (r/a)²)^ν` gives, by Ampère's law, a closed
form for the safety factor:

    q(r) = q0 (ν+1) r² / (1 − (1 − r²)^(ν+1)),     q(0)=q0,  q(a)=(ν+1)q0.

A perturbation `~ exp(i(mθ + kz))`, `k = −n/R`, resonates on the **rational surface**
`q(r_s) = m/n`. Because the instability drive enters the outer equation only through
the ratio `μ0 J_z' / (B_θ (1 − nq/m))`, the overall `B_θ` scale (i.e. `R`) cancels —
**stability depends only on the q-profile and the mode numbers** (we set a=B_z=R=1).

## What was built (`plasmaplay/cylinder_mhd.py`)

- `screw_pinch_q`, `b_theta`, `rational_surface` — the equilibrium and the q=m/n
  surface (bisection on the monotone q).
- `delta_prime_cylinder(m, n, q0, ν)` — integrates the outer **Newcomb equation**

      ψ'' + ψ'/r − (m²/r²)ψ − [μ0 J_z'(r)/(B_θ(r)(1−nq/m))] ψ = 0

  inward from the axis (regular `ψ~r^m`) and from the wall (`ψ(a)=0`) to the rational
  surface, returning **Δ′ = [ψ'/ψ]₊ − [ψ'/ψ]₋**. Δ′ > 0 ⇒ tearing-unstable.
- `internal_kink_unstable(q0)` — the m=1/n=1 criterion: a q=1 surface exists ⟺
  q(0) < 1. `internal_kink_xi(r, r1)` — the ideal m=1 eigenfunction (a rigid top-hat
  core displacement falling to zero at the q=1 surface).
- `fkr_growth_rate(Δ′, S)` — the Furth–Killeen–Rosenbluth resistive-layer law
  `γτ_A = 0.55 Δ′^{4/5} S^{−3/5}`. The layer physics is local, so it carries over
  unchanged from the slab; only Δ′ becomes the cylindrical value.

## Validation (falsifiable — `tests/test_cylinder_mhd.py`, 14 tests)

- q(0)=q0, q(a)=(ν+1)q0, q monotone; the q=m/n surface sits where q actually = m/n.
- **The m=1 internal kink is unstable exactly when q(0) < 1** (tested at q0 = 0.7,
  0.85, 0.95 → unstable; 1.05, 1.3 → stable) — the sawtooth trigger.
- **The sign of Δ′ predicts tearing stability** (a broad-current m=2 surface → Δ′>0;
  a high-m edge surface → Δ′<0), and **Δ′ falls with mode number m** (higher m is
  more wall-stabilized) — both robust to the layer-skip gap.
- **γ ∝ S^(−3/5)**: a decade in Lundquist number drops γ by 10^(−0.6) (to 1e-6); γ=0
  for a tearing-stable (Δ′≤0) surface.

## Deliverable

`outputs/kink_eigenmode.gif` — the m=1 internal kink: the eigenfunction `ξ_r(r)` and
`q(r)` beside the poloidal cross-section, the hot core displaced into the
characteristic crescent and growing. Regen: `python gif_gallery.py kink_eigenmode`;
`python experiments/10_tokamak_stability/run.py --save` for the still.

## Scope boundary (stated honestly)

- **The absolute Δ′ is resolution-dependent** near the singular layer (the outer
  solution carries a logarithmic term that cancels only in the gap→0 limit); the
  **sign and the m-ordering** — the physically meaningful, charter-mandated outputs —
  are gap-robust, and the tests assert only those, not an absolute Δ′ value.
- The growth rate is the **constant-ψ FKR layer estimate**, not a full cylindrical
  resistive eigenvalue solve (that was done for the slab in T4, which pins the
  S^(−3/5) slope by eigenvalue, not formula). The cylindrical solve is a possible
  refinement.
- Cylinder only: **no toroidal coupling** of poloidal harmonics, no shaping, no
  ballooning modes — those are the toroidal rungs beyond B-track.

## References
- Newcomb, *Ann. Phys.* 10, 232 (1960) — the cylindrical stability equation.
- Furth, Killeen & Rosenbluth, *Phys. Fluids* 6, 459 (1963) — tearing modes & Δ′.
- Wesson, *Tokamaks* — the kink/tearing/sawtooth picture and the q(0)<1 trigger.

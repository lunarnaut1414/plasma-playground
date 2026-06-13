# Building toward a 3-D tokamak

A practical guide to extending the playground into a **3-D tokamak** — what is
genuinely reachable on a laptop, in what order, reusing what you've already
built, with a falsifiable check at every step.

## Scope: read this first

"Full 3-D tokamak simulation" spans wildly different problems. Be clear which one
you mean:

| Target | Real codes | In scope here? |
|--------|-----------|----------------|
| Magnetic field + single-particle orbits in a real equilibrium | (analysis tools) | ✅ **yes** — this guide |
| Field-line topology, islands, stochasticity (3-D fields) | POINCARÉ tools, SIESTA | ✅ **yes** — this guide |
| 3-D MHD stability (kinks, tearing, ELMs, disruptions) | JOREK, NIMROD, M3D-C1 | ◐ a *reduced* taste (T4); real version = research code |
| Turbulence & transport (gyrokinetics) | GENE, XGC, CGYRO, stella | ❌ no — supercomputer research codes |
| Integrated whole-device | TRANSP, IMAS | ❌ no |

The reachable goal: **"I built a tokamak's 3-D magnetic field, traced banana
orbits and the q-profile in it, broke axisymmetry to make magnetic islands, and
watched flux surfaces go stochastic."** That is a real, validated, laptop-scale
3-D tokamak playground. Self-consistent turbulence/disruption physics is a
different category of software (teams, decades, 10³–10⁵ cores) — for that, the
answer is to learn one of the codes above.

> A tokamak *equilibrium* is **axisymmetric (2-D)** by symmetry. "3-D" means
> breaking that symmetry — instabilities, coil ripple, applied perturbations,
> turbulence. T0–T2 below build the (2-D-symmetric) field and orbits in 3-D
> space; T3 is where real 3-D field structure begins.

## The key idea — you're one bridge away

Experiment 04 already solved for the poloidal flux ψ(R, Z). The **full magnetic
field** of an axisymmetric tokamak follows directly:

```
B = ∇φ × ∇ψ + F(ψ) ∇φ        (∇φ = φ̂/R)

  B_R = -(1/R) ∂ψ/∂Z
  B_Z =  (1/R) ∂ψ/∂R
  B_φ =  F(ψ)/R               (F = R B_φ, the toroidal-field function)
```

So a ~30-line bridge — differentiate the experiment-04 ψ, add the toroidal field,
wrap it as a callable `B(x)` in Cartesian (X, Y, Z) — turns the 2-D equilibrium
into a **3-D-in-space field that drops straight into the pushers and tracers you
already have**. Everything downstream (T1, T2) is then mostly *reuse*.

---

## The ladder

Each rung: the goal, the physics, what to build (and which existing kernel it
reuses), the validation check, and rough compute.

### T0 — The 3-D equilibrium field (the bridge)
- **Goal:** a callable `tokamak_field(x) -> (3,) B` built from an experiment-04 equilibrium.
- **Build:** new `plasmaplay/tokamak.py`:
  - `equilibrium_field(R, Z, psi, F_of_psi)` — finite-difference B_R, B_Z; B_φ = F/R; bilinear-interpolate onto arbitrary (R, Z); convert (X,Y,Z) ↔ (R, φ, Z).
  - Choose `F(ψ)`: simplest is the vacuum value `F = R₀ B₀` (constant); a Solov'ev FF′ = const gives `F(ψ) = √(F₀² + 2·FF′·(ψ−ψ_b))`.
- **Reuses:** `solvers.grad_shafranov_solve` (experiment 04) for ψ.
- **Validation:** **∇·B ≈ 0** numerically everywhere; on the magnetic axis B is purely toroidal; |B| ∝ 1/R to leading order (the 1/R tokamak field). 
- **Compute:** instant (interpolation).

### T1 — Field-line topology & the q-profile
- **Goal:** trace field lines, make a Poincaré section at φ = const, extract the safety factor q(ψ).
- **Build:** almost nothing new — point the existing tracer at the T0 field.
- **Reuses:** `diagnostics.trace_field_line`, `poincare_section`, `rotational_transform` (experiment 05, **as-is**). q = 1/ι = toroidal turns per poloidal turn.
- **Validation:** field lines lie on nested closed surfaces (Poincaré = nested curves); q increases outward; for a large-aspect-ratio circular equilibrium q(r) ≈ r B_φ /(R B_θ) matches the analytic estimate within a few %.
- **Compute:** seconds–minutes.

### T2 — Particle orbits: banana orbits & drift surfaces
- **Goal:** trapped "banana" orbits and passing orbits in the real field; the trapped/passing boundary; the outward grad-B drift and why the poloidal field cancels it.
- **Build:** nothing new — launch particles into the T0 field.
- **Reuses:** `pushers.boris_push` and `guiding_center.gc_push` (experiments 01/02, **as-is**). Project orbits onto the (R, Z) plane.
- **Physics:** |B| is larger on the inboard side (small R), so particles with small v∥/v mirror-**trap** and bounce, tracing a banana in the poloidal plane; the guiding center drifts radially off the flux surface by the banana width.
- **Validation:** trapped fraction ≈ √(2ε) (ε = r/R₀); banana width scales like q·ρ/√ε; passing particles stay (nearly) on a flux surface. μ conserved along the orbit.
- **Compute:** seconds per orbit (Boris); guiding-center is cheaper for long times.

### T3 — Break axisymmetry: magnetic islands & stochasticity (real 3-D)
- **Goal:** add a non-axisymmetric perturbation and watch flux surfaces tear into **magnetic islands** at rational (q = m/n) surfaces, then go **stochastic** when islands overlap.
- **Build:** add a resonant helical perturbation to the T0 field, e.g. δB from
  `δψ = δ · cos(mθ − nφ)` (a divergence-free perturbation; simplest is a small
  radial field `B_r ∝ δ cos(mθ − nφ)`). Then re-run the T1 Poincaré.
- **Reuses:** the T1 Poincaré machinery, now showing island chains.
- **Physics:** islands form only at the **resonant surface** q = m/n; island width W ∝ √(δ); when neighboring island chains overlap (Chirikov criterion), field lines become chaotic — the loss of good surfaces that limits tokamak performance.
- **Validation:** island appears at the predicted q = m/n radius; W ∝ √(perturbation amplitude); increasing δ drives the Poincaré plot from islands → stochastic sea. This is the first genuinely **3-D** field structure.
- **Compute:** minutes (many field-line punctures).

### T4 — (stretch) Reduced 3-D MHD instabilities
- **Goal:** evolve a simplified MHD model in the torus to see a **tearing mode** or **internal kink** grow — a self-consistent 3-D instability, not a prescribed perturbation.
- **Build:** a reduced-MHD solver (Strauss reduced MHD: evolve poloidal flux ψ and vorticity/stream function U, two coupled scalar PDEs, in (r, θ, φ) with the toroidal mode n as a Fourier index). This is a real PDE project — finite differences in r, spectral in θ/φ, implicit or small-timestep explicit, with resistivity for tearing.
- **Reuses:** the elliptic solver (`solvers`) for the stream-function inversion; FFT diagnostics for mode growth.
- **Validation:** linear tearing-mode growth rate scales like the resistive-MHD prediction γ ∝ S^(−3/5) (S = Lundquist number); the m/n mode is unstable only when Δ′ > 0. Match a published linear growth rate.
- **Compute:** minutes–hours. **This is the hard rung** — expect it to be a project, and a simplified one compared to JOREK/NIMROD.

### Beyond T4 — where the laptop stops
Self-consistent **3-D extended MHD** (disruptions, ELMs) and **gyrokinetic
turbulence/transport** are the domain of established codes. They are open source
and worth learning if this is your direction:
- MHD: **JOREK**, **NIMROD**, **M3D-C1**
- Gyrokinetics: **GENE**, **XGC**, **CGYRO**, **stella**
- Equilibrium (a better T0): **FreeGS** (you already use it on the experiment-04 roadmap), **EFIT**

Driving one of those on a real case is the honest path to "a full 3-D tokamak
simulation."

---

## New kernels this adds (small)

Most of the ladder is *reuse*. The genuinely new code is modest:
- `plasmaplay/tokamak.py` — the equilibrium-field bridge (T0) + the resonant
  perturbation (T3).
- a poloidal-projection / banana-width diagnostic (T2) and a q-profile helper (T1).
- (only at T4) a reduced-MHD solver — a real new component.

Each ships with its validation check, same as every other kernel in
[`FUNDAMENTALS.md`](FUNDAMENTALS.md).

## Suggested order & first step

T0 → T1 → T2 is a tight, high-payoff sequence that mostly wires existing kernels
together and gets you orbits + q-profile + Poincaré in a *real* equilibrium —
genuinely "a tokamak in 3-D space." T3 adds the first true 3-D field structure.
T4 is an open-ended stretch.

**Start at T0:** build `equilibrium_field` from the experiment-04 ψ, validate
∇·B ≈ 0 and the 1/R field, then immediately reuse the experiment-05 Poincaré and
experiment-01/02 pushers on it. That one bridge unlocks T1 and T2 at once.

# T1 — field-line topology: the q-profile & toroidal Poincaré section

Rung **T1** of [`3D_TOKAMAK_GUIDE.md`](3D_TOKAMAK_GUIDE.md). Builds on the T0
field (`plasmaplay/tokamak.py::equilibrium_field`) and adds the two diagnostics
that define a tokamak's magnetic topology:

- `toroidal_poincare(B_func, start_RZ, n_punctures, ds)` — punctures of a field
  line through the φ = 0 half-plane.
- `safety_factor(B_func, start_RZ, axis_RZ, n_poloidal, ds)` — q = toroidal
  turns / poloidal turns.

Visualised by `tokamak_t1_viz.py` (→ `outputs/tokamak_t1.png`); validated by the
`# T1` tests in `tests/test_tokamak.py`.

---

## Design choices (and why)

### 1. Section in φ, not in a Cartesian plane
The stellarator diagnostic `diagnostics.poincare_section` sections on a fixed
Cartesian *z*-plane. That is wrong for a tokamak: a tokamak field line winds in
the toroidal angle φ while only *oscillating* in z, so it would cross z = const
twice per poloidal turn at unrelated φ — not a clean section. The natural section
is the **half-plane φ = 0** (one puncture per completed toroidal turn). Hence a
new `toroidal_poincare` rather than reusing the Cartesian one. On a good flux
surface its punctures trace the surface's poloidal cross-section.

### 2. Angles accumulated modulo 2π
Both `toroidal_poincare` and `safety_factor` accumulate φ (about z) and θ (about
the magnetic axis) by **summing wrapped increments** `_wrapped(Δangle) =
angle(exp(iΔ))`, never by differencing raw `atan2` values. This is robust across
the atan2 branch cut: the running totals grow smoothly past ±π instead of
jumping by 2π. q is then simply `Δφ_total / Δθ_total` once the line has wound
`n_poloidal` times poloidally.

### 3. Arc-length tracing with RK4
Field lines are integrated as `dx/ds = B/|B|` (unit speed) with the shared
`integrators.rk4_step`. Unit speed makes `ds` a geometric step (a length), so the
same `ds` behaves consistently regardless of |B|.

---

## Performance: one vector interpolator, not three

**Change:** `equilibrium_field` originally built three scalar
`RegularGridInterpolator`s (one each for B_R, B_Z, B_φ) and called all three per
field evaluation. The per-call Python overhead of `RegularGridInterpolator`
dominates field-line tracing (millions of evals). Stacking the three components
into one vector-valued interpolator (`values.shape = (nR, nZ, 3)`, one call
returning `(B_R, B_Z, B_φ)`) cut the work ~3× **everywhere** — pushers, tracers,
Poincaré — for no loss of accuracy (identical interpolant).

**Gotcha noted while doing it:** a vector-valued `RegularGridInterpolator` called
with a single point tuple `interp((r, z))` returns shape **`(3,)`** directly, not
`(1, 3)`. Indexing `[0]` (as one would for a scalar interpolator over a batch)
silently breaks. Tests caught it immediately.

Combined with test-parameter tuning (below), the tokamak suite went **113 s →
~27 s**.

---

## A real bug the visualization exposed: rational-surface degeneracy

The first T1 benchmark equilibrium was a large-aspect circular one with
ψ = C·r², which has a **uniform** safety factor q = a·B₀/(Rc·Bθ0) = 2 on *every*
surface. q = 2 is **rational** (= 2/1), so every field line **closes on itself**
after 2 toroidal turns and punctures the φ = 0 plane at only **2 discrete
points** (θ = 0 and π). The picture showed dots/arcs, not filled circles.

Worse, this made `test_poincare_is_a_closed_nested_surface` **pass trivially**:
its check is "all punctures lie at one minor radius" (`std/mean < 1e-3`), and two
coincident points satisfy that without ever exercising "the line fills a closed
curve." A green test that wasn't testing the thing.

**Fix.** Use a **sheared** circular equilibrium for both the picture and the
test, where q varies with radius and is therefore irrational on essentially every
surface, so a single traced line **densely fills** its circle:

```
choose  q(r) ≈ q₀ (1 + (r/a)²)
   =>   dψ/dr = R Bθ ≈ B₀ r / [q₀ (1 + (r/a)²)]      (large aspect, B_φ ≈ B₀)
   =>   ψ(r) = (B₀ a² / 2q₀) · ln(1 + (r/a)²)
```

The test now also asserts the punctures **spread around ≈2π in poloidal angle**
(`ptp(poloidal) > 5.0`), which is the property that was silently untested. The
uniform-q = 2 equilibrium is *kept* — it's the right tool for the *analytic q*
benchmark (`test_safety_factor_matches_large_aspect_analytic`), where a closed
rational orbit is fine because `safety_factor` measures q correctly regardless.

**Lesson:** rational vs irrational q is not a detail — it determines whether a
Poincaré section is points or a curve. Benchmark equilibria for *topology* tests
should be sheared (irrational); benchmark equilibria for *scalar-q* tests can be
uniform.

This sheared profile is also a gift to T3: its q crosses **q = 1, 3/2, 2**, which
are exactly the rational surfaces where a resonant perturbation will grow
magnetic islands.

---

## Validation (all in `tests/test_tokamak.py`)

| Check | Assertion | Result |
|-------|-----------|--------|
| q matches analytic large-aspect estimate | `|q−2|/2 < 3%` at r = 0.4, 0.7, 1.0 | q ≈ 2.005, err ≈ 0.5% |
| q is a converged number | `|q(4 turns) − q(10 turns)| < 2%` | passes |
| Poincaré = densely-filled closed nested circle | `std/mean < 2e-3` **and** poloidal spread > 5 rad; inner surface nested inside outer | passes |
| Real Solov'ev equilibrium has shear | q increases monotonically outward, q(axis) > 0.5 | q rises 1.78 → 3.19 |

### `ds` choice in the tests
For these smooth analytic equilibria the field lines lie on exact circles, so RK4
arc-length tracing is converged well before `ds` matters: q at ds = 0.02, 0.05,
0.08 is identical to 5 digits (err 0.005, the geometric O(ε²) term, not a stepping
error). The tests use **ds = 0.05–0.06**, ~2.5× faster than 0.02 with no accuracy
loss. Don't mistake this for a general license to use coarse `ds`: on a *stochastic*
field (T3) the trajectory is sensitive and a finer step is warranted.

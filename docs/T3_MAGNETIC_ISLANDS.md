# T3 — breaking axisymmetry: magnetic islands & stochasticity

Rung **T3** of [`3D_TOKAMAK_GUIDE.md`](3D_TOKAMAK_GUIDE.md), and the first rung
with genuinely **3-D field structure**. T0–T2 lived in a 2-D-symmetric
(axisymmetric) equilibrium drawn in 3-D space; T3 *breaks* that symmetry with a
non-axisymmetric perturbation and watches the nested flux surfaces tear.

New kernel: `plasmaplay/tokamak.py::helical_perturbation` (+ `superpose`).
Visualised by `tokamak_t3_viz.py` (→ `outputs/tokamak_t3.png`); validated by the
`# T3` tests in `tests/test_tokamak.py`.

---

## The central design choice: perturb the *flux*, stay divergence-free

The guide suggests "the simplest is a small radial field B_r ∝ δ cos(mθ − nφ)".
A bare radial field like that is **not divergence-free** — and a non-solenoidal B
poisons field-line tracing (the lines aren't tangent to a real field). So instead
of perturbing B directly, perturb the **poloidal flux**, exactly as the
equilibrium field is built:

```
Ψ(R,Z,φ) = ψ(R,Z) + δψ_h,     δψ_h = δ · env(r) · cos(mθ − nφ)
B = ∇Ψ × ∇φ + F ∇φ
```

This is divergence-free **identically**, because `∇·(∇Ψ × ∇φ) ≡ 0` for *any*
scalar field Ψ — even one that depends on φ. (It's a `∇·(∇A×∇B)=0` identity.)
Keeping F = R₀B₀ constant means the toroidal term `F∇φ` is also divergence-free.

Working out `∇Ψ × ∇φ` in cylindrical components gives a small, clean surprise:

```
B_R = -(1/R) ∂Ψ/∂Z,   B_Z = (1/R) ∂Ψ/∂R,   and the ∂Ψ/∂φ piece contributes 0
```

— the φ-derivative lands along `φ̂ × φ̂ = 0`, so the helical perturbation enters
**only** through the R and Z derivatives of δψ_h (evaluated at the local φ), with
**δB_φ = 0**. The implementation just finite-differences δψ_h in R and Z. The
`test_perturbation_stays_divergence_free` check confirms ∇·B/|B| < 10⁻³ stays put
with the perturbation on (numerically ≈ 4×10⁻⁶).

`env(r)` is a radial envelope (default constant). The viz/tests pass a **Gaussian
centred on the resonant surface**, for two reasons: it localises a single island
chain so chains don't blur together, and it avoids the `∂θ/∂R ∼ 1/r` blow-up of
the poloidal angle near the magnetic axis.

---

## Why a single mode can't be stochastic

A *single* helical mode m/n leaves a conserved helical flux (the system is still
integrable — one ignorable coordinate). So one mode, at *any* amplitude, makes a
bigger island, never chaos. **Stochasticity requires ≥ 2 resonances** whose island
separatrices overlap (the Chirikov criterion). The viz therefore drives **two**
modes (3/2 and 2/1) on a sheared equilibrium tuned (q₀ = 1.2) so those two
rational surfaces sit at mid-radius, and raises the amplitude until the islands
between r(q=3/2) and r(q=2) overlap into a stochastic layer.

This is a real physics point worth stating plainly: the chaos is not in any one
resonance — it lives in the *overlap* of neighbouring ones.

---

## Measuring island width without fooling yourself

Island width is measured as the radial spread (`max−min` minor radius) of the
puncture pattern of a field line started at the chain's **X-point**, which traces
the separatrix and so spans the full island.

A pitfall found while building this: a line started **on the resonant surface at
an O-point poloidal angle** sits near the *elliptic* fixed point, and its nearby
contour grows **∝ δ (linear)**, not ∝ √δ — measuring there gave a misleading
~linear scaling. Starting at the X-point (the *hyperbolic* fixed point, here the
top of the chain at θ = π/2) recovers the textbook **W ∝ √δ**.

A second pitfall: the resonance must be at the surface you think it is. At aspect
ratio 5 the large-aspect estimate q ≈ q₀(1+(r/a)²) is off by ~10%, so the q=2
surface had to be placed by **measuring** q(r), not by the formula. With q₀ = 1.6
the measured q(0.5) = 2.010 — then the scaling came out clean.

### Validated scaling (Rc=5, q₀=1.6, q=2 at r=0.5, Gaussian envelope)
| δ | island width W | ratio per 2× δ |
|------|------|------|
| 2.5e-4 | 0.141 | — |
| 5.0e-4 | 0.199 | ×1.41 |
| 1.0e-3 | 0.282 | ×1.42 |
| 2.0e-3 | 0.403 | ×1.43 |

`√2 = 1.414` — **W ∝ √δ** to ~1 %.

---

## Performance note — aspect ratio drives the cost

Field-line tracing cost scales with the **toroidal circumference** 2πR₀: every
puncture is one full toroidal turn of arc length. The T1/T2 large-aspect
benchmark used R₀ = 10, but T3 needs *many* punctures across *many* seeds and
amplitudes — at R₀ = 10 a single width-scaling probe ran for minutes. Dropping to
**R₀ = 5** (still comfortably large-aspect) halved every trace and made the rung
tractable. The T3 tests use few punctures (≈70) and coarse ds (0.07) on top of
that, since island *width* converges long before puncture count does.

---

## Validation (in `tests/test_tokamak.py`)

| Check | Assertion | Result |
|-------|-----------|--------|
| Resonance is at q = m/n | measured q at r = 0.5 is 2.0 ± 0.05 | q = 2.010 |
| Perturbation stays divergence-free | ∇·B/|B| < 10⁻³ with the mode on | ≈ 4×10⁻⁶ |
| Island opens, centred on the resonant surface | W > 0.1 and ⟨r⟩ = 0.5 ± 0.06; off-resonance surface narrower | passes |
| Island width scales as √δ | W(2e-3)/W(5e-4) ∈ (1.8, 2.2) (= √4) | ≈ 2.03 |

Stochasticity (islands → chaos via overlap) is **shown** in the visualization
rather than asserted by a scalar test — it's a qualitative, parameter-sensitive
transition; the quantitative falsifiable claims are the four above.

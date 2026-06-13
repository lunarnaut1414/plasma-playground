# T4 — the linear resistive tearing mode (reduced MHD)

Rung **T4** of [`3D_TOKAMAK_GUIDE.md`](3D_TOKAMAK_GUIDE.md), the stretch rung: the
first **self-consistent instability** of the ladder. T0–T3 prescribed fields; T4
takes an equilibrium *current sheet* and lets resistivity tear it into a magnetic
island on its own, growing at a rate set by the physics.

Kernel: `plasmaplay/tearing.py`. Viz: `tokamak_t4_viz.py`
(→ `outputs/tokamak_t4.png`). Tests: `tests/test_tearing.py`.

## Scope — why the slab, not the torus

The guide frames T4 as a reduced-MHD solver in (r, θ, φ) — a real project. The
honest laptop-scale version that is **rigorously validatable** is the **slab
(Harris-sheet) tearing mode**: equilibrium `B_y(x) = tanh(x/a)` with a strong
constant guide field (the reduced-MHD ordering). The slab keeps every piece of
tearing physics — Δ', the resistive layer, the `S^{−3/5}` growth law, the
reconnected island — while having an **exact analytic Δ'** the toroidal problem
lacks. Matching a known number beats a plausible-looking toroidal eigenmode.
(JOREK / NIMROD / M3D-C1 are the real toroidal codes; this is the teaching core.)

## Two pieces, two validations

### 1. The ideal outer region & Δ' — `delta_prime_slab`
Away from the thin resistive layer the marginal ideal equation is the Newcomb
equation `ψ'' = [k² + B_y''/B_y] ψ`, and for the tanh sheet `B_y''/B_y =
−(2/a²) sech²(x/a)` — a **Pöschl–Teller** potential with the exact decaying
solution `ψ(x) = e^{−k|x|}(1 + tanh(|x|/a)/(ka))`. Its jump in logarithmic
derivative across the sheet is the tearing stability index

```
Δ' = [ψ'(0⁺) − ψ'(0⁻)]/ψ(0) = (2/a)(1/(ka) − ka),
```

**unstable (Δ' > 0) only for k a < 1.** The numeric integrator reproduces this to
**4 digits** at every k tested (including the marginal Δ'(ka=1)=0 and the stable
ka>1 branch). This is the strongest check in the whole ladder.

### 2. The resistive layer & the `S^{−3/5}` law — `tearing_growth_rate`
The full linearised reduced-MHD eigenvalue problem (fields ψ₁, φ₁ for a mode
∝ exp(i k y + γ t)), with φ = i·φ̂ to make it real:

```
γ ψ   = −k B_y φ̂ + (1/S)(ψ'' − k² ψ)
γ Δ*φ̂ =  k B_y Δ*ψ − k B_y'' ψ,     Δ* ≡ d²/dx² − k²
```

solved as a sparse generalised eigenproblem `γ M v = A v`. Furth–Killeen–
Rosenbluth predict `γ τ_A ∝ Δ'^{4/5} S^{−3/5}` in the constant-ψ regime; the
computed slope of `ln γ` vs `ln S` is **−0.605** (ka=0.8), bang on **−3/5**.

## Two numerical traps found (and how they were fixed)

These cost real debugging time and are the heart of why the solver works:

1. **The tearing eigenvalue is *interior* to the spectrum, not extremal.** A first
   attempt with `eigs(..., which='LR')` (largest real part) returned `γ ≈ 660`
   scaling like `1/S` — a **spurious grid-scale mode**, not the tearing mode. The
   physical γ is small (10⁻³–10⁻⁴) and buried among faster damped modes. Fix:
   **shift-invert** (`eigs(sigma=…, which='LM')`) targeting the FKR estimate
   `γ ≈ 0.55 Δ'^{4/5} S^{−3/5}`, then take the most-unstable nearly-real eigenvalue
   near that shift.

2. **Under-resolving the resistive layer flattens the slope.** The layer width is
   `δ ∼ a S^{−2/5}`, so it *shrinks* with S; at S = 10⁵ it is ~0.01a. A coarse
   grid (`dx ≈ 0.05`) can't see it and the fitted slope drifted to ~−0.4. Fix: a
   **fine sparse grid** (N ≈ 3000), affordable only because shift-invert on sparse
   banded matrices is fast. Also choose **small Δ'** (ka ≈ 0.8): large Δ' (ka=0.5,
   Δ'=3) sits in the `S^{−1/3}` large-Δ' regime, *not* the constant-ψ `S^{−3/5}`
   one, so the exponent reads short.

3. (Cosmetic) the shift-invert eigenvector carries grid-scale ripple; the viz
   smooths ψ₁ with a short moving average before drawing the island. The growth
   rate and Δ' are untouched by this.

## Validation (in `tests/test_tearing.py`)

| Check | Assertion | Result |
|-------|-----------|--------|
| Δ' matches analytic tanh sheet | `|Δ'_num − Δ'_ana|/|Δ'_ana| < 1%` at ka = 0.3…0.9 | exact to 4 digits |
| Threshold at ka = 1 | Δ'(0.6) > 0, Δ'(1.4) < 0, |Δ'(1.0)| < 0.02 | passes |
| Unstable mode grows, slowly | γ(ka=0.8, S=3e4) ∈ (0, 0.1) | passes |
| Growth rate scales as S^(−3/5) | slope of ln γ–ln S ∈ (−0.68, −0.52) | −0.605 |

## Where this stops
This is the *linear* tearing mode. The nonlinear saturation (Rutherford regime,
island growth to finite width), toroidal coupling, and multi-mode disruption
dynamics are the domain of the real extended-MHD codes. T4 is the validated
linear core they all rest on.

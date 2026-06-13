# T2 — particle orbits: banana orbits, trapping & the μ invariant

Rung **T2** of [`3D_TOKAMAK_GUIDE.md`](3D_TOKAMAK_GUIDE.md). No new kernels: this
rung is launching particles into the T0 field with the existing pushers and
reading off the orbit physics.

- `plasmaplay.pushers.boris_push` (experiment 01) — gyro-resolved orbit.
- `plasmaplay.guiding_center.gc_push` (experiment 02) — the drifting guiding
  centre, including the parallel mirror force `m dv∥/dt = −μ ∇∥|B|`.

Visualised by `tokamak_t2_viz.py` (→ `outputs/tokamak_t2.png`); validated by the
`# T2` tests in `tests/test_tokamak.py`.

## The physics

On a flux surface |B| = F/R is **larger on the inboard side** (small R). Each
surface is therefore a magnetic mirror: a particle launched at the outboard
midplane with pitch λ = v∥/v small enough cannot climb the |B| hill toward the
inboard side. Conserving energy and μ, the turning (bounce) point is where

```
v∥² = v²(1 − (1−λ²) B/B_min) = 0   ⇒   trapped if  |λ| < λ_c = √(1 − B_min/B_max).
```

For the circular |B| ∝ 1/R surface, `B_min/B_max = (R₀−r)/(R₀+r)`, so

```
λ_c = √(2ε/(1+ε)),   ε = r/R₀,   ≈ √(2ε).
```

A trapped particle bounces back and forth, and its guiding centre — pushed by the
grad-B and curvature drifts (vertical on the outboard midplane) — drifts radially
off the surface, tracing a crescent **banana** in the poloidal plane. A passing
particle keeps its sign of v∥ and circulates the whole poloidal angle, staying
(to a drift width) on its surface.

For an **isotropic** velocity distribution, λ = cosθ is uniform on [−1, 1], so the
**trapped fraction equals λ_c itself**: f_t = √(2ε/(1+ε)) ≈ √(2ε). Measuring the
boundary pitch *is* measuring the trapped fraction.

## Design choices

### gc_push for the banana shape, Boris for μ
A full *gyro-resolved* banana would need the Boris step to resolve the gyro-period
(`dt ≈ T_c/30`) for the *entire bounce*. With a 1 keV proton in B ≈ 1 T,
`T_c ≈ 70 ns` while the bounce period is `~1 ms` — about **4×10⁵ Boris steps per
banana**. That's too slow for an interactive figure or a test. So:
- the **banana shape** (panel A) and the bounce / trapping classification come
  from `gc_push`, which steps on the *bounce* timescale (thousands of steps), and
- `boris_push` is used only for the **μ-invariance** panel, over ~500 gyro-orbits
  — enough to demonstrate adiabatic invariance without tracing a whole banana.

### A bracket test in the suite, full bisection only in the viz
Locating λ_c by bisection costs ~9 guiding-centre orbits *per radius* (~25 s for
the 5-radius viz panel). Too slow to do for every test run. The suite instead
**brackets** the prediction: a pitch 15 % below λ_c must trap, 15 % above must
pass, checked at r = 0.3, 0.5, 0.7. That validates `λ_c = √(2ε/(1+ε))` at three
radii (and that it grows with ε) in ~13 s, while the precise bisection lives in
the figure where wall-clock doesn't matter.

## A measurement pitfall worth recording

The obvious "banana width" proxy `R_max − R_min` over an orbit is **wrong for
classification**: a *passing* particle circulates the whole flux surface, so its
`R_max − R_min ≈ 2r` (the full surface diameter, ~1 m here) — far *larger* than a
trapped banana's radial excursion. `R_max − R_min` mixes poloidal *position* with
radial *drift*. The robust trapped/passing signal is instead:
- **v∥ changes sign** (the particle bounced) → trapped; and
- the **poloidal angle stays bounded** (`ptp(θ) < π`) for a banana vs **advances
  past 2π** for a circulating passing particle.

Both are what the tests assert.

## Validation (in `tests/test_tokamak.py`)

| Check | Assertion | Result |
|-------|-----------|--------|
| Trapping boundary = mirror prediction | 0.85 λ_c traps, 1.15 λ_c passes, at ε = 0.03, 0.05, 0.07 | passes |
| Trapped bounces, passing circulates | trapped: v∥ reverses & `ptp(θ)<π`; passing: v∥ steady & Δθ > 2π | passes |
| μ is an adiabatic invariant | `(μ_max−μ_min)/μ̄ < 1%` over ~500 gyro-orbits (Boris) | spread ≈ 0.09% |

Measured boundary (viz, bisection) vs analytic, ε = 0.03 → 0.09:
`λ_c = [0.255, 0.304, 0.345, 0.380, 0.411]` vs `√(2ε/(1+ε)) = [0.241, 0.293,
0.336, 0.374, 0.406]` — within ~6 % at the smallest ε, < 2 % at larger ε. The
residual is the expected finite-banana-width correction to the zero-orbit-width
estimate, which is relatively largest where ε (and the mirror well) is smallest.

# B2 — Nonlinear 2-D reduced MHD: a tearing mode into a magnetic island

The second rung of the MHD-instability track (NIGHT.md Track B). B1 gave the
*linear* stability of a current profile; B2 builds the **nonlinear** solver that
follows a tearing mode past the exponential phase, where it reconnects field and
grows a finite **magnetic island**. This is the laptop-scale cousin of JOREK/NIMROD.

**Status: linear phase VALIDATED (committed); nonlinear saturation is the follow-on
(B2b).** The solver runs nonlinearly, but only the linear-growth physics is asserted
by tests this rung — see the scope note.

## The model (`plasmaplay/reduced_mhd.py`)

The **Strauss reduced-MHD equations** in a 2-D slab (the periodic-y "straightened"
tokamak around one rational surface), evolving the poloidal flux ψ and the vorticity
U = ∇²φ (φ the velocity stream function, **v** = ẑ×∇φ):

    ∂ψ/∂t = −[φ, ψ] + η(∇²ψ − J_eq)
    ∂U/∂t = −[φ, U] + [ψ, ∇²ψ] + ν∇²U

with the Poisson bracket `[a,b] = a_x b_y − a_y b_x` (advection / Lorentz drive) and
`J = ∇²ψ`. The equilibrium is the **Harris sheet** `B_y0 = tanh(x)`, `J_eq = sech²x`
— the same sheet whose linear Δ′ and growth rate are analytically anchored in
`tearing.py` (T4).

Numerics: finite differences in x (Dirichlet walls), spectral (FFT) in periodic y,
an **FFT + tridiagonal elliptic solve** for φ from U (one Helmholtz solve per ky),
and an SSP-RK2 step. Normalisation: lengths in the sheet width a=1, time in τ_A,
η = 1/S. Diagnostics: the m=1 reconnected flux at the neutral line and the island
width `W = 4√(ψ_rec / B_y0'(0))`.

## Validation (falsifiable — `tests/test_reduced_mhd.py`)

- **Elliptic inversion** `∇²φ = U` recovers φ to round-off (~1e-15).
- **Instability threshold:** a seeded mode **grows for ka < 1** (Δ′ > 0) and **decays
  for ka > 1** (Δ′ < 0) — the tearing criterion, reproduced dynamically.
- **FKR scaling γ ∝ S^(−3/5):** measured from direct simulation, the growth-rate
  exponent is ≈ −0.59 (FKR −0.6) across a factor-4 range in Lundquist number — the
  same resistive-layer law the slab eigenvalue gives in T4.

## Scope boundary (stated honestly)

- **Absolute growth coefficient.** The simulated growth rate tracks the T4 *eigenvalue*
  in scaling but is ~0.54× its absolute value — an O(1) convention/discretization
  difference (initial-value vs eigenvalue, finite Δ′-layer resolution). The tests
  assert the **scaling and threshold**, which are convention-independent, *not* the
  absolute rate. The eigenvalue's own absolute value is itself FKR-seeded.
- **Nonlinear saturation (B2b).** The solver integrates the full nonlinear brackets,
  but the **Rutherford saturation** (island width W(t) following dW/dt ∝ Δ′(W) and
  saturating at Δ′(W_sat)=0) is not yet asserted by a falsifiable test — that, and the
  `tearing_island_saturation.gif` deliverable, are the follow-on rung.
- **Slab, single helicity.** A periodic-y slab around one rational surface (the
  standard tearing testbed), not the full annular (r,θ) with toroidal harmonic
  coupling — that is a later refinement.

## References
- Strauss, *Phys. Fluids* 19, 134 (1976) — the reduced-MHD equations.
- Furth, Killeen & Rosenbluth, *Phys. Fluids* 6, 459 (1963) — linear tearing & Δ′.
- Rutherford, *Phys. Fluids* 16, 1903 (1973) — nonlinear island saturation (B2b).
- Biskamp, *Nonlinear Magnetohydrodynamics* — the numerical island problem.

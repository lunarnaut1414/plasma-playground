# A3 / F3 — Transport on the real (D-shaped) equilibrium

The rung that moves the burning-plasma transport off the circular toy column and
onto the nested, shaped, outboard-shifted flux surfaces of a real Grad–Shafranov
equilibrium (NIGHT.md Track A, rung A3). It is the bridge that makes exp 09 a
"real machine" rather than a cylinder.

## The idea: keep transport 1-D, carry the geometry in two metrics

A real plasma sits on the flux surfaces ψ(R,Z) = const of an equilibrium. You do
**not** need a 2-D transport solve to run on it — transport is fast along a field
line and slow across surfaces, so it stays 1-D in a flux label ρ. The shaped
geometry enters through two flux-surface-averaged metric coefficients:

- **V'(ρ) = dV/dρ** — the volume between neighbouring flux surfaces;
- **⟨|∇ρ|²⟩** — the flux-surface average of |∇ρ|².

The transport divergence becomes `(1/V') ∂_ρ( V' ⟨|∇ρ|²⟩ n χ ∂_ρ T )`. When the
surfaces are circular this is *exactly* the cylindrical operator (V' ∝ ρ,
⟨|∇ρ|²⟩ = 1/a²) — the consistency check that pins the implementation.

## What was built

- `plasmaplay/equilibrium_metrics.py`
  - `flux_surface_metrics(R, Z, psi)` — extracts V'(ρ), ⟨|∇ρ|²⟩, V(ρ) from a
    gridded ψ using the **volume-derivative identity** for a flux-surface average,
    `⟨X⟩ = (dI_X/dρ)/(dV/dρ)`, evaluated by binning grid cells into flux-label
    shells (ρ = √ψ_n, dV = 2πR dR dZ). No fragile contour tracing — robust on a
    coarse grid. Also returns the 2-D ψ_n / ρ grids for rendering.
  - `confinement_time_ipb98(...)` — the ITER IPB98(y,2) H-mode confinement scaling.
- `plasmaplay/transport.py :: FluxSurfaceTransport1D` — subclasses `Transport1D`,
  overriding only the diffusion operator (V'/⟨|∇ρ|²⟩ finite-volume) and the volume
  average. Interpolates the metric arrays onto its uniform ρ grid. Adds
  `plasma_volume()` and `energy_confinement_time(P_loss)` diagnostics.
- `plasmaplay/animate.py :: animate_poloidal_field` — animates a 2-D field on the
  *real* (R,Z) grid (vs `animate_cross_section`'s circular revolve), so a shaped
  plasma renders with its true geometry. T(ρ,t) is mapped onto the grid via ρ(R,Z).

## Validation (falsifiable — `tests/test_equilibrium_metrics.py`, 7 tests)

- **Circular-limit metrics are analytic:** ⟨|∇ρ|²⟩ = 1/a² (mid-radius, <3%),
  V' ∝ ρ (V'/ρ constant to <5%), total V = 2π²R₀a² (<3%).
- **A real Grad–Shafranov solve** (Solov'ev source, reusing the exp-04 solver) puts
  the magnetic axis **outboard** of R₀ — the Shafranov shift (+0.28 m) — and yields a
  vertically elongated plasma (κ ≈ 1.48).
- **The flux-surface solver reduces to the cylindrical `Transport1D`** to <0.2% when
  fed circular metrics — the operator is correct.
- **A burn on the real geometry** builds a monotonically peaked profile with >10 MW
  of alpha self-heating.
- **IPB98(y,2) reproduces the ITER baseline** τ_E = 3.67 s (published ~3.7 s) and has
  the right Ip / P monotonicities.

## Deliverable

`outputs/burn_dshaped_cross_section.gif` — the headline "watch it burn" movie, now on
the **real D-shaped flux surfaces** (heating ramp → peaked ~28 keV burning core,
75 MW alpha). Regen: `python gif_gallery.py burn_dshaped_cross_section`. Also
`python experiments/09_burning_plasma/run.py --mode dshaped --save`.

## Scope boundary (stated honestly)

- **Fixed equilibrium.** The geometry is solved once; the **self-consistent Picard
  re-solve** of the equilibrium as the pressure profile evolves is *not* done here —
  it is the next sub-rung (A3b). The metrics are therefore those of the vacuum-ish
  Solov'ev source, not of the self-consistent burning pressure.
- **No β-limit (yet).** As in A2, the 1-D model has no β-limit, so the burn is
  bistable (ignite → run away hot, or stay cold). The showcase uses a sustained
  auxiliary-heated, sub-ignition operating point (χ = 0.6) to sit at a realistic
  ~28 keV. The β-limit / Greenwald / L-H operating modes are rung F3.5 (A4).
- **τ_E vs IPB98.** The IPB98 helper is validated against ITER independently. The
  toy device's simulated τ_E and its IPB98 value agree only within a factor of a
  few — χ is set for the showcase, not fit to the scaling. No overclaim of a match.
- **Shaping.** The elongation (κ ≈ 1.5) comes naturally from the Solov'ev solve on
  a tall box; strong triangularity would need a Cerfon–Freidberg analytic form.

## References
- ITER Physics Basis, *Nucl. Fusion* 39, 2175 (1999) — IPB98(y,2) scaling.
- Hirshman & Jardin, *Phys. Fluids* 22, 731 (1979) — flux-coordinate transport.
- Solov'ev / Cerfon & Freidberg, *Phys. Plasmas* 17, 032502 (2010) — analytic equilibria.

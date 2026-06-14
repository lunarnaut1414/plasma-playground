# NIGHT_LOG.md — append-only progress log for the autonomous overnight run

> First thing a resuming context reads. Newest entries at the bottom. One block per
> milestone. Format:
>
> ```
> ## <rung id> — DONE | PARTIAL | FAILED — <commit hash or "uncommitted">
> - built: <what>
> - validation: <the number achieved, or why it failed>
> - gif: <path in outputs/>
> - gotcha: <the single most useful thing learned>
> ```
>
> Resume protocol: read the last entry → `git log --oneline -15` → confirm
> `pytest -q` is green → start the next unstarted rung in NIGHT.md §5. Never redo a
> DONE rung.

---

## Starting state (handoff at charter creation, 2026-06-14)
- Branch `master`, working tree clean before the night begins.
- DONE before the night: 3-D tokamak field ladder T0–T4 (`tokamak.py`, `tearing.py`);
  burning-plasma transport exp 09 F0+F2 (`transport.py`, `tests/test_transport.py`).
- Full suite green at handoff (138 + 2 exp-09 smoke tests). `ruff check .` clean.
- First rung to start: **G1** (visualization/gif foundation) in NIGHT.md §5.

<!-- append milestone entries below this line -->

## G1 — DONE — 6e38e74
- built: `plasmaplay/animate.py` (gif foundation): `make_frames`, `torus_surface`,
  `animate_profiles`, `animate_cross_section`, `animate_torus_3d` (PillowWriter, pure
  Python). `gif_gallery.py` registry + `docs/G1_ANIMATION.md` memo. Wired into
  `plasmaplay/__init__.py`.
- validation: 1-D Gaussian-diffusion reference — mass drift 9.08e-03 (conserved),
  peak-decay-law err 2.22e-16, torus invariant <1e-9. 7 new tests; full suite
  **147 passed**, ruff clean.
- gif: `outputs/_smoke_diffusion.gif` (460K, 90 frames). Regen: `python gif_gallery.py smoke_diffusion`.
- gotcha: a prior session had written these files but never committed and left an
  unused `pytest` import that failed ruff — always run `ruff check .` before commit.
  Use `MPLBACKEND=Agg` for headless gif gen.
- next: **A1 (F1)** — He-ash + dilution + β-limit in exp 09 (`burn_0d_ignition.gif`).

## A1 (F1) — DONE — 1f3053b
- built: `transport.burn_0d_ash` — three coupled ODEs (fuel-ion n_DT, He-ash n_He,
  energy W) with fusion-born ash, fuel dilution (n_e = n_DT + 2 n_He), Z_eff-raised
  bremsstrahlung, and a soft β-limit capping the operating point. Wired into exp 09
  `run.py` (`--mode ash`) and `gif_gallery.py` (`burn_0d_ignition`).
- validation: steady T = 14.3 keV (lands in the 10–25 keV burning band ✓), ash
  balance n_He = τ_He*·R_fus = 0.991 ✓, β pinned at its 3.96% limit, f_He = 5.2%.
  Full suite **154 passed**, ruff clean.
- gif: `outputs/burn_0d_ignition.gif` (357K, regen: `python gif_gallery.py burn_0d_ignition`);
  PNG stills `burn_0d_ash.png`, `burn_0d_ignition.png`.
- gotcha: a prior session left the A1 code/tests green but UNCOMMITTED with no gif —
  the priming pass generated the gif and committed it. `gif_gallery.py` writes to
  repo-root `outputs/` (gitignored), not the experiment's `outputs/`.
- next: **A2 (F2)** — two-temperature (Te, Ti) + heating mix in exp 09
  (`burn_1d_two_temperature.gif`).

## A2 (F2.5) — DONE — c0bbd3c
- built: `TwoTempTransport1D` (subclasses `Transport1D`) evolving Te, Ti, n with
  separate χ_e/χ_i and a Braginskii electron-ion exchange term. New kernels in
  `transport.py`: `coulomb_logarithm`, `collision_frequency_ei` (NRL Spitzer ν_e),
  `equipartition_power` (Q_Δ), `equipartition_time`, `two_temperature_relax_0d`.
  Fusion/α-power use Ti, brem uses Te, NBI→ions, α→electrons (f_αe=0.85).
- validation: τ_eq(n=1e20, Te=10 keV) = 231 ms matches NRL/Spitzer (~230 ms);
  τ_eq ∝ Te^1.5/n_e; 0-D split (12,4)keV relaxes to energy-conserving mean (8 keV)
  at exactly 1/τ_eq; 1-D beam-heated discharge steadies at Ti0=24.2 > Te0=12.8 keV
  (Ti/Te=1.89). **161 passed** (7 new), ruff clean.
- gif: `outputs/burn_1d_two_temperature.gif` (657K, Te(ρ,t)+Ti(ρ,t) panels; regen
  `python gif_gallery.py burn_1d_two_temperature`); PNG `burn_1d_two_temperature.png`.
  Memo `docs/A2_TWO_TEMPERATURE.md`; `run.py --mode twotemp` prints all three checks.
- gotcha: the 1-D model has NO β-limit (that's A4), so with low χ it runs away to
  ~130 keV. The showcase is deliberately a beam-heated SUB-ignition case (χ_e=0.8,
  χ_i=0.4, n=8e19) sitting at a realistic ~24 keV. Also: τ_eq uses the *electron*
  temperature in ν_e; test the difference-decay rate from the FIRST step (instantaneous)
  and predict at Te0, else the window-average drifts ~12% and the test fails spuriously.
- next: **A3 (F3)** — transport on the real Grad–Shafranov equilibrium (exp 04
  `solvers.grad_shafranov_solve`); deliverable `burn_dshaped_cross_section.gif`.

## A3 (F3) — DONE (fixed-equilibrium); Picard re-solve deferred to A3b — 7a553fa
- built: `plasmaplay/equilibrium_metrics.py` — `flux_surface_metrics` extracts
  V'(ρ), ⟨|∇ρ|²⟩ from gridded ψ(R,Z) by the volume-derivative identity (cell
  binning, NO contour tracing); `confinement_time_ipb98` (IPB98(y,2)).
  `transport.FluxSurfaceTransport1D` (subclasses Transport1D) runs transport on
  those metrics: (1/V')∂_ρ(V'⟨|∇ρ|²⟩ nχ ∂_ρT). `animate.animate_poloidal_field`
  renders a field on the real (R,Z) cross-section. Wired: `run.py --mode dshaped`,
  gallery `burn_dshaped_cross_section`. Reuses exp-04 `grad_shafranov_solve`.
- validation: circular-limit metrics analytic (⟨|∇ρ|²⟩=1.000=1/a², V'∝ρ, total
  V=59.18 vs 59.22 torus); Solov'ev solve → Shafranov shift +0.28 m, κ=1.48; the
  flux-surface solver REDUCES to cylindrical Transport1D to <0.2%; IPB98 → ITER
  τ_E=3.67 s (pub ~3.7). **168 passed** (7 new in `tests/test_equilibrium_metrics.py`),
  ruff clean.
- gif: `outputs/burn_dshaped_cross_section.gif` (1.0M, core ~28 keV / 75 MW alpha on
  real D-shaped flux surfaces); PNG `burn_dshaped_cross_section.png`. Memo
  `docs/A3_REAL_EQUILIBRIUM.md`.
- gotcha: (1) no β-limit in 1-D → burn is BISTABLE (ignite→runaway ~65 keV, or die);
  the showcase uses a sustained aux-heated SUB-ignition point (χ=0.6, n=5e19) at a
  realistic ~28 keV. (2) numpy `ndarray.ptp()` is gone → use `np.ptp()`. (3) sim τ_E
  (0.63 s) vs IPB98 (~0.2 s) differ ~3× — χ is showcase-tuned, NOT fit to IPB98;
  reported honestly (IPB98 validated vs ITER only). (4) the volume-derivative
  flux-average identity (bin |∇ρ|²·dV into ρ-shells / bin dV) is robust on coarse
  grids — far better than tracing contours.
- next: **A3b** (self-consistent Picard equilibrium re-solve as pressure evolves)
  OR jump to **A4 (F3.5)** — Greenwald density limit, L→H transition, radiative
  collapse (`operating_modes.gif`). A4 also adds the missing 1-D β-limit. Recommend
  A4 next (higher showcase value; A3b is a refinement).

## A4 (F3.5) — DONE — 0639869
- built: `plasmaplay/operating_limits.py` — `greenwald_density` (n_G=Ip/(πa²)),
  `lh_power_threshold` (Martin 2008), `confinement_factor_lh` (smooth L→H ×2
  bifurcation), `confinement_factor_greenwald` (density-limit collapse to a floor).
  Added a `tau_factor(t,n_e,T,p_heat_density)` state hook to `burn_0d_ash`, and a
  **soft β-limit to the 1-D `Transport1D`** (B/beta_limit/beta_stiffness → raises χ
  above Troyon β; caps the A2/A3 runaway). `animate.animate_operating_space` draws
  multiple tracks on one (n,T) plane. Wired: `run.py --mode modes`, gallery
  `operating_modes`.
- validation: n_G(15MA,2m)=1.19e20 + scalings; P_LH≈52 MW at ITER point (pub ~50);
  H-mode >2× hotter than L-mode across threshold; **over-fuel collapse past n_G is
  REVERSIBLE** (22.7→0.7→21.6 keV); 1-D β-limit pins ⟨β⟩≈4% (vs 11% runaway).
  **177 passed** (9 new: 8 in `test_operating_limits.py` + 1-D β-limit in
  `test_transport.py`), ruff clean.
- gif: `outputs/operating_modes.gif` (310K, L/H/disruption sweeping n-T with the
  Greenwald line + burning band); PNG `operating_modes.png`. Memo
  `docs/A4_OPERATING_MODES.md`.
- gotcha: the 0-D burn is bistable, so the L/H/disruption scenarios need careful
  tuning — L-mode uses weak CONTINUOUS heating below P_LH (kick-then-off reverts to
  L when P drops); H-mode needs sustained heating ABOVE P_LH (alphas alone may not
  hold it after a kick). Greenwald collapse must be gentle (fuel→4.5e19, lands
  n~1.1 n_G) or n runs to ~1e22. Device S (plasma area) ≈ 4π²R₀a·√((1+κ²)/2).
- next: Track A is essentially complete (F0–F3.5). Options: **A3b** (self-consistent
  Picard equilibrium re-solve — refinement), or move to **Track B (MHD instabilities)**
  — B1 cylindrical linear MHD (internal kink / tearing on a real q(r), reuse
  `tearing.py`), deliverable `kink_eigenmode.gif`. Recommend **Track B1** (new physics,
  higher value than A3b; opens the MHD half toward the Track-C coupling).

## B1 (MHD track) — DONE — 4f1e6e5
- built: **new experiment 10** (`experiments/10_tokamak_stability/`) + new kernel
  `plasmaplay/cylinder_mhd.py` — the straight-tokamak (periodic cylinder) linear
  stability. `screw_pinch_q` (q(r) from (1-r²)^ν current; q(0)=q0, q(a)=(ν+1)q0),
  `rational_surface` (q=m/n by bisection), `delta_prime_cylinder` (outer Newcomb
  eqn → tearing index Δ′), `internal_kink_unstable`+`internal_kink_xi` (m=1 q(0)<1
  sawtooth trigger + top-hat eigenfunction), `fkr_growth_rate` (FKR S^-3/5, reused
  from slab/T4). Registered `tearing` + `cylinder_mhd` in `__init__.py`.
- validation: q(0)=q0, q(a)=(ν+1)q0; **m=1 kink unstable iff q(0)<1** (0.7/0.85/0.95
  unstable, 1.05/1.3 stable); **sign of Δ′ predicts tearing stability** + Δ′ falls
  with m (gap-robust); **γ∝S^-3/5** (10^-0.6/decade). **191 passed** (14 new in
  `tests/test_cylinder_mhd.py`), ruff clean.
- gif: `outputs/kink_eigenmode.gif` (715K, ξ_r(r)+q(r) panel + the m=1 core-shift
  crescent growing); PNG `kink_eigenmode.png`. Memo `docs/B1_CYLINDER_MHD.md`.
- gotcha: the **absolute Δ′ is resolution-dependent** near the singular layer (the
  outer ψ has a log term that only cancels as gap→0; values grow ~5.8→8.2 as gap
  4e-3→1e-3) — so tests assert only the **sign and m-ordering** (the charter's gate),
  NOT an absolute Δ′. R/B_θ scale cancels in the Newcomb drive term → stability
  depends only on the q-profile. solve_ivp t_eval must be sorted in the integration
  direction (decreasing when integrating inward from the wall).
- next: **B2** — nonlinear 2-D reduced MHD (`plasmaplay/reduced_mhd.py`): evolve ψ &
  vorticity in (r,θ), watch a tearing island grow and SATURATE (Rutherford dW/dt∝Δ′(W));
  deliverable `tearing_island_saturation.gif`. Then B3 sawtooth cycle, then Track C
  (couple a sawtooth/tearing event into the exp-09 burn). NOTE: B2 is a bigger rung
  (2-D nonlinear PDE) — may need a validated partial across two wakes.

## B2a (MHD track) — PARTIAL (linear phase DONE; saturation = B2b) — 3e10d4b
- built: `plasmaplay/reduced_mhd.py` — `ReducedMHD`: the Strauss reduced-MHD eqns
  (ψ, vorticity U=∇²φ) on a 2-D slab, x finite-difference + y spectral (FFT), a
  **vectorized** FFT+tridiagonal elliptic solve for φ from U (`_thomas_vec` over all
  ky at once — ~2× faster), SSP-RK2. Harris sheet B_y0=tanh(x) (same as T4
  `tearing.py`). Diagnostics: m=1 reconnected flux + island width W=4√(ψ_rec).
- validation (3 tests, `tests/test_reduced_mhd.py`): elliptic inversion exact (~1e-15);
  seeded mode grows for ka<1, decays for ka>1; **γ∝S^-3/5** measured exponent −0.583
  (FKR −0.6) by direct simulation, across a factor-4 in S. **194 passed**, ruff clean.
- gif: `outputs/tearing_island.gif` (920K, the sheet tearing into an island; named
  honestly — NOT "saturation"); PNG `tearing_island.png` (`run.py --island`). Memo
  `docs/B2_REDUCED_MHD.md`.
- gotcha: (1) the simulated **absolute growth rate is ~0.54× the T4 eigenvalue** — an
  O(1) convention/discretization difference (IVP vs eigenvalue; the eigenvalue is
  itself FKR-seeded). So tests assert the **S^-3/5 scaling + threshold**, NOT the
  absolute γ. Never forced the match. (2) tearing.py eigenvalue is **inviscid** → set
  Pm=0 to compare. (3) resistive layer δ~S^-2/5 must be resolved (δ/dx≳3): nx≈224 over
  Lx=4 at S=400-1600. (4) the per-ky python tridiagonal loop was the bottleneck —
  vectorizing across ky cut test time ~2×.
- next: **B2b** — nonlinear **Rutherford saturation**: run the solver long, show island
  width W(t) saturates (dW/dt→0) and follows dW/dt∝Δ′(W); validate + make
  `tearing_island_saturation.gif`. (A ~700-step-to-t700 run was started this wake but
  not completed; the solver is ready, just needs a long run + a falsifiable saturation
  check.) Then **B3** sawtooth cycle, then **Track C** (couple MHD event into exp-09 burn).

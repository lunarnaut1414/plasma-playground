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

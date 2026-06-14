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

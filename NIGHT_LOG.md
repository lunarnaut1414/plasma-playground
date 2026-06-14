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

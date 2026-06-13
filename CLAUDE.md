# CLAUDE.md — guide for Claude Code working in this repo

Read this first. It's the high-signal map; the linked docs have the detail.

## What this repo is

`plasma-playground` — a collection of **self-contained plasma physics simulation
experiments** for self-teaching (stellarator, tokamak, MHD space drive, PIC, waves).
Pedagogy over performance. Python + the open-source scientific/plasma stack. Developed
on an Apple-Silicon Mac (M2 Max, 96 GB), single machine, mostly CPU.

## Orient yourself (read in this order)

1. `README.md` — human-facing overview + experiment roadmap.
2. `docs/FIDELITY.md` — the **F0→F4 fidelity ladder** every experiment is organized around.
3. `docs/FUNDAMENTALS.md` — the **shared kernels + validation test suite to build first** (the engineering plan).
4. `experiments/NN_*/PLAN.md` — per-experiment fidelity ladder (what to build, in what order, how to validate each rung).
5. `docs/3D_TOKAMAK_GUIDE.md` — cross-experiment roadmap toward a 3-D tokamak (T0–T4), reusing existing kernels.

## Layout

```
plasmaplay/                 # the shared package — reusable, validated kernels
experiments/NN_name/        # one self-contained experiment each:
  ├── run.py                #   entry point (argparse; --save writes to ./outputs/)
  ├── README.md             #   what it shows + how to run + concepts
  └── PLAN.md               #   fidelity ladder F0→F4 with per-rung validation
docs/                       # FIDELITY.md, FUNDAMENTALS.md
outputs/                    # generated figures/animations (gitignored — never commit)
tests/                      # pytest validation suite (see FUNDAMENTALS.md §2)
```

## The `plasmaplay` API (current surface)

- `constants` — SI constants + `e, c, m_e, m_p`; `thermal_velocity(T_eV, mass)`,
  `gyrofrequency(q, B, m)`, `gyroradius(v_perp, q, B, m)`.
- `fields` — field models, each a **factory returning a callable** `f(position) -> (3,) ndarray`:
  `uniform_B(Bz)`, `uniform_E(Ex,Ey,Ez)`, `zero_field()`, `magnetic_mirror(B0, mirror_ratio, length)`.
- `pushers` — `boris_push(position, velocity, charge, mass, E_func, B_func, dt, n_steps)`
  → `(t, positions, velocities)`. Energy-conserving; the reference particle integrator.
- `plotting` — `ensure_outputs_dir(__file__)`, `plot_trajectory_3d(positions, ...)`.

Planned modules (see `docs/FUNDAMENTALS.md` for the build order): `integrators.py`,
`solvers.py` (Poisson + elliptic), `fvm.py` (finite-volume MHD), `pic.py`, `collisions.py`,
`dispersion.py`, `diagnostics.py`.

## Conventions (follow these when adding code)

- **SI units everywhere.** Use `astropy.units` / `plasmapy` only for cross-checks, not in hot loops.
- **Fields are callables** `f(pos) -> (3,) ndarray`. Any field must drop into any pusher/tracer unchanged.
- **Kernels are pure** — arrays in, arrays out, no global state. Put them in `plasmaplay/`, not in an experiment.
- **Experiments are self-contained** — a `run.py` + `README.md` + `PLAN.md` per folder. Shared logic goes up into `plasmaplay/`, never copy-pasted between experiments.
- **Validation-first** — every kernel has a falsifiable test (a number it must hit) before an experiment relies on it. No "the plot looks plasma-ish."
- **Docstrings explain the physics**, not just the code — match the voice in `pushers.py`/`fields.py` (e.g. *why* Boris over RK4).
- **Performance**: start in plain NumPy. Reach for `numba` only on proven hot loops (PIC pushes, Biot–Savart, finite-volume). JAX optional; CPU is the reliable path on macOS.

## How experiments depend on each other

Not silos — a web. Honor these couplings (reuse, don't duplicate):

- 01 ⟶ 02: the Boris pusher and fields feed the guiding-center comparison.
- 04 ⟶ 01/02: a real equilibrium B-field becomes the F3 field for banana orbits.
- 03 ⟶ 08: PIC field output is FFT'd to *measure* the dispersion relation.
- 06 ⟶ 07: the MHD finite-volume solver is the engine of the space-drive thruster.
- shared: ω–k FFT, Poisson solver, conservation monitors are used by several.

## Common commands

```bash
# Env — local .venv via uv (preferred). ACTIVATE before installing, else uv
# targets an active base conda env instead of the venv.
uv venv --python 3.11 .venv && source .venv/bin/activate && uv pip install -e ".[dev]"
# conda alternative: conda env create -f environment.yml && conda activate plasma-playground

# Always work inside the activated .venv (the dev machine has base conda active by
# default). If unsure, use the explicit interpreter: .venv/bin/python ...

# Run an experiment
python experiments/01_single_particle_motion/run.py [--save]
# Headless figure check (no GUI): prefix MPLBACKEND=Agg

# Tests / lint
pytest                      # validation suite (FUNDAMENTALS.md §2); 37 passing
ruff check .                # lint (line length 100, py311)
```

## When asked to "build experiment NN" or "add fidelity rung FX"

1. Open that experiment's `PLAN.md`; find the rung and its **Validation** line.
2. Check `docs/FUNDAMENTALS.md` — is the needed kernel built & validated? If not, build+test the kernel in `plasmaplay/` first.
3. Implement the rung in the experiment's `run.py`; keep reusable parts in `plasmaplay/`.
4. Make the validation check pass and **print/plot it** so it's visible when run.
5. Update the experiment `PLAN.md` status marker and `README.md` roadmap if a rung is newly done.

## Git

Public GitHub repo `lunarnaut1414/plasma-playground` (branch `master`). Commit/push only
when the user asks. `outputs/` and figures are gitignored by design.

# plasma-playground

A personal collection of **fun plasma physics simulation experiments** — a place
to self-teach (even if superficially) how you model and simulate magnetized
plasmas: stellarators, tokamaks, MHD drives for space propulsion, and whatever
else looks interesting.

Built on Python and the open-source scientific/plasma ecosystem. Everything here
runs comfortably on a laptop CPU (developed on an M2 Max); heavier experiments
note where a GPU or more cores would help.

> Philosophy: each experiment is small, self-contained, and explains the physics
> it demonstrates. Favor understanding over performance. Break things, change the
> knobs, see what happens.

## Quickstart

```bash
# Option A — local venv via uv (recommended: fast, isolated, reads pyproject.toml)
uv venv --python 3.11 .venv
source .venv/bin/activate            # activate FIRST so tools target the venv,
uv pip install -e ".[dev]"           #   not an active base conda env
# (add extras as needed: ".[dev,notebooks]", ".[equilibria]", ".[stellarator]")

# Option B — conda (handy for the compiled stack on Apple Silicon)
conda env create -f environment.yml
conda activate plasma-playground

# Run the first experiment + the validation suite
python experiments/01_single_particle_motion/run.py
pytest
```

> Heads-up: if a base conda env is active, `uv` will install into *it* unless you
> `source .venv/bin/activate` first (or pass `uv pip install --python .venv/bin/python ...`).
> The `.venv/` is gitignored — it stays local to your machine.

## Repository layout

```
plasma-playground/
├── plasmaplay/        # shared utilities (constants, field models, pushers, plotting)
├── experiments/       # one folder per self-contained experiment, each with a README
│   └── 01_single_particle_motion/
├── notebooks/         # exploratory Jupyter notebooks
└── outputs/           # generated figures/animations (gitignored)
```

## Roadmap of experiments

A rough self-teaching path, from a single particle up to full devices. Add,
reorder, and skip freely. **Each experiment has a `PLAN.md`** laying out a
*fidelity ladder* (F0 → F4): how to grow that experiment from an analytic sanity
check to a research-grade simulation, one assumption at a time. The shared ladder
is defined in [`docs/FIDELITY.md`](docs/FIDELITY.md).

| # | Experiment | Core idea | Plan | Likely tools |
|---|------------|-----------|------|--------------|
| 01 | Single-particle motion ✅ | gyro-orbits, drifts, mirror trapping | [plan](experiments/01_single_particle_motion/PLAN.md) | numpy (Boris pusher) |
| 02 | Guiding-center drifts ✅ | grad-B & curvature drifts, invariants | [plan](experiments/02_guiding_center_drifts/PLAN.md) | numpy / plasmapy |
| 03 | Many particles / PIC 1D ✅ | electrostatic waves, Landau damping | [plan](experiments/03_pic_1d/PLAN.md) | numpy + numba |
| 04 | Tokamak equilibrium ✅ | Grad–Shafranov, flux surfaces, q-profile | [plan](experiments/04_tokamak_equilibrium/PLAN.md) | [FreeGS](https://github.com/freegs-plasma/freegs) |
| 05 | Stellarator field lines ✅ | 3D fields, Poincaré plots, rotational transform | [plan](experiments/05_stellarator_field_lines/PLAN.md) | [simsopt](https://github.com/hiddenSymmetries/simsopt) |
| 06 | Ideal MHD basics ✅ | continuity/momentum/induction, Alfvén waves | [plan](experiments/06_ideal_mhd/PLAN.md) | numpy / scipy |
| 07 | MHD accelerator (space drive) ✅ | J×B body force, magnetic nozzle concept | [plan](experiments/07_mhd_space_drive/PLAN.md) | numpy / scipy |
| 08 | Plasma waves & dispersion ✅ | cold/warm plasma dispersion relations | [plan](experiments/08_plasma_waves_dispersion/PLAN.md) | [PlasmaPy](https://docs.plasmapy.org) |
| 09 | Burning plasma (transport) ✅ | ignition → steady burn → fuel injection; He-ash & β-limit; two-temperature Te/Ti; real D-shaped equilibrium; operating modes (L/H/disruption) | [plan](experiments/09_burning_plasma/PLAN.md) | numpy (ASTRA-like) |

✅ = implemented. Each plan marks which fidelity rungs are done.

### The 3-D tokamak ladder (T0 → T4) ✅ built

[`docs/3D_TOKAMAK_GUIDE.md`](docs/3D_TOKAMAK_GUIDE.md) lays out a laptop-reachable
path from the experiment-04 equilibrium to a **3-D tokamak** playground. The whole
ladder is now built and validated (`plasmaplay/tokamak.py`, `plasmaplay/tearing.py`,
`tokamak_t{0..4}_viz.py`, `tests/test_tokamak.py` + `tests/test_tearing.py`):

| Rung | What it shows | Validation | Memo |
|------|---------------|------------|------|
| T0 | 3-D equilibrium B-field from ψ(R,Z) | ∇·B≈0, \|B\|∝1/R, axis purely toroidal | (in `tokamak.py`) |
| T1 | field-line topology: Poincaré + q-profile | q matches large-aspect analytic; nested surfaces | [T1](docs/T1_QPROFILE_POINCARE.md) |
| T2 | banana orbits, trapping, μ invariant | trapped fraction √(2ε/(1+ε)); μ steady to 0.1% | [T2](docs/T2_BANANA_ORBITS.md) |
| T3 | magnetic islands & stochasticity (real 3-D) | island at q=m/n; W∝√δ; ∇·B kept ~1e-6 | [T3](docs/T3_MAGNETIC_ISLANDS.md) |
| T4 | linear resistive tearing mode (reduced MHD) | exact analytic Δ'; γ∝S^(−3/5) (slope −0.605) | [T4](docs/T4_TEARING_MODE.md) |

Honest scope boundaries (nonlinear/turbulent/whole-device) and where research
codes (JOREK / GENE / XGC) take over are spelled out in the guide. Run any rung's
figure with e.g. `MPLBACKEND=Agg python tokamak_t1_viz.py`.

## The open-source plasma stack (what's out there)

- **[PlasmaPy](https://docs.plasmapy.org)** — formulary, particle tracking,
  dispersion solvers; great for sanity-checking your own code.
- **[FreeGS](https://github.com/freegs-plasma/freegs)** — free-boundary
  Grad–Shafranov solver for tokamak equilibria.
- **[simsopt](https://github.com/hiddenSymmetries/simsopt)** /
  **[DESC](https://desc-docs.readthedocs.io)** — stellarator optimization and
  3D MHD equilibria (DESC is JAX-based).
- **astropy.units** — carry physical units through your math so you catch
  mistakes.
- **numba / JAX** — speed up the hot loops; both work on Apple Silicon.

## A note on compute (M2 Max, 96 GB)

Plenty for 1D/2D PIC, field-line tracing, and most equilibrium solves on CPU.
The 96 GB of unified memory is the standout — you can hold large 2D/3D grids in
RAM. JAX has an experimental Metal backend, but for now the CPU path is the most
reliable on macOS; reach for `numba` to accelerate explicit Python loops.

## License

MIT — see `pyproject.toml`. Have fun.

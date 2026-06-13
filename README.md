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
# Option A — conda (recommended on Apple Silicon for the compiled stack)
conda env create -f environment.yml
conda activate plasma-playground

# Option B — uv / pip
python -m venv .venv && source .venv/bin/activate
pip install -e ".[notebooks]"

# Run the first experiment
cd experiments/01_single_particle_motion
python run.py
```

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
reorder, and skip freely.

| # | Experiment | Core idea | Likely tools |
|---|------------|-----------|--------------|
| 01 | Single-particle motion ✅ | gyro-orbits, drifts, mirror trapping | numpy (Boris pusher) |
| 02 | Guiding-center drifts | grad-B & curvature drifts, invariants | numpy / plasmapy |
| 03 | Many particles / PIC 1D | electrostatic waves, Landau damping | numpy + numba |
| 04 | Tokamak equilibrium | Grad–Shafranov, flux surfaces, q-profile | [FreeGS](https://github.com/freegs-plasma/freegs) |
| 05 | Stellarator field lines | 3D fields, Poincaré plots, rotational transform | [simsopt](https://github.com/hiddenSymmetries/simsopt) |
| 06 | Ideal MHD basics | continuity/momentum/induction, Alfvén waves | numpy / scipy |
| 07 | MHD accelerator (space drive) | J×B body force, magnetic nozzle concept | numpy / scipy |
| 08 | Plasma waves & dispersion | cold/warm plasma dispersion relations | [PlasmaPy](https://docs.plasmapy.org) |

✅ = implemented.

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

# The fidelity ladder

Every experiment in this playground is organized as a **fidelity ladder**: a
sequence of versions of the same experiment, each one more physically faithful
(and more expensive) than the last. The point is pedagogical — you climb the
ladder one rung at a time, and at each rung you understand *exactly* which
simplifying assumption you just relaxed and what new physics appeared because of
it.

You do **not** need to finish every rung. Climbing two or three rungs of an
experiment and moving on is a perfectly good way to learn. The ladder just makes
the "what would make this more real?" question explicit and answerable.

## The five levels

| Level | Name | What it is | Typical assumption being relaxed |
|-------|------|------------|----------------------------------|
| **F0** | Analytic / sanity | Closed-form formulas or a 0-D estimate. No real simulation — you compute known results to anchor everything else. | — (this *is* the reference) |
| **F1** | Minimal numerical | The simplest possible simulation. Reduced dimensionality, prescribed (non-self-consistent) fields, linear, smallest useful method. | "Solve it numerically instead of on paper." |
| **F2** | Self-consistent / textbook | The standard method you'd find in a textbook chapter. Fields respond to the plasma; nonlinear effects appear; usually 1-D or 2-D. | "Let the plasma and fields affect each other." |
| **F3** | Realistic geometry / library | Real device geometry and/or a community library. 2-D or 3-D, real coil sets or equilibria, validated building blocks. | "Use a real machine's shape and proper tools." |
| **F4** | Research-grade | Coupled physics, performance work, and validation against real data or reference codes. Where a hobby project meets the literature. | "Make it good enough to compare with reality." |

Not every experiment has a meaningful version at all five levels — some top out at
F3, some have two distinct F2 variants. Each plan lists only the rungs that make
sense for that topic and says so.

## What each rung specifies

For every rung, the experiment's `PLAN.md` answers the same six questions:

- **Models** — what physics is actually represented.
- **Assumes** — the simplifications that define this rung (and that the next rung relaxes).
- **Method & tools** — numerical method and which Python libraries do the work.
- **You'll learn** — the concept this rung is meant to teach.
- **Validation** — the concrete check that tells you the code is correct (compare to an F0 formula, a reference solution, or a published number).
- **Compute (M2 Max)** — rough cost on the dev machine: seconds/minutes/hours, CPU vs. would-benefit-from-acceleration, memory.

## A note on validation

Every rung must be **falsifiable**. F0 gives you a number; F1+ must reproduce it
(or a published equivalent) before you trust anything new the simulation tells
you. "The plot looks plasma-ish" is not validation. This habit — always have a
check — is most of what separates a simulation you can learn from from one that
quietly lies to you.

## How this maps to compute on the M2 Max (96 GB)

- **F0–F2** are almost always **seconds on a single CPU core**. Pure NumPy/SciPy.
- **F2–F3** is where `numba` (JIT the hot loops) and the 96 GB of unified memory
  start to matter — 1-D/2-D PIC, field-line tracing, equilibrium solves.
- **F4** is where you'd want either real parallelism (`numba` parallel, `jax`)
  or an established external code (Smilei, WarpX, DESC, FreeGS). Some F4 rungs
  are explicitly "run a real research code and compare," not "write your own."

See each experiment's `PLAN.md` for the concrete ladder.

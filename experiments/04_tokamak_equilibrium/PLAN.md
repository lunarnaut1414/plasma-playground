# 04 — Tokamak equilibrium — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** not started. First *device-scale* experiment — "I simulated a tokamak."

## The question

What shape does a confined plasma actually settle into? Given a balance between
plasma pressure pushing out and magnetic force pushing in, solve for the nested
magnetic flux surfaces — the equilibrium that every tokamak operates around.

## Why it matters (real devices)

The Grad–Shafranov equation is *the* equation of tokamak design. Every shot on
every tokamak is planned around a GS equilibrium: where the plasma boundary sits,
where the X-point and divertor strike points land, what the safety factor q
profile looks like. This is the most directly "real engineering" experiment in
the playground.

## Prerequisites

Basic MHD force balance (∇p = J×B) — a paragraph of experiment 06 F0 is enough.
2-D finite differences / sparse linear solves. Concept of a flux surface.

## Fidelity ladder

### F0 — Analytic Solov'ev equilibrium
- **Models:** the Solov'ev solution — a closed-form exact solution of Grad–Shafranov for special (linear) profile choices.
- **Assumes:** specific p'(ψ) and FF'(ψ) that make GS linear and solvable on paper.
- **Method & tools:** evaluate the analytic ψ(R,Z) on a grid; matplotlib contours.
- **You'll learn:** what a flux function ψ is and that its contours are nested flux surfaces; a ground-truth equilibrium to test your solver against.
- **Validation:** this *is* the reference for F1.
- **Compute:** instant.

### F1 — Fixed-boundary Grad–Shafranov solver
- **Models:** solve the nonlinear GS PDE Δ*ψ = -μ₀R²p'(ψ) - FF'(ψ) inside a fixed boundary, with chosen pressure and current profiles.
- **Assumes:** plasma boundary is *prescribed* (no external coils yet); axisymmetric.
- **Method & tools:** finite-difference Δ* operator + sparse solve (`scipy.sparse`), Picard iteration on the nonlinearity.
- **You'll learn:** how an equilibrium is actually computed; the role of the two free profile functions; what the safety factor q(ψ) is and how to extract it.
- **Validation:** with Solov'ev profiles your solver must reproduce the F0 analytic ψ to grid accuracy.
- **Compute:** seconds; pure SciPy, single core.

### F2 — Free-boundary equilibrium with coils (FreeGS)
- **Models:** put real poloidal-field coils outside the plasma; let the boundary, X-point, and divertor *emerge* from the coil currents.
- **Assumes:** axisymmetric; you specify coil currents (or constrain shape and solve for them).
- **Method & tools:** **FreeGS** (`pip install freegs`) — the right tool; don't reinvent the coil Green's functions.
- **You'll learn:** how shaping coils make elongated, diverted plasmas; what "free-boundary" means; X-points and the separatrix.
- **Validation:** reproduce a FreeGS tutorial equilibrium; check q₀, q₉₅, elongation against expected values.
- **Compute:** seconds–minutes.

### F3 — A real machine
- **Models:** build a free-boundary equilibrium matching a real device's coil set and parameters — e.g. MAST-U, DIII-D, or ITER.
- **Assumes:** axisymmetric; published coil geometry and target plasma parameters.
- **Method & tools:** FreeGS with a real machine description; compute q-profile, β, shaping factors.
- **You'll learn:** how close a few-hundred-line script can get to a real tokamak cross-section; the parameter space designers work in.
- **Validation:** boundary shape, q₉₅, and β within a sensible margin of published values for the device.
- **Compute:** minutes.

### F4 — Time evolution or stability
- **Models:** either (a) couple the equilibrium to a simple transport step and evolve the profiles in time, or (b) feed the equilibrium into a stability analysis (ballooning/kink criterion, or a Mercier check).
- **Assumes:** depends on path; this is where equilibrium meets the next physics layer.
- **Method & tools:** your own transport loop on top of FreeGS, or an external stability code; optionally an EFIT-style reconstruction from synthetic diagnostics.
- **You'll learn:** that equilibrium is only the start — real plasmas evolve and can go unstable; what sets operational limits (β-limit, q-limit).
- **Validation:** stability boundary matches a known scaling; reconstructed equilibrium recovers the input within diagnostic error.
- **Compute:** minutes–hours.

## Diagnostics you'll reuse
Flux-surface contour plots, q(ψ) profiles, pressure/current profiles, separatrix & X-point location, β and shaping metrics.

## Key references
- Freidberg, *Ideal MHD* (Grad–Shafranov, Solov'ev).
- Wesson, *Tokamaks* — the reference handbook.
- FreeGS docs and examples: github.com/freegs-plasma/freegs.

## Stretch goals
- Reproduce the iconic "D-shaped" diverted plasma and label the strike points.
- Scan coil currents and watch the X-point move.
- Use your equilibrium's B-field as the F3 field for experiments 01/02 (banana orbits in a *real* equilibrium).

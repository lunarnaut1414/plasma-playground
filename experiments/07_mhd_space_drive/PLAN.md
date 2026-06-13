# 07 — MHD drive for space propulsion — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** not started. The "why I started this repo" experiment. Builds on 06.

## The question

Can we accelerate plasma to produce thrust using the J×B (Lorentz) body force?
Model an electromagnetic plasma accelerator — a magnetoplasmadynamic (MPD)
thruster or a magnetic-nozzle drive — and compute thrust, exhaust velocity, and
efficiency. This is electric propulsion at its most powerful.

## Why it matters (real devices)

MPD thrusters and magnetic-nozzle drives (VASIMR, applied-field MPD, helicon
thrusters) are leading candidates for high-power in-space propulsion — far higher
thrust density than ion engines, with specific impulse chemical rockets can't
touch. The physics is exactly the J×B force from experiment 06 applied to make
exhaust. This is where the playground's plasma physics meets a real engineering
goal.

## Prerequisites

Experiment 06 (ideal MHD, the J×B force). Rocket basics: thrust, exhaust velocity,
specific impulse Isp, the rocket equation. Helpful: experiment 01 (drifts) for
the magnetic-nozzle detachment question.

## Fidelity ladder

### F0 — Analytic / scoping
- **Models:** thrust from J×B; the self-field MPD scaling (Maecker formula) T ≈ (μ₀/4π) I² [ln(r_a/r_c) + ¾]; exhaust velocity, Isp = v_e/g₀, input power, and thrust efficiency.
- **Assumes:** lumped 0-D; current I and geometry prescribed.
- **Method & tools:** arithmetic; the rocket equation.
- **You'll learn:** the headline trade — thrust vs. Isp vs. power; why MPD wants megawatts; ballpark numbers for a real thruster.
- **Validation:** reproduce published MPD thrust/power numbers (e.g. ~N-class thrust at ~MW) within a factor.
- **Compute:** instant.

### F1 — 1-D channel acceleration model
- **Models:** a 1-D plasma flow through a channel with prescribed current density J and magnetic field B; integrate the J×B acceleration to get the exit velocity and thrust.
- **Assumes:** prescribed (non-self-consistent) J and B; quasi-1-D, steady, single fluid.
- **Method & tools:** NumPy ODE/quadrature along the channel.
- **You'll learn:** how the body force converts electrical input into directed kinetic energy; where the acceleration happens along the channel.
- **Validation:** integrated thrust matches the F0 J×B estimate for the same I and geometry.
- **Compute:** seconds.

### F2 — Self-consistent 1-D/quasi-2-D MHD accelerator
- **Models:** couple the flow to the electromagnetic field — solve resistive MHD in an accelerator channel so the current distribution and B are computed, not imposed; get thrust and efficiency self-consistently.
- **Assumes:** 1-D or quasi-2-D, resistive single-fluid MHD, simple equation of state.
- **Method & tools:** extend the experiment-06 MHD solver with Ohm's law + applied current boundary conditions; `numba`.
- **You'll learn:** that current distribution and back-EMF set the real efficiency; resistive dissipation as a loss channel.
- **Validation:** efficiency and thrust trends vs. current match MPD experimental scaling (e.g. thrust ∝ I²).
- **Compute:** seconds–minutes.

### F3 — 2-D axisymmetric thruster / magnetic nozzle
- **Models:** 2-D axisymmetric MHD of either an applied-field MPD thruster or a magnetic-nozzle drive; include the diverging applied field and the question of plasma *detachment* from the nozzle.
- **Assumes:** axisymmetric; single-fluid (optionally Hall term); given coil/electrode geometry.
- **Method & tools:** a 2-D axisymmetric MHD solver (your F3 from experiment 06, adapted) or an open plasma-propulsion code; `numba`/parallel.
- **You'll learn:** how a magnetic nozzle converts thermal/azimuthal energy into axial thrust; the detachment problem (how does plasma leave the field?); thrust vectoring.
- **Validation:** nozzle thrust gain and plume divergence match magnetic-nozzle theory/experiments; detachment location consistent with the literature.
- **Compute:** minutes–hours; 96 GB RAM helps the 2-D grid.

### F4 — Research-grade thruster model
- **Models:** add the physics real designs need — finite resistivity & Hall effect, a real propellant equation of state (argon/xenon/lithium), electrode/wall losses; compare against a published thruster's measured performance.
- **Assumes:** the chosen high-fidelity model; validation against experimental data.
- **Method & tools:** extended MHD (Hall-MHD) solver or an established multi-fluid/PIC-hybrid code; real propellant data.
- **You'll learn:** the gap between idealized MHD and a thruster that actually works; where the losses really are; why high-power EP is hard.
- **Validation:** computed thrust, Isp, and efficiency within experimental error of a real device (e.g. a lithium-fed applied-field MPD or VASIMR data point).
- **Compute:** hours; this is the genuinely hard rung.

## Diagnostics you'll reuse
Thrust / Isp / efficiency vs. current and power curves, velocity and current-density fields along the channel, magnetic-nozzle field-line + plume plots, energy-balance (input power → kinetic + losses) breakdown.

## Key references
- Jahn, *Physics of Electric Propulsion* — the classic EP text (MPD chapters).
- Goebel & Katz, *Fundamentals of Electric Propulsion* (JPL).
- Review articles on magnetic-nozzle detachment (e.g. Ahedo & Merino) and applied-field MPD scaling.

## Stretch goals
- Build a small design tool: input power & propellant → predicted thrust/Isp across the ladder's models.
- Compare your J×B drive's Isp/thrust envelope against ion and Hall thrusters on a T–Isp chart.
- Use experiment-01 guiding-center drifts to study single-ion behavior in your magnetic nozzle.

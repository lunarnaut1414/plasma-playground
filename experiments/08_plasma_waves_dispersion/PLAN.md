# 08 — Plasma waves & dispersion — Plan & fidelity ladder

> Fidelity ladder defined in [`docs/FIDELITY.md`](../../docs/FIDELITY.md).
> **Status:** F2 + F3 implemented (`run.py`) — kinetic Langmuir dispersion +
> Landau damping, and ω–k measured from a PIC run. Ties the playground together.

## The question

A plasma supports a zoo of waves — Langmuir, ion-acoustic, Alfvén, whistler,
Bernstein, and more — each with its own dispersion relation ω(k). Can we compute
these dispersion relations across models (cold → warm → kinetic), and then
*measure* them directly out of the simulations from experiments 03 and 06?

## Why it matters (real devices)

Waves are how plasmas are heated (ICRH, ECRH in tokamaks), diagnosed (reflectometry,
interferometry), and how they go unstable. The dispersion relation is the
fingerprint that tells you which wave you're looking at. This experiment is also
the unifier: it connects the kinetic (03) and fluid (06) pictures by showing they
predict the *same* waves in their overlapping regime.

## Prerequisites

Fourier analysis (this experiment lives in ω–k space). Helpful: experiment 03
(to measure dispersion from PIC) and experiment 06 (for the MHD waves).
The cold-plasma F0 of experiment 03 overlaps here.

## Fidelity ladder

### F0 — Analytic / sanity
- **Models:** the standard dispersion relations — Langmuir ω² = ω_pe² + 3k²v_th², ion-acoustic, EM waves with cutoff at ω_pe, and the magnetized cold-plasma relations (cyclotron resonances).
- **Assumes:** linear, often cold or simple Maxwellian.
- **Method & tools:** PlasmaPy `plasmapy.formulary` and `plasmapy.dispersion` for cross-checks.
- **You'll learn:** the named waves and their limiting behaviors (cutoffs, resonances).
- **Validation:** reference curves the later rungs reproduce.
- **Compute:** instant.

### F1 — Cold-plasma dispersion solver
- **Models:** numerically root-find ω(k) for the cold magnetized-plasma dispersion relation at arbitrary propagation angle; plot dispersion curves and the CMA / Clemmow–Mullaly–Allis diagram.
- **Assumes:** cold plasma (no thermal/kinetic effects); linear.
- **Method & tools:** NumPy/SciPy root-finding; the cold-plasma dielectric tensor.
- **You'll learn:** how cutoffs and resonances organize wave propagation; the full cold-plasma wave landscape in one diagram.
- **Validation:** curves pass through the F0 cutoffs/resonances at the right frequencies.
- **Compute:** seconds.

### F2 — Kinetic dispersion (Landau damping built in)  ✅ implemented
- **Models:** solve the *kinetic* (Vlasov) dispersion relation using the plasma dispersion function Z(ζ); get complex ω → real frequency **and** Landau damping rate.
- **Assumes:** Maxwellian, linear; unmagnetized or magnetized depending on how far you push.
- **Method & tools:** SciPy + `scipy.special.wofz` for Z(ζ); complex root-finding; PlasmaPy's kinetic dispersion solvers as a check.
- **You'll learn:** that damping falls out of the *same* dispersion relation as the frequency — the kinetic content that fluid theory misses; why warm plasma ≠ cold plasma.
- **Validation:** Im(ω) matches the experiment-03 F0 Landau rate; Re(ω) matches Bohm–Gross.
- **Compute:** seconds.

### F3 — Measure dispersion from a simulation  ✅ implemented
- **Models:** take field data E(x, t) from your experiment-03 PIC (or experiment-06 MHD) run, 2-D FFT it to the ω–k plane, and watch the dispersion relation light up as a ridge of power.
- **Assumes:** the simulation already ran; this is post-processing.
- **Method & tools:** NumPy 2-D FFT; overlay the F1/F2 analytic curves on the simulated ω–k power spectrum.
- **You'll learn:** the deeply satisfying moment when a self-consistent simulation reproduces the theoretical dispersion you derived independently — closing the loop on the whole playground.
- **Validation:** the simulated ω–k ridge falls on the analytic curve; thermal corrections (Bohm–Gross slope) are visible.
- **Compute:** seconds (post-processing); the sim itself is experiment 03/06.

### F4 — Full hot magnetized dispersion
- **Models:** the general hot, magnetized, multi-species dispersion relation at arbitrary angle — Bernstein modes, cyclotron harmonics, the works.
- **Assumes:** Maxwellian species, linear; full kinetic dielectric tensor with Bessel-function sums.
- **Method & tools:** implement the hot-plasma dielectric tensor, or drive an established solver (WHAMP/DSHARK-style); careful root tracking across parameter space.
- **You'll learn:** the complete linear wave physics of a magnetized plasma — the basis of RF heating and many instabilities.
- **Validation:** recover Bernstein-mode and cyclotron-harmonic structure; match a published dispersion diagram for a known regime.
- **Compute:** minutes; root-tracking is the hard part, not raw compute.

## Diagnostics you'll reuse
ω(k) dispersion curves, CMA diagram, complex-ω (damping/growth) plots, ω–k power spectra from simulation data with analytic overlays.

## Key references
- Stix, *Waves in Plasmas* — the definitive reference.
- Swanson, *Plasma Waves*. Gurnett & Bhattacharjee, *Introduction to Plasma Physics*.
- PlasmaPy `plasmapy.dispersion` docs; Fried & Conte, *The Plasma Dispersion Function*.

## Stretch goals
- Build one figure overlaying cold (F1), kinetic (F2), and simulated (F3) dispersion for the same wave.
- Identify which wave heats a tokamak at the ion-cyclotron resonance.
- Reuse the ω–k FFT tool on the experiment-06 Orszag–Tang run to pull out all three MHD wave speeds at once.

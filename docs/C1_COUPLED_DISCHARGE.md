# C1 — The event-coupled discharge (transport burn + sawteeth)

The headline rung (NIGHT.md Track C): the two halves of the project — the **transport
burn** (Track A, exp 09, seconds) and the **MHD instabilities** (Track B, exp 10,
microseconds) — finally run *together*. A sawtooth crash fires **during** the burning
discharge and redistributes the profiles, the laptop-scale cousin of a JOREK-style
"flight simulator".

## The staged two-timescale coupling (named honestly)

Transport evolves on the energy-confinement time τ_E (~seconds); the m=1 sawtooth
reconnection happens on the Alfvén/resistive-kink time, **instantaneous** on the
transport scale (τ_E/τ_A ~ 10⁶). So we do **not** march one monolithic code — we
stage them: step the transport burn, and when an MHD threshold trips, apply the MHD
event as an instantaneous profile redistribution, then continue transport. This is
the physically-correct reduction; fully self-consistent extended MHD across a whole
discharge is JOREK/NIMROD territory (named at the boundary).

## How it is wired (reusing validated pieces)

- The burn is the F2 `transport.Transport1D` (now with the A4 soft β-limit, so the
  core sits at a realistic ~25 keV).
- The MHD trigger is `sawtooth.q_from_temperature` — a **Spitzer-ohmic q-profile**
  (current J ∝ T^{3/2} peaks in the hot core, lowering q(0)). As the burn peaks the
  core, q(0) falls; when it crosses the kink threshold, `sawtooth.sawtooth_event`
  fires a Kadomtsev crash (`crash_profiles`) that flattens n and T inside the mixing
  radius, **conserving both particles and thermal energy**, then transport resumes.

No new instability physics — Track C is the *integration* of the B3a reconnection
operator into the Track-A burn.

## Validation (falsifiable — `tests/test_sawtooth.py`, the 3 Track-C tests)

- **`q_from_temperature`**: a peaked (hot-core) T gives a lower q(0) than a flat one,
  and q rises monotonically to the prescribed edge value — the coupling that lets the
  burn drive its own sawteeth.
- **`crash_profiles`**: a crash conserves the particle content **and** the thermal
  energy exactly (to 1e-10) while flattening n and T inside the mixing radius.
- **`sawtooth_event`**: fires on a hot, kink-unstable core (q(0) < trigger) and leaves
  a cool, kink-stable one **untouched** — the "events-off recovers Track-A" regression
  in kernel form.

And in the driver (`run.py --mode coupled`, printed): the burning core drives q(0) to
~0.93 → **179 sawtooth crashes** over the discharge; the steady H-mode core T0
sawtooths between **23–26 keV**; energy is conserved across every crash to grid
resolution (~5e-3 at the mixing boundary); and **turning events off gives 0 crashes
and the pure Track-A burn** (the sawteeth shift the core by up to ~16 keV — the
coupling is real, not a no-op).

## Deliverable

`outputs/tokamak_discharge_full.gif` — the integrated discharge: the poloidal
cross-section T(ρ,t) beside the core-T0 sawtooth trace and q(0) crossing the kink
threshold, through **ignition → burning H-mode with periodic sawteeth → pellet fuel
injection → settling**. Regen: `python gif_gallery.py tokamak_discharge_full`;
`python experiments/09_burning_plasma/run.py --mode coupled --save` for the still.

## Scope boundary (stated honestly)

- **Staged, not self-consistent.** The MHD event is an instantaneous Kadomtsev
  redistribution applied to the transport profiles; there is no in-situ 2-D MHD solve
  during the burn. This is the correct two-timescale reduction, not full extended MHD.
- **The q-profile is a model** (Spitzer current from the instantaneous T), not a
  self-consistent current-diffusion equation; it sets *when* and *where* the crash
  fires, which is what the coupling needs. The crash trigger carries a small margin
  (q(0) < 0.93) so crashes are finite-amplitude and well separated (as in B3).
- **Energy conservation across a crash is ~5e-3** (full-domain), limited by the
  mixing-boundary grid cell; the kernel conserves the inside-region integral exactly.
- Only the **sawtooth** event is wired; a tearing/island event (Δ′>0 → island
  flattening) is the natural next coupling (C1b), and the 3-D torus render is C2.

## References
- NIGHT.md Track C — the staged two-timescale coupling charter.
- Kadomtsev, *Sov. J. Plasma Phys.* 1, 389 (1975) — the sawtooth crash.
- Hawryluk, in *Physics of Plasmas Close to Thermonuclear Conditions* (1980) — the
  transport-with-events modelling approach (ASTRA/TRANSP lineage).

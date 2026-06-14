# NIGHT.md — autonomous overnight work charter

> **You (Claude) are running unattended overnight with permissions skipped.** This
> file is your directive. Read it fully, then read `CLAUDE.md`, `docs/FIDELITY.md`,
> `docs/3D_TOKAMAK_GUIDE.md`, `experiments/09_burning_plasma/PLAN.md`, and the two
> memory files (`tokamak-3d-ladder-progress`, `burning-plasma-transport`). Then
> resume from `NIGHT_LOG.md` (create it if absent) and keep climbing the ladder
> below until morning. Work continuously. Leave the repo green at every commit.

---

## 1. The mission (in the user's words, corrected)

The user — new to plasma physics, deferring to your expertise — asked for:

> "an MHD/CFD sim of the tokamak: ignition, steady state, fuel injection, and other
> instabilities … each iteration with more and more fidelity … ideally a `.gif` for
> each operation mode showcasing this 3-D tokamak sim. Fill in the gaps. If you
> finish early, extend the workflow to a 3-D stellarator."

**The honest reframing (state this in your memos so it stays correct):** what they
want is **two coupled simulations on two timescales**, not one CFD run:

| The user said | What it physically is | Timescale | Where it lives |
|---|---|---|---|
| ignition / steady state / fuel injection | **transport** — energy & particle balance | seconds (τ_E) | exp 09 (`plasmaplay/transport.py`) |
| "other instabilities" (sawtooth, tearing, kink, ELM) | **MHD/CFD** — plasma *fluid motion* | µs–ms (τ_A) | exp 06 `fvm.py` + T4 `tearing.py` |
| "3-D tokamak" | the field **geometry** that hosts both | static | T0–T3 (`plasmaplay/tokamak.py`) |

MHD ≈ CFD + magnetic field, so "MHD/CFD" is fair for the instability side. The burn
arc is transport, not CFD. The grand deliverable is a **staged coupling**: run the
transport burn, and let MHD events (a sawtooth crash, a tearing island) fire *during*
it and redistribute the profiles — the laptop-scale cousin of a JOREK "flight
simulator." Fully self-consistent 3-D extended MHD is the frontier; **do not pretend
to do it** — build the reduced, *validated* version and name the boundary.

---

## 2. Iron rules (you are unattended — these are not optional)

1. **Validation-first, always.** Every kernel gets a falsifiable test (a number it
   must hit) *before* an experiment relies on it. "The plot looks plasma-ish" is a
   bug, not a result. If you cannot validate a rung, **do not ship it as working** —
   write a `docs/` memo explaining why and move on.
2. **Never lie in a result.** If a check fails, report the failure in the log and
   memo. Faking or hand-waving a number is the worst possible outcome of the night.
3. **Keep the tree green.** Never commit with failing `pytest` or `ruff check .`
   errors. Run both before every commit.
4. **Commit locally after each green milestone; never push.** Use the message
   convention in `CLAUDE.md`. Branch off `master` is fine; staying on `master` with
   local commits is also fine (the user's established pattern). Do **not** `git push`,
   force-push, rebase published history, or touch the remote.
5. **Stay in the repo + the memory dir.** Do not modify anything outside
   `/Users/julianwang/plasma-playground/` except the auto-memory directory. No
   network installs unless a rung explicitly needs one and it's pure-Python and
   safe (`pip install` into the `.venv` only — prefer NOT to; reach for it only for
   `imageio`/`pillow` if gif-writing needs it).
6. **No destructive shell.** No `rm -rf`, no overwriting files you didn't create
   without reading them first, no `git clean -fdx`. `outputs/` is gitignored — write
   gifs there freely.
7. **Reuse, don't duplicate.** Shared logic goes up into `plasmaplay/`. Experiments
   are thin drivers. Honor the dependency web in `CLAUDE.md`.
8. **Checkpoint relentlessly** (see §6) so a fresh context can resume seamlessly.

---

## 3. Environment quick-reference

```bash
source .venv/bin/activate                 # ALWAYS; base conda is active by default
python -m pytest -q                       # full suite (~100 s); keep it green
ruff check .                              # lint: line length 100, py311
MPLBACKEND=Agg python experiments/09_burning_plasma/run.py --save   # headless figs
```

Gifs: use `matplotlib.animation.FuncAnimation` + `PillowWriter` (pure Python, no
external binary). If Pillow is missing, `uv pip install pillow` into the `.venv`.
Keep gifs sane: ~80–120 frames, dpi ≈ 90, target < ~8 MB each. Save to the
experiment's `outputs/` (gitignored). Also save a representative **PNG montage**
(committable is fine — but note outputs/ is gitignored, so montages for the README
go in the experiment folder root or a `docs/figs/` you create and un-ignore if you
want them tracked; default: keep everything in outputs/ and reference by path).

---

## 4. Architecture & reuse map (build on what exists)

Already built and validated — **reuse, do not rewrite**:
- `plasmaplay/transport.py` — Bosch-Hale `reactivity_dt`, `fusion_power_density`,
  `bremsstrahlung_density`, `burn_0d` (F0), `Transport1D` (F2), `gaussian_deposition`.
- `plasmaplay/tokamak.py` — `equilibrium_field`, `toroidal_poincare`,
  `safety_factor`, `helical_perturbation`, `superpose` (T0–T3).
- `plasmaplay/tearing.py` — slab linear tearing: `delta_prime_*`, `tearing_growth_rate`.
- `plasmaplay/solvers.py` — `grad_shafranov_solve` (the 2-D equilibrium), Poisson.
- `plasmaplay/fvm.py` — 1-D ideal-MHD Godunov (HLL); `diagnostics.py` — FFT/Poincaré.

New modules you will likely add (validate each): `plasmaplay/animate.py` (gif
helpers), `plasmaplay/reduced_mhd.py` (2-D Strauss reduced MHD), and extensions to
`transport.py` (ash, two-temperature, equilibrium coupling).

---

## 5. THE LADDER — climb in order; each rung has a validation gate and a deliverable

> Pacing: each rung = kernel+test → experiment wiring → validation print → figure →
> **gif** → update PLAN/README → `docs/` memo → local commit → log entry. Do them
> one at a time and commit green. The gif is part of "done," not optional.

### Track G — Visualization foundation (DO THIS FIRST; everything else needs it)
- **G1. `plasmaplay/animate.py`.** Helpers: `animate_profiles(times, fields, ...)`
  (radial profiles over time → gif), `animate_cross_section(times, profile_t, ...)`
  (poloidal disk heatmap over time → gif), `animate_torus_3d(...)` (rotating 3-D
  flux-surface / field render). All return/save gifs via PillowWriter.
  - **Validation:** a deterministic synthetic field (e.g. a diffusing Gaussian)
    produces a gif whose frame count, shape, and final-frame values match an
    analytic expectation; unit-test the frame-data builder (not the gif bytes).
  - **Deliverable:** `outputs/_smoke_diffusion.gif` (proof the pipeline works).

### Track A — Burning-plasma transport: climb fidelity (exp 09)
- **A1 (F1). He-ash + dilution + β-limit.** Add a helium-ash species produced by
  fusion (`dn_He/dt = R_fus − n_He/τ_He*`), fuel dilution (n_DT = n_e − 2n_He),
  and a soft β-limit that caps the operating point. *Gap you're filling:* this is
  why real machines need ash pumping and continuous fuelling, and why the burn
  point sits at ~15–25 keV not ~80 keV.
  - **Validation:** steady ash fraction f_He ≈ τ_He*/τ_E·(burn rate) prediction;
    Q falls as dilution rises; with the β-limit the 0-D operating point lands in
    the 10–25 keV burning band.
  - **Deliverable gif:** `burn_0d_ignition.gif` — phase-space (n,T) trajectory
    igniting onto the POPCON operating point, with ash building up.
- **A2 (F2). Two-temperature (Te, Ti) + heating mix.** Split into electron and ion
  channels with collisional equipartition (Spitzer ν_ei); alpha power to electrons,
  NBI to ions, RF selectable; separate χ_e, χ_i. *Gap:* real plasmas run Ti≠Te.
  - **Validation:** equipartition time matches the Spitzer formula; relaxes to
    Te=Ti when collisional; ITER-like Ti/Te ratio under alpha+NBI heating.
  - **Deliverable gif:** `burn_1d_two_temperature.gif` — Te(ρ,t) and Ti(ρ,t) panels.
- **A3 (F3). Transport on the REAL equilibrium.** Replace the circular ρ with
  flux-surface-averaged coordinates from `solvers.grad_shafranov_solve` (exp 04):
  V'(ψ), ⟨|∇ρ|²⟩ metrics, a D-shaped boundary. Re-solve the equilibrium as pressure
  evolves (Picard outer loop). *This is the bridge that makes it a "real machine."*
  - **Reuses:** exp 04 GS solver, `tokamak.py` field bridge.
  - **Validation:** flux geometry matches the equilibrium; τ_E consistent with an
    ITER98 / IPB98(y,2) H-mode scaling for the chosen Ip, B, P, n, a, R.
  - **Deliverable gif:** `burn_dshaped_cross_section.gif` — the **D-shaped**
    cross-section igniting and burning (the headline "watch it burn" movie).
- **A4 (F3.5). Operational limits & modes.** Add the Greenwald density limit, the
  L→H transition (a confinement bifurcation), and a **radiative-collapse failure
  mode** (over-fuel past n_G → radiation > heating → the burn dies). *Gap:* shows
  the operating *window*, not just one happy path.
  - **Validation:** density limit at n_G = Ip/(πa²); L-H power threshold scaling;
    collapse is reversible by backing off the fuelling.
  - **Deliverable gif:** `operating_modes.gif` — a scan/montage of L-mode, H-mode,
    and a disruption-by-over-fuelling.

### Track B — MHD instabilities: the fluid-motion half (exp 06 / T4 → new)
- **B1. Cylindrical linear MHD stability.** Lift T4 from slab to a periodic
  cylinder ("straight tokamak"): the m/n internal kink and tearing on a real
  q(r) profile (reuse the T1/T3 q-profile). Solve the Newcomb/tearing eigenvalue.
  - **Validation:** m=1 internal kink unstable when q(0)<1 (the sawtooth trigger);
    Δ′(q=m/n surface) sign predicts tearing stability; growth rate γ∝S^(−3/5) holds
    on the cylinder as it did in the slab.
  - **Deliverable gif:** `kink_eigenmode.gif` — the unstable eigenfunction.
- **B2. Nonlinear 2-D reduced MHD (`plasmaplay/reduced_mhd.py`).** Implement the
  Strauss reduced-MHD equations — evolve poloidal flux ψ and vorticity U (stream
  function φ) in (r,θ) for a chosen toroidal n — with resistivity. Watch a tearing
  mode grow and **saturate** into a magnetic island (Rutherford regime).
  - **Reuses:** `solvers` elliptic inversion for φ from U; FFT diagnostics for mode
    growth; T3 Poincaré to *see* the island.
  - **Validation:** linear phase matches B1 growth rate; nonlinear island width
    follows the Rutherford equation dW/dt ∝ Δ′(W) and saturates at Δ′(W_sat)=0.
  - **Deliverable gif:** `tearing_island_saturation.gif` — flux contours tearing
    into an island chain and saturating.
- **B3. The sawtooth cycle (Kadomtsev).** Using B1/B2: when q(0)<1, trigger an
  m=1 reconnection that flattens the core, conserving helical flux; q(0) relaxes
  back >1, then resistive diffusion re-peaks the current and it repeats.
  - **Validation:** reconnection conserves helical flux (to numerical tol); the
    crash flattens T inside q=1; sawtooth period scales with resistive time.
  - **Deliverable gif:** `sawtooth_cycle.gif` — the core temperature crashing and
    recovering in a sawtooth.

### Track C — The coupling (the integrated dream movie)
- **C1. Event-coupled discharge.** Run the Track-A transport burn; monitor q(0)
  and Δ′ from the evolving profiles; when an MHD threshold trips (q(0)<1 → sawtooth;
  Δ′>0 at a rational surface → tearing), apply the corresponding profile
  redistribution (Kadomtsev flatten / island flattening), then continue transport.
  *This is the staged two-timescale coupling — name it honestly as such.*
  - **Validation:** energy/particle conserved across each event (minus the
    physically-correct redistribution); sawtooth-period and island onset track the
    inputs sensibly; turning events off recovers pure Track-A.
  - **Deliverable gif:** **`tokamak_discharge_full.gif`** — the headline:
    ignition → burning H-mode with periodic sawteeth → pellet fuel injection →
    a tearing island appearing → settling. Cross-section + time-trace, side by side.
- **C2. 3-D showcase.** Render the discharge on the 3-D torus: flux surfaces colored
  by T(ρ,t), with a magnetic island / sawtooth visible, rotating.
  - **Deliverable gif:** `tokamak_3d_discharge.gif`.

### Track E — Stellarator (only if A–C are solidly done; stretch)
- **E1.** Build/confirm a 3-D stellarator vacuum field (exp 05 has field lines +
  Poincaré). Make a clean `stellarator_field` callable and a flux-surface extractor.
  - **Validation:** nested flux surfaces in the Poincaré; rotational transform ι
    from external coils only (stellarators carry ~no net plasma current).
  - **Deliverable gif:** `stellarator_flux_surfaces.gif` — the twisty surfaces rotating.
- **E2.** Run the Track-A transport model on stellarator flux surfaces. *Gaps to
  teach:* stellarators are **inherently steady-state** (no current drive, no
  disruptions, no sawteeth — contrast with the tokamak!), so the "operation modes"
  collapse to startup → steady burn → fuelling; the relevant instability physics is
  different (no current-driven kinks; neoclassical/turbulent transport dominates).
  - **Validation:** steady burn with zero loop voltage; confinement vs the tokamak.
  - **Deliverable gif:** `stellarator_burn.gif`.

> **If you finish *everything*:** deepen validation (convergence studies, compare
> growth rates / scalings against published numbers in the memos), optimize the
> proven hot loops (field-line tracer batching, reduced-MHD step) with numba, and
> assemble a single `gallery.py` that regenerates every gif. Then stop and write a
> closing summary in `NIGHT_LOG.md`.

---

## 6. Checkpointing & resumption protocol (critical for an unattended run)

You **will** lose context across the night. The survival mechanism is **NOT
self-compaction** — you cannot read your own context gauge, and you cannot reset
your own window on demand. Only the harness can compact, and it fires automatically
near the limit (~92–95%), not at any percentage you choose. So **do not try to
"compact at 50–70%"** — you can't detect that, and any plan that depends on it will
silently fail. Instead, make context loss a **non-event** by checkpointing to disk
at every milestone boundary (a real, detectable event), so resumption is
deterministic no matter when — or whether — the window resets.

**Checkpoint at the END of every rung, BEFORE starting the next** (never mid-rung):

1. **`NIGHT_LOG.md`** (repo root) — append-only. Write:
   `## [rung id] — DONE/PARTIAL/FAILED — <commit hash>` then 3–6 bullets: what was
   built, the validation number achieved (or why it failed), the gif path, and the
   single most useful gotcha. This is the first thing the next context reads.
2. **Memory** — update `burning-plasma-transport` (and add new memory files for
   reduced-MHD / coupling / stellarator as those start) with durable facts and
   gotchas, per the memory rules in the system prompt. Keep `MEMORY.md` index current.
3. **Repo memos** — `docs/<NAME>.md` per substantial rung, justifying design
   choices and recording the validation, matching the existing `docs/T*_*.md` voice.
4. **Commit locally** (green only). The commit + log entry together are the
   checkpoint. After this point, the repo + log + memory are the source of truth —
   **assume your conversation history may vanish at any moment.**

**Keep each turn lean** so you reach milestone boundaries before the window fills:
re-read files from disk instead of holding large contents in context; delegate
broad searches to sub-agents (their file dumps don't land in your context, only
their conclusions); don't paste big outputs you can regenerate.

**When the harness auto-compacts** (or a fresh session starts), resume:
`NIGHT_LOG.md` last entry → `git log --oneline -15` → confirm `pytest -q` green →
pick up at the next unstarted rung. **Never redo a DONE rung.**

**Strongly preferred runner: a fresh-session loop.** The most robust way to run the
whole night is one milestone per wake-up, each starting from a near-empty context
(read log → do ONE rung → checkpoint → schedule next). This sidesteps context
filling entirely instead of managing it. If launched that way, do exactly one rung
per invocation and stop after the checkpoint. If instead running as one long
session, rely on the per-milestone checkpoints above and let harness auto-compact
handle the window.

---

## 7. Definition of done (per rung — all must hold before you commit)

- [ ] Kernel code in `plasmaplay/` (pure, SI, docstring explains the *physics*).
- [ ] A test in `tests/` that asserts a **specific number / scaling** (falsifiable).
- [ ] Experiment `run.py` wiring that **prints the validation result** when run.
- [ ] The **gif** (and a PNG still) generated in `outputs/`.
- [ ] `pytest -q` fully green **and** `ruff check .` clean.
- [ ] `PLAN.md` status marker + `README.md` roadmap updated.
- [ ] A `docs/` memo for non-trivial rungs.
- [ ] Local commit with a descriptive message; `NIGHT_LOG.md` + memory updated.

---

## 8. Physics gaps I'm filling for you (so the sim is honest, not just pretty)

- **Q, breakeven, ignition.** Q = P_fusion/P_aux. Breakeven Q=1, "burning plasma"
  Q≳5, ignition Q=∞ (no external heating needed). The night's burn should report Q.
- **The Lawson triple product** nTτ_E ≳ 3×10²¹ keV·s·m⁻³ for DT ignition — the
  single number that decides "does it light."
- **Greenwald density limit** n_G = Ip/(πa²): fuel past it and you disrupt. Track A4.
- **L–H transition:** above a power threshold, confinement bifurcates upward (the
  H-mode pedestal). Real machines live in H-mode. Track A4.
- **Two-timescale separation** τ_E/τ_A ~ 10⁶ — *why* transport and MHD are separate
  codes, and why the coupling in Track C is staged, not monolithic.
- **Stellarator contrast:** no plasma current ⇒ no disruptions, no sawteeth,
  inherently steady-state; the price is 3-D-shaped coils and neoclassical transport.
- **What we are NOT doing (say so):** self-consistent 3-D extended MHD, gyrokinetic
  turbulence, real RMP/ELM control, disruption runaway electrons. Those are
  JOREK/NIMROD/GENE/XGC — name them at the scope boundary in the memos.

---

## 9. The deliverable the user actually wants

A set of **`.gif`s, one per operation mode**, showcasing the (reduced-D, validated,
3-D-rendered) tokamak: **ignition**, **steady burn (H-mode)**, **fuel injection
(pellet)**, **sawtoothing**, **tearing/islands**, and the **full coupled discharge**
— plus a 3-D torus render. Then the **stellarator** equivalents if time remains.
Every gif must be backed by a passing validation test; a pretty gif with no check
behind it is a failure by the rules of this repo.

Good night. Climb carefully, validate everything, commit green, and log as you go.

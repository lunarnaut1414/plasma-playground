# REVIEW.md — gif-gallery review & polish charter (loop one rung per wake)

A self-paced review loop, in the spirit of `NIGHT.md`, but pointed *backward* at what we
already built. The job: **look at every output gif, judge it on four axes, and make it
correct, well-paced, legible, and pretty** — without ever faking a number or breaking the
green tree.

Resume from `REVIEW_LOG.md`. One rung (≈ one gif) per wake-up. Checkpoint at every rung.

---

## Why this loop exists (the user's concerns, verbatim intent)

1. **"Some gifs don't look right."** — physics or rendering artifacts; the picture
   disagrees with the validated number.
2. **"Too short of a sim time."** — the *physical* duration `t_end` ends before the
   phenomenon finishes (island still growing, burn not settled, only one rotation).
3. **"The timestep looks too short."** — pacing: too few saved frames, too-fast fps, or a
   **frame stride that aliases fast events** (e.g. ~179 sawtooth crashes sampled onto 100
   frames → flicker). Distinguish the *numerical* `dt`, the *physical* `t_end`, the
   *frame stride*, and the *playback fps* — they are four different knobs.
4. **"I want the gifs to look pretty and eye-catching."** — aesthetics is a real
   deliverable here, not gold-plating. Pretty *and* honest.

---

## The four review axes (score each gif 1–5 on all four)

| Axis | The question | Red flags |
|------|--------------|-----------|
| **A. Correct** | Does the picture match the validated test/number? | profile non-monotone where physics says monotone; energy/particles visibly not conserved across a crash; ι twist absent; colorbar implies values the sim never reached |
| **B. Paced** | Is `t_end` long enough, and does playback show the whole story at a watchable rate? | ends mid-transient; fast events alias/flicker (stride > event spacing); >half the frames are a static hold; <3 s or >12 s total playback |
| **C. Legible** | Can a newcomer read it in 5 s? | unlabeled axes; no units; no time readout; no annotation of *what* moment (ignition, crash, pellet); cramped layout |
| **D. Pretty** | Is it eye-catching and polished? | default matplotlib look; jet/clashing cmap; muddy background; jaggy/low dpi; no easing on the camera; abrupt loop seam; chart-junk |

**Verdict per gif:** `KEEP` (≥4 on all axes) · `TWEAK` (one axis weak, cheap fix) ·
`REWORK` (≥2 axes weak, needs real changes) · `DEFER` (real but out of scope — log it).

---

## Iron rules (inherited from NIGHT.md, plus review-specific)

- **Never fake a number, never fake a frame.** If you prettify, the underlying sim must
  still produce the same validated quantities. Cosmetic ≠ dishonest: a *schematic*
  animation (e.g. a hand-cranked amplitude) is allowed **only if its title/caption says
  so**. If a gif implies a time-integration it isn't doing, that's an A-axis defect.
- **Keep the tree green.** After any code change: `pytest` and `ruff check .` must pass
  before you commit. Gif code lives in `gif_gallery.py` / `plasmaplay/animate.py`.
- **Validation-first still holds.** If a fix changes the *physics* path (not just style),
  there must be a test pinning the new behavior. Pure style changes (cmap, dpi, labels,
  fps, easing) don't need a new test but must not break existing ones.
- **`outputs/` is gitignored** — never commit gifs/pngs. Commit the *code* that makes
  them. Regenerate to verify.
- **Commit per rung** once green. Message ends with the Co-Authored-By trailer. **Never
  push** unless the user explicitly asks.
- **One rung per wake. CRITICAL: do NOT end the turn with only prose — the loop continues
  ONLY if you actually call ScheduleWakeup with the same prompt.** Forgetting it stops the
  loop. Re-arm every turn until the queue is done.
- **Don't end on a broken gif.** If a rework is mid-flight at rung's end, either finish it
  or revert to the last green state before checkpointing.

---

## How to actually LOOK at a gif (you can't watch it — sample frames and Read them)

The Read tool renders PNGs visually. So dump representative frames and read them:

```python
# scripts/_dump_frames.py  (write once in R0; reuse every rung)
import sys, os
from PIL import Image, ImageSequence
g = sys.argv[1]                                  # e.g. outputs/tearing_island_saturation.gif
im = Image.open(g); frames = list(ImageSequence.Iterator(im))
n = len(frames); name = os.path.splitext(os.path.basename(g))[0]
out = "outputs/_review"; os.makedirs(out, exist_ok=True)
picks = sorted(set([0, n//8, n//4, n//2, 3*n//4, 7*n//8, n-1]))   # weight the ends
for f in picks:
    frames[f].convert("RGB").save(f"{out}/{name}_{f:03d}.png")
print(name, "n_frames=", n, "saved", picks)
```

Run `python scripts/_dump_frames.py outputs/<name>.gif`, then **Read each
`outputs/_review/<name>_NNN.png`** and judge axes A–D from what you actually see. To check
pacing/aliasing, also compute event spacing vs stride (see per-gif seeds below).

---

## House visual style (the "pretty" target — apply consistently)

Make the gallery look like one designed set, not 11 defaults. Put shared helpers in
`plasmaplay/animate.py` (e.g. `apply_house_style()`), don't copy-paste.

- **Colormaps (perceptually uniform, on-theme):** temperature → `inferno`/`magma`;
  flux/signed fields → `RdBu_r`; density/neutral → `viridis`/`cividis`. Never `jet`.
  Keep one cmap per physical quantity across all gifs.
- **Background & frame:** dark figure face (`#0e1116`) with light text reads as "plasma"
  and makes hot cmaps pop; or a clean white set — pick ONE and apply to all. Kill the
  default gray spines; thin them or remove.
- **Type:** a consistent title + a persistent **time/parameter readout** (`t = 12.3 s`),
  axis labels with **units**, and a one-line caption of the physics moment. Legible font
  sizes (≥10). Annotate key events inline (ignition band, sawtooth crash tick, pellet
  arrow).
- **Motion:** ease the camera/parameter (smoothstep, not linear) and **make it loop
  seamlessly** (return to start, or ping-pong) so there's no jarring snap. Hold the final
  state ~0.5 s before looping.
- **Crispness:** dpi ≥ 110 for the hero gifs (3-D discharge, tearing, stellarator);
  anti-alias lines (`lw≥1.2`, `alpha` for depth). Watch file size — keep heroes < ~6 MB
  (trim frames or dpi if needed).
- **Pacing default:** aim **5–10 s** total playback. fps 16–24 for smooth motion; lower
  only for deliberately stepwise stories. **Frame stride must be ≤ the spacing of the
  fastest event you want visible** — otherwise crashes/oscillations alias.

---

## The review queue (rungs)

Do R0 first, then one gif per rung (heroes first — they carry the gallery). R-final last.

- **R0 — setup:** write `scripts/_dump_frames.py`; (re)generate all 11 gifs fresh
  (`MPLBACKEND=Agg python gif_gallery.py all`); record the validation line each prints
  into `REVIEW_LOG.md`; confirm `pytest`/`ruff` green as a baseline.
- **R1 — `tokamak_3d_discharge`** (hero; suspected sawtooth aliasing).
- **R2 — `tearing_island_saturation`** (hero; suspected slow-start / late growth).
- **R3 — `stellarator_flux_surfaces`** (hero; pure camera spin — easing + seamless loop).
- **R4 — `tokamak_discharge_full`** (hero; the 2-panel flight-sim).
- **R5 — `stellarator_burn`** (the contrast 2-panel).
- **R6 — `kink_eigenmode`** (check: schematic amplitude ramp — is it labeled as such?).
- **R7 — `operating_modes`** (45 s sim — watchable? L/H/disruption legible?).
- **R8 — `burn_dshaped_cross_section`**.
- **R9 — `burn_1d_two_temperature`**.
- **R10 — `burn_0d_ignition`**.
- **R11 — `_smoke_diffusion`** (the smoke test — lowest priority).
- **R-final — synthesis:** apply `apply_house_style()` consistently across any gif that
  missed it; regenerate `all`; update `gallery.png` montage if present; write the
  REVIEW_LOG closing summary (scores table + what changed + any DEFERs). Confirm green.

Scale effort to verdict: a `KEEP` rung is just "looked, scored, logged." Spend the budget
on `REWORK`s and on the shared house-style helper.

---

## Per-gif seed notes (grounded in current params — verify, don't trust blindly)

All saved at dpi=90 currently (→ bump heroes to ≥110). Frame counts confirmed from the
files. fps from `gif_gallery.py`.

- **`tokamak_3d_discharge`** — `dt=2e-3, t_end=22 s`, **stride = steps/100 → 100 frames**,
  fps=14. Sim fires **~179 sawtooth crashes**. **100 frames ≪ 179 crashes ⇒ crashes alias
  → flicker.** Fix options: (a) raise saved frames so stride ≤ crash spacing, (b) render a
  short *zoomed* window (e.g. 8–13 s) at fine stride so individual crashes are visible,
  (c) carry a "crash happened since last frame" flag and flash a marker so the eye reads
  the event instead of aliasing. Bump dpi, dark theme, add `t=` readout + pellet arrow.
- **`tearing_island_saturation`** — `dt=0.012, t_end=300 τ_A`, 91 frames @16. Island is
  resistive/slow; **likely flat for the first half then grows late** → looks static then
  sudden. Fix: trim `t_end` to where W actually turns over, or **non-uniform sampling**
  (denser frames through the growth knee). Annotate the dW/dt-peak frame. RdBu_r is good;
  add island-width marker synced to the right panel.
- **`stellarator_flux_surfaces`** — 72 frames @16, **geometry static, only the camera
  azimuths 0→360 linearly.** Make the spin **ease + loop seamlessly** (azim already wraps
  0→360, good — verify no seam), maybe add a slow elevation bob. Richer line shading
  (depth alpha), dark bg → the colored surfaces pop. This is the gallery's beauty shot.
- **`tokamak_discharge_full`** — `t_end=22 s`, 110 frames @14, 2 panels. Check the burn
  has *settled* by 22 s (if T0 still drifting, extend). Sawtooth ticks on the T0 trace
  should be visible, not aliased (same stride concern as R1 but on a 1-D trace it reads
  better). Label ignition band / pellet line (already has axvspan/axvline — style them).
- **`stellarator_burn`** — 100 frames @14, 2 panels (stellarator cross-section + T0
  contrast). Cross-section is `broadcast_to` a 1-D profile → concentric rings (correct,
  but flat-looking). Could add subtle l=2 elliptical shaping to *look* like a stellarator
  while labeling it as illustrative. Verify the "0 vs N sawteeth" headline still prints.
- **`kink_eigenmode`** — 90 frames @18. **Amplitude is `0.33*(linspace^1.5)` and phase
  only sweeps 1.5π — a hand-cranked schematic, NOT a γ-driven time-integration.** A-axis
  call: either (a) re-label clearly as "schematic m=1 displacement (illustrative)", or
  (b) drive the amplitude from the actual linear growth rate so it's physical. Pretty: the
  crescent on `inferno` is already nice; full 2π phase would loop seamlessly.
- **`operating_modes`** — `t_end=45 s`, 111 frames @20 → ~5.5 s playback (good). Verify
  L→H→disruption are each annotated and the Greenwald/β lines are labeled.
- **`burn_dshaped_cross_section`** — 100 frames @20. D-shape should be visibly D-shaped;
  inferno temperature; add `t=` + Q or fusion-power readout.
- **`burn_1d_two_temperature`** — 100 frames @20. Two curves Te/Ti — make the
  ion-hotter-than-electron gap unmistakable (fill_between the gap, label equipartition).
- **`burn_0d_ignition`** — 101 frames @20. 0-D ignition; show the Lawson/operating-point
  context, not just a rising line.
- **`_smoke_diffusion`** — 90 frames @20. Internal smoke test; lowest priority, but a
  clean Gaussian-spreading look is a cheap win.

---

## Fix protocol (per rung)

1. `python scripts/_dump_frames.py outputs/<name>.gif`; Read the sampled PNGs.
2. Score A–D; write the verdict + the specific defect(s) into `REVIEW_LOG.md`.
3. If `TWEAK`/`REWORK`: edit `gif_gallery.py` (and/or `plasmaplay/animate.py` for shared
   style). Prefer pushing reusable polish into `animate.py` (`apply_house_style`, eased
   camera, seamless-loop helper, crash-flash marker) so every gif benefits.
4. If the change touches the **physics/data path**, add/adjust a test so the validated
   number is still pinned. Pure-style changes: no new test, but run the suite.
5. Regenerate just that gif (`MPLBACKEND=Agg python gif_gallery.py <name>`), re-dump,
   re-Read → confirm the defect is gone and nothing regressed visually.
6. `pytest` + `ruff check .` green → commit (code only; gifs stay gitignored).
7. **Call ScheduleWakeup with this same prompt to advance to the next rung.**

---

## Logging format (`REVIEW_LOG.md`, append-only)

Per rung:

```
## R<k> — <gif_name> — <KEEP|TWEAK|REWORK|DEFER> — <commit hash or "no code change">
Scores: A<n> B<n> C<n> D<n>
Saw: <what the sampled frames actually showed>
Defect: <the specific problem, with the number — e.g. "179 crashes / 100 frames → alias">
Did: <what changed, or "scored only">
Verify: <regenerated; pytest N passed; ruff clean>
```

End with a **CLOSING SUMMARY**: the 11×4 score table, total gifs reworked, the shared
house-style helper added, file-size/dpi deltas, and any `DEFER`red items for a later loop.

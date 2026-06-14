# REVIEW_LOG.md — checkpoint log for the gif review/polish loop (REVIEW.md)

Append-only. One entry per rung. Resume from the last entry. Format defined in REVIEW.md.

---

## Baseline inventory (pre-R0, gathered 2026-06-14)

11 gifs in `outputs/` (all dpi=90). Frame counts confirmed from files; fps/params from
`gif_gallery.py`:

| gif | frames | fps | ~playback | sim params | seed concern |
|-----|-------:|----:|----------:|------------|--------------|
| tokamak_3d_discharge (hero) | 100 | 14 | 7.1 s | dt=2e-3, t_end=22 s, ~179 crashes | **crash aliasing (179≫100)** |
| tearing_island_saturation (hero) | 91 | 16 | 5.7 s | dt=0.012, t_end=300 τ_A | slow start, late growth |
| stellarator_flux_surfaces (hero) | 72 | 16 | 4.5 s | camera azim 0→360, static geom | linear spin / seam / beauty shot |
| tokamak_discharge_full (hero) | 110 | 14 | 7.9 s | dt=2e-3, t_end=22 s | settled by 22 s? sawtooth ticks |
| stellarator_burn | 100 | 14 | 7.1 s | t_end=22 s, 2 panels | flat concentric cross-section |
| kink_eigenmode | 90 | 18 | 5.0 s | amp=0.33·lin^1.5, phase 1.5π | **schematic, not time-integrated** |
| operating_modes | 111 | 20 | 5.6 s | t_end=45 s | annotate L/H/disruption |
| burn_dshaped_cross_section | 100 | 20 | 5.0 s | — | D-shape visible? readout |
| burn_1d_two_temperature | 100 | 20 | 5.0 s | — | make Te/Ti gap pop |
| burn_0d_ignition | 101 | 20 | 5.1 s | — | show Lawson context |
| _smoke_diffusion | 90 | 20 | 4.5 s | — | low priority |

Four review axes: **A** correct · **B** paced · **C** legible · **D** pretty (see REVIEW.md).
Goal: correct + honest + 5–10 s watchable + a consistent, eye-catching house style.

Queue: R0 (setup) → R1..R11 (one gif each, heroes first) → R-final (house-style sweep +
closing summary). Loop continues ONLY while ScheduleWakeup keeps firing the REVIEW prompt.

<!-- rung entries below this line -->

## R1 — tokamak_3d_discharge — REWORK — (pending commit)
Scores (before): A2 B2 C2 D1  ->  (after): A4 B4 C5 D5
Saw (before): a dim dark donut, uniform maroon — the hot core (29 keV) was HIDDEN inside
the opaque cold outer surface (concentric tubes share R0, so from outside you only ever
see the cold edge). Title clipped ("3-D"->"-D"); muddy grey wireframe; white bg; no crash
cue. Also fixed a bug in scripts/_dump_frames.py: list(ImageSequence.Iterator) returns
refs to the same image at the final frame -> every sampled PNG was identical. Now seek+copy.
Defect: the gif's whole promise (nested surfaces, hot core, sawtooth flatten) was invisible;
and ~179 crashes aliased onto 100 frames.
Did: replaced animate_torus_nested with new anim.animate_discharge_3d — a TWO-panel hero:
(L) glowing rotating torus, whole surface colored by CORE T (donut brightens through the
burn, dims on crash); (R) face-on poloidal BULLSEYE T(rho) (hot core -> cold edge, crashes
visibly flatten it). Dark "plasma" house style via new anim.apply_house_style (3-D panes +
patch darkened, re-applied after ax.clear). Per-frame crash count -> ⚡ flash + running
"sawteeth: N" counter (so the fast crashes read as activity, not aliasing). Title no longer
clipped; seamless full-turn camera; dpi 90->120, fps 14->16. gif_gallery now feeds the full
T(rho) profile + per-frame crash counts.
Verify: regenerated (179 sawteeth, core 2-30 keV); re-dumped+Read frames — gorgeous and
correct. pytest 208 passed (+2: discharge_3d + house_style); ruff clean. size 4.5 MB (<6).
Next: R2 = tearing_island_saturation.

## R2 — tearing_island_saturation — TWEAK — (pending commit)
Scores (before): A5 B4 C3 D2  ->  (after): A5 B4 C5 D5
Saw: physics is textbook (X-point reconnects, island grows, W(t) S-curve bends into
saturation) and pacing was actually fine (island visible by t~73, growth spread across
300 tau_A — my "flat then sudden" worry was wrong). BUT default white matplotlib clashing
with the new dark gallery, and the flux panel had 30 black contour lines = a busy barcode;
the two panels weren't tied together.
Defect: D-axis (default look, barcode contours, white bg) + C-axis (no saturation-onset
marker, no W<->island linkage).
Did: pure-style TWEAK. Dark house style (extended anim.apply_house_style to also style 2-D
ticks/spines/labels, re-applied each frame after ax.clear). Contour guide lines 30->12,
light + alpha 0.30 (de-barcoded). Tied the panels: a gold bar at the island O-point
(y=Ly/2) whose height IS the current W(t) — grows with the island, links to the left curve.
Marked the dW/dt-peak (saturation onset) with a line+label; filled+brightened the W(t)
curve, gold current-point marker. dpi 90->120.
Verify: regenerated (W->2.23, dW/dt peak at t=133, ratio 0.28 — unchanged physics);
re-dumped+Read — clean island, the W-bar tracks growth. pytest 208 passed; ruff clean.
No new test (style-only; apply_house_style already covered). size 2.8 MB.
Next: R3 = stellarator_flux_surfaces.

## R3 — stellarator_flux_surfaces — REWORK — (pending commit)
Scores (before): A5 B3 C3 D2  ->  (after): A5 B4 C4 D5
Saw: thin pastel loops floating tiny in a big WHITE void — pale, small, hard to read the
twist; clashed with the dark gallery. Physics correct (field lines traced on nested vacuum
surfaces, iota from geometry) but the beauty shot undersold it badly.
Defect: D-axis (the gallery's showpiece looked the weakest) + wasted frame (zlim ±2.2 vs
data ±0.6 -> donut floated small).
Did: rendering-only REWORK. Dark house style + dark 3-D panes; warm-core->cool-edge vivid
palette (gold/orange/pink/cyan); halo+core double-stroke GLOW per line (lw 3.4 alpha .10
under lw 1.15 alpha .95); 3->4 nested surfaces, 5->8 lines each (32 lines, denser weave);
tightened zlim ±2.2->±1.1 + box_aspect z .42 so it FILLS the frame; elev 26->24; full-turn
camera (seamless); dpi 90->120, fps 16->18.
Verify: regenerated (4 surfaces / 32 lines / 6 periods, zero net current — unchanged);
re-dumped+Read — glowing, centered, the iconic twisty stellarator. pytest 208; ruff clean.
No new test (render-only; helical_stellarator field + iota covered in test_fields). 4.6 MB.
Next: R4 = tokamak_discharge_full.

## R4 — tokamak_discharge_full — TWEAK — (pending commit)
Scores (before): A5 B5 C4 D2  ->  (after): A5 B5 C5 D5
Saw: the 2-panel flight-sim narrative is excellent and correct — ignition ramp (yellow
band) -> burning H-mode with T0 sawtoothing (red zigzag) + q(0)~1 (the 10*q line) ->
pellet at t=14 cools then re-climbs. Sawteeth read fine on the 1-D trace (no aliasing like
the 3-D had). BUT default white (clashing with the dark gallery) and the inferno
cross-section had NO colorbar (core temperature unreadable).
Defect: D-axis (white bg, white box around the circular cross-section) + C-axis (missing
colorbar).
Did: pure-style TWEAK. Dark house style on both panels (re-applied per frame); added a
T[keV] colorbar to the cross-section + a subtle rho=1 boundary circle, hid the box spines;
brightened the trace (T0 #ff5a5a, q(0)x10 cyan), labeled ignition band / pellet line / q=1;
styled the legend dark. dpi 90->120, fps 14->16.
Verify: regenerated (179 sawteeth, q(0) min 0.94 — unchanged); re-dumped+Read — clean dark
2-panel with colorbar. pytest 208; ruff clean. No new test (style-only). size 3.7 MB.
ALL 4 HEROES DONE (R1 3-D discharge, R2 tearing, R3 stellarator, R4 discharge). Next: R5 =
stellarator_burn (the contrast 2-panel).

## R5 — stellarator_burn — TWEAK — (pending commit)
Scores (before): A5 B5 C4 D2  ->  (after): A5 B5 C5 D5
Saw: the steady-state contrast is crisp and correct (smooth stellarator T0 vs sawtoothing
tokamak, pellet drop+reclimb on both). BUT default white (clashing), no colorbar, and the
cross-section used VIRIDIS for temperature — inconsistent with the house convention
(inferno=T, as in R4's tokamak cross-section: same quantity must use the same cmap).
Defect: D-axis (white) + C-axis (no colorbar) + cmap-convention break (viridis T).
Did: pure-style TWEAK. Dark house style both panels; cross-section viridis->INFERNO (T
convention) + T[keV] colorbar + rho=1 boundary, box spines hidden; trace brightened
(stellarator cyan #22d3ee bold, tokamak grey), legend now QUANTIFIES the contrast
(title "crashes": tokamak 179 vs stellarator 0); labeled heating band + pellet line; fixed
colorbar/ylabel overlap with wspace=0.32. dpi 90->110, fps 14->16.
Verify: regenerated (0 vs 179 sawteeth — unchanged); re-dumped+Read — clean, consistent,
the 179-vs-0 contrast explicit. pytest 208; ruff clean. No new test (style-only). 2.2 MB.
Next: R6 = kink_eigenmode (check the schematic-amplitude honesty call).

## R6 — kink_eigenmode — TWEAK — (pending commit)
Scores (before): A4 B4 C4 D2  ->  (after): A5 B4 C5 D5
Saw: left panel is EXACT (m=1 top-hat eigenfunction xi_r + q(r) crossing 1 at r1=0.548);
right shows the real crescent displacement. The amplitude readout said "amp=" (not "t=") so
it wasn't faking a clock, but it never disclosed that the growth ENVELOPE is schematic
(0.33*linspace^1.5, no gamma) — the charter requires a schematic to say so. Also default
white.
Honesty call: driving a real ideal-kink gamma isn't available (only resistive FKR is, in
R2), and ideal-kink saturation is nonlinear — so the correct move is to LABEL it schematic
(eigenfunction exact, envelope illustrative), NOT invent a growth rate. Did that.
Did: TWEAK. Added explicit "illustrative growth envelope · eigenfunction exact" footnote +
docstring note; readout now "displacement = X.XX a". Dark house style both panels (axL
styled once; axR re-styled per frame). Envelope x^1.5 -> smoothstep, phase 1.5π -> full 2π
(seamless rotation). Colors -> house palette (xi_r red, q cyan). dpi 90->120.
Verify: regenerated (q(0)=0.85<1 -> UNSTABLE, r1=0.548 — unchanged); re-dumped+Read — dark,
honest, the crescent clear. pytest 208; ruff clean. No new test (style+labeling; kink
physics covered in test_cylinder_mhd). size 1.1 MB.
Next: R7 = operating_modes.

## R7 — operating_modes — TWEAK — (pending commit)
Scores (before): A5 B5 C4 D2  ->  (after): A5 B5 C5 D5
Saw: the operating-window physics is clear + correct (L-mode stays cool, H-mode ignites
into the burning band ~22 keV, over-fuel crosses the Greenwald limit n_G and the
confinement collapses to ~0.6 keV / n_e/n_G=1.20). Well-paced (45 s, 5.5 s playback). Only
flaw: default white (clashing) and the disruption track was dark grey "0.35" (near-invisible
on dark).
Defect: D-axis (white) + a track color that wouldn't survive the dark theme.
Did: TWEAK in the SHARED renderer — added dark=True to anim.animate_operating_space
(house-style bg, gold burning band, #ff5a5a Greenwald vline, light head edges, dark legend).
Brightened the gallery track colors (L cyan, H red, disruption #c0c5cf). dpi 90->110.
Added the previously-MISSING smoke test for animate_operating_space (dark path, band+vline).
Verify: regenerated (n_G 2.23e20, P_LH~16 MW, L 5.1 / H 22.3 / disrupt 0.6 keV — unchanged);
re-dumped+Read — dark, the disruption clearly diving past n_G. pytest 209 (+1); ruff clean.
size 0.49 MB.
Next: R8 = burn_dshaped_cross_section.

## R0 — setup — done — 374bed3 (charter) / baseline green
Did: wrote scripts/_dump_frames.py; regenerated all 11 gifs fresh (gif_gallery.py all).
Baseline validation lines (each gif's own print):
- _smoke_diffusion: mass drift 9.08e-3 conserved; peak-law err 2.2e-16
- burn_0d_ignition: steady T=14.3 keV, beta 4.17% (limit 3.96%), f_He 5.2%, ash bal 0.991
- burn_1d_two_temperature: Ti0=24.2, Te0=12.8 keV, Ti/Te=1.89, tau_eq(core)=274 ms
- burn_dshaped: Shafranov +0.28 m, kappa 1.48; core 28.0 keV, P_alpha 75 MW; tau_E sim 0.63 s vs IPB98 0.19 s (honest: chi set for showcase)
- operating_modes: n_G=2.23e20; P_LH~16 MW; L 5.1 keV / H 22.3 keV / disruption 0.6 keV (n/n_G 1.20)
- kink_eigenmode: q(0)=0.85<1 -> kink UNSTABLE, q=1 at r1=0.548
- tearing_island_saturation: W->2.23 sheet widths; dW/dt peaks 1.24e-2 at t=133 tau_A, falls to 3.48e-3 (ratio 0.28, saturating)
- tokamak_discharge_full: 179 sawteeth; q(0) min 0.94; core T0 23-26 keV
- tokamak_3d_discharge: 179 sawteeth; core(rho=0.12) T 2-29 keV
- stellarator_flux_surfaces: l=2 vacuum, 3 surfaces/15 lines/6 periods, zero net current
- stellarator_burn: stellarator 0 sawteeth vs tokamak 179 (same burn)
Verify: pytest 206 passed (175 s); ruff clean.
Note: confirms the R1 aliasing hypothesis numerically — tokamak_3d_discharge fires 179
crashes onto 100 saved frames (stride ~1.8 crashes/frame). Next: R1.

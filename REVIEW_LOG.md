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

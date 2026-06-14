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

"""gif_gallery.py — regenerate every animated showcase .gif.

One function per animation, registered in GALLERY. Each prints its validation
number(s) and writes a .gif (and, for operation-mode showcases, a PNG still) to
`outputs/`. The night's autonomous build (see NIGHT.md) adds one entry per rung.
(Static-figure montage lives in the separate `gallery.py`.)

Usage:
    python gif_gallery.py                 # list available gifs
    python gif_gallery.py <name> [...]    # regenerate the named gif(s)
    python gif_gallery.py all             # regenerate everything
"""

from __future__ import annotations

import sys

import numpy as np

from plasmaplay import animate as anim, transport as tr

OUT = "outputs"


def smoke_diffusion():
    """G1 foundation smoke: a 1-D diffusing Gaussian. Validates the gif pipeline
    against the analytic diffusion solution (mass conservation + peak-decay law)."""
    x = np.linspace(-8, 8, 161)
    times = np.linspace(0, 6, 90)
    D, s0 = 0.7, 1.0

    def gauss(t):
        s2 = s0**2 + 2 * D * t
        return 1.0 / np.sqrt(2 * np.pi * s2) * np.exp(-x**2 / (2 * s2))

    frames = anim.make_frames(gauss, times)
    mass = np.trapezoid(frames, x, axis=1)
    peak = frames.max(axis=1)
    exp_peak = peak[0] / np.sqrt(1 + 2 * D * times / s0**2)
    print(f"  [smoke_diffusion] mass drift = {np.abs(mass / mass[0] - 1).max():.2e} "
          f"(conserved); peak-law err = {np.abs(peak / exp_peak - 1).max():.2e}")
    out = anim.animate_profiles(x, frames, times, path=f"{OUT}/_smoke_diffusion.gif",
                                ylabel="c(x,t)", title="G1 smoke: 1-D diffusion",
                                fps=20, dpi=90)
    print(f"  wrote {out}")


def burn_0d_ignition():
    """A1 (F1): 0-D burn with He ash + beta-limit. Phase-space (n,T) ignition track
    igniting onto the beta-limited operating point, the point colored by ash
    fraction as ash builds up. Validates the steady ash balance n_He = tau_he*R."""
    beta_lim = tr.troyon_limit(3.0, 7.0, 1.0, 5.3)   # ~3.96%

    def paux(t):
        return 5.0e5 if t < 5.0 else 0.0

    r = tr.burn_0d_ash(1.0e20, 5.0, tau_E=3.0, p_aux=paux, B=5.3, tau_p=6.0,
                       tau_he=10.0, fuel_rate=1.0e20 / 6.0, beta_limit=beta_lim,
                       t_end=40.0)
    R = tr.reaction_rate_dt(r["n_DT"][-1], r["T"][-1])
    print(f"  [burn_0d_ignition] steady T = {r['T'][-1]:.1f} keV, "
          f"beta = {r['beta'][-1]*100:.2f}% (limit {beta_lim*100:.2f}%), "
          f"f_He = {r['f_He'][-1]*100:.1f}%, ash balance = "
          f"{r['n_He'][-1]/(10.0*R):.3f}")
    # subsample to a sane frame count
    s = slice(0, None, max(1, r["t"].size // 100))
    out = anim.animate_phase_track(
        r["n_e"][s], r["T"][s], r["t"][s], path=f"{OUT}/burn_0d_ignition.gif",
        color=r["f_He"][s] * 100, xlabel=r"$n_e$ [m$^{-3}$]", ylabel="T [keV]",
        clabel="ash fraction [%]", title="0-D ignition onto the burning point",
        ylim=(0, 20), fps=20, dpi=90)
    print(f"  wrote {out}")


def burn_1d_two_temperature():
    """A2 (F2.5): 1-D two-temperature burn. Neutral beams heat the ions, fusion
    alphas heat the electrons, and collisional (Spitzer) equipartition couples the
    two channels — so the beam-heated plasma settles at Ti > Te. The gif overlays
    Te(rho,t) and Ti(rho,t). Validates the equipartition time against the formula."""
    sim = tr.TwoTempTransport1D(a=1.0, n_grid=129, chi_e=0.8, chi_i=0.4, D=0.06,
                                mu_i=2.5, Te_edge=0.1, Ti_edge=0.1, n_edge=2e19)
    sim.set_state(Te=2.0, Ti=2.0, n=8e19)
    n_target, tau_p = 8.0e19, 6.0
    hold = tr.gaussian_deposition(sim.rho, 0.0, 0.4)
    nbi = tr.gaussian_deposition(sim.rho, 0.0, 0.35)
    dt, t_end = 4e-3, 12.0
    nsteps = int(round(t_end / dt))
    stride = max(1, nsteps // 100)
    times, te_fr, ti_fr = [], [], []
    for k in range(nsteps):
        t = sim.t
        p_i = 6.0e5 * (0.3 + 0.7 * min(t / 4.0, 1.0))    # NBI ion heating ramps, stays on
        p_e = 1.0e5                                       # modest RF/ohmic electron heating
        sim.step(dt, p_aux_i_total=p_i, p_aux_e_total=p_e, aux_i_profile=nbi,
                 frac_alpha_e=0.85, fuel_total=n_target / tau_p, fuel_profile=hold)
        if k % stride == 0:
            times.append(sim.t); te_fr.append(sim.Te.copy()); ti_fr.append(sim.Ti.copy())
    te_fr, ti_fr, times = np.array(te_fr), np.array(ti_fr), np.array(times)
    d = sim.diagnostics()
    tau_eq = tr.equipartition_time(d["n0"], d["Te0"])
    print(f"  [burn_1d_two_temperature] steady Ti0 = {d['Ti0']:.1f} keV, "
          f"Te0 = {d['Te0']:.1f} keV, Ti/Te = {d['Ti0']/d['Te0']:.2f}, "
          f"tau_eq(core) = {tau_eq*1e3:.0f} ms")
    frames = np.stack([te_fr, ti_fr], axis=1)            # (n_t, 2, n_rho)
    out = anim.animate_profiles(
        sim.rho, frames, times, path=f"{OUT}/burn_1d_two_temperature.gif",
        labels=[r"$T_e$ (electrons)", r"$T_i$ (ions, NBI-heated)"],
        xlabel=r"$\rho = r/a$", ylabel="T [keV]",
        title="Two-temperature burn: ions hotter than electrons", fps=20, dpi=90)
    print(f"  wrote {out}")


GALLERY = {
    "smoke_diffusion": smoke_diffusion,
    "burn_0d_ignition": burn_0d_ignition,
    "burn_1d_two_temperature": burn_1d_two_temperature,
}


def main(argv):
    if not argv:
        print("Available gifs:")
        for name in GALLERY:
            print(f"  {name}")
        print("\nUsage: python gif_gallery.py <name>[ ...] | all")
        return
    names = list(GALLERY) if argv == ["all"] else argv
    for name in names:
        if name not in GALLERY:
            print(f"  ! unknown gif '{name}' (skipped)")
            continue
        print(f"[{name}]")
        GALLERY[name]()


if __name__ == "__main__":
    main(sys.argv[1:])

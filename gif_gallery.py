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

import matplotlib

matplotlib.use("Agg", force=False)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter  # noqa: E402

from plasmaplay import (  # noqa: E402
    animate as anim, cylinder_mhd as cm, equilibrium_metrics as em,
    operating_limits as ol, reduced_mhd as rm, sawtooth as sw, transport as tr,
)
from plasmaplay.solvers import grad_shafranov_solve  # noqa: E402

OUT = "outputs"

# A small ITER-like toy device, shared by the operating-mode scenarios.
_R0, _A, _B, _IP, _KAPPA = 3.0, 1.0, 5.3, 7.0, 1.5
_S = 4 * np.pi ** 2 * _R0 * _A * np.sqrt((1 + _KAPPA ** 2) / 2)
_VOL = 2 * np.pi ** 2 * _R0 * _A ** 2
_N_G = ol.greenwald_density(_IP, _A)


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


def _dshaped_equilibrium(n=141):
    """A shaped Solov'ev Grad-Shafranov equilibrium and its flux-surface metrics.

    Reuses the validated `grad_shafranov_solve` (exp 04). Returns R, Z, psi, the
    metrics dict, and the chosen boundary flux. The plasma boundary is taken at a
    fraction of psi_axis so the last closed surface sits inside the rectangular box.
    """
    R0 = 3.0
    R = np.linspace(R0 - 1.3, R0 + 1.3, n)
    Z = np.linspace(-1.9, 1.9, n)
    RR, _ = np.meshgrid(R, Z, indexing="ij")
    psi = grad_shafranov_solve(R, Z, -(RR ** 2 + 1.0), boundary=0.0)
    psi_b = 0.06 * psi.max()
    m = em.flux_surface_metrics(R, Z, psi, n_rho=80, psi_bnd=psi_b)
    return R, Z, psi, m, psi_b


def burn_dshaped_cross_section():
    """A3 (F3): burning-plasma transport on a REAL Grad-Shafranov equilibrium.

    Runs 1-D flux-surface-averaged transport on the V'(rho), <|grad rho|^2> metrics
    of a shaped Solov'ev equilibrium, then maps T(rho, t) back onto the actual (R,Z)
    flux surfaces — the headline 'watch the D-shaped plasma burn' movie. Validates
    the Shafranov shift and the IPB98(y,2) confinement time."""
    R, Z, psi, m, _ = _dshaped_equilibrium()
    rho_grid, psi_n = m["rho_grid"], m["psi_n"]
    inside = (psi_n >= 0.0) & (psi_n <= 1.0)
    i_ax, _ = np.unravel_index(np.argmax(psi), psi.shape)
    kappa = np.ptp(Z[(psi_n <= 1).any(axis=0)]) / np.ptp(R[(psi_n <= 1).any(axis=1)])

    fs = tr.FluxSurfaceTransport1D(1.0, m["Vprime"], m["grad_rho2"], rho_metric=m["rho"],
                                   n_grid=129, chi=0.6, D=0.05, T_edge=0.1, n_edge=2e19)
    fs.set_state(T=2.0, n=5e19)
    dt, t_end = 4e-3, 16.0
    nsteps = int(round(t_end / dt))
    stride = max(1, nsteps // 100)
    times, field_fr = [], []
    for k in range(nsteps):
        paux = 5.0e5 * min(fs.t / 4.0, 1.0)               # NBI ramp, then held on
        fs.step(dt, p_aux_total=paux, fuel_total=5e19 / 6.0,
                fuel_profile=tr.gaussian_deposition(fs.rho, 0.0, 0.4))
        if k % stride == 0:
            T2d = np.interp(rho_grid.ravel(), fs.rho, fs.T).reshape(rho_grid.shape)
            times.append(fs.t); field_fr.append(T2d)
    field_fr, times = np.array(field_fr), np.array(times)

    p_alpha = fs._vol_avg(tr.fusion_power_density(fs.n, fs.T, "alpha")) * fs.plasma_volume()
    p_in = p_alpha + 5.0e5 * fs.plasma_volume()
    tau_sim = fs.energy_confinement_time(p_in)
    tau_98 = em.confinement_time_ipb98(Ip_MA=7.0, B=5.3, n19=5.0, P_MW=p_in / 1e6,
                                       R=3.0, a=1.0, kappa=kappa, M=2.5)
    print(f"  [burn_dshaped] Shafranov shift = {R[i_ax]-3.0:+.2f} m, elongation kappa = {kappa:.2f}")
    print(f"  [burn_dshaped] steady core T0 = {fs.T[0]:.1f} keV, P_alpha = {p_alpha/1e6:.0f} MW")
    print(f"  [burn_dshaped] tau_E: sim = {tau_sim:.2f} s, IPB98(y,2) for these params "
          f"= {tau_98:.2f} s (within a factor of a few; chi is set for the showcase, "
          f"not fit to the scaling — IPB98 is validated against ITER = 3.7 s in the tests)")
    out = anim.animate_poloidal_field(
        R, Z, field_fr, times, path=f"{OUT}/burn_dshaped_cross_section.gif",
        mask=inside, clabel="T [keV]", cmap="inferno", vmax=float(field_fr.max()),
        title="Burn on the real D-shaped equilibrium", fps=20, dpi=90)
    print(f"  wrote {out}")


def _modes_tau_factor(n0):
    """Confinement multiplier for this device: L->H bifurcation x Greenwald collapse."""
    p_lh = ol.lh_power_threshold(n0 / 1e20, _B, _S)

    def factor(t, n_e, T, p_heat_density):
        return (ol.confinement_factor_lh(p_heat_density * _VOL / 1e6, p_lh)
                * ol.confinement_factor_greenwald(n_e, _N_G))
    return factor


def _modes_scenario(n0, p_aux, tau_E, fuel_rate, t_end=45.0):
    return tr.burn_0d_ash(n0, 3.0, tau_E=tau_E, p_aux=p_aux, B=_B, tau_p=6.0,
                          tau_he=10.0, fuel_rate=fuel_rate, beta_limit=0.04,
                          tau_factor=_modes_tau_factor(n0), t_end=t_end, dt=1e-3)


def operating_modes():
    """A4 (F3.5): the tokamak operating WINDOW, not one happy path. Three 0-D burns
    sweep an (n_e, T) operating diagram: an L-mode (heating below the L->H threshold,
    stays cool), an H-mode (above threshold, ignites into the beta-limited burning
    band), and a disruption (over-fuelled past the Greenwald limit -> the confinement
    collapses and the burn dies). Validates n_G and the L->H power threshold."""
    p_lh = ol.lh_power_threshold(0.7, _B, _S)
    print(f"  [operating_modes] n_G = {_N_G:.2e} m^-3, L->H threshold P_LH ~ {p_lh:.0f} MW")
    lmode = _modes_scenario(5e19, lambda t: 1.2e5, 1.0, 5e19 / 6)
    hmode = _modes_scenario(7e19, lambda t: 3e5, 1.8, 7e19 / 6)
    disrupt = _modes_scenario(7e19, lambda t: 3e5, 1.8,
                              lambda t: 7e19 / 6 if t < 22 else 4.5e19)
    for name, r in (("L-mode", lmode), ("H-mode", hmode), ("disruption", disrupt)):
        print(f"    {name:11s}: end T0 = {r['T'][-1]:5.1f} keV, "
              f"n_e/n_G = {r['n_e'][-1]/_N_G:.2f}, beta = {r['beta'][-1]*100:.1f}%")

    sl = slice(0, None, max(1, lmode["t"].size // 110))
    tracks = [
        {"x": lmode["n_e"][sl], "y": lmode["T"][sl], "label": "L-mode (sub-threshold)",
         "color": "tab:blue"},
        {"x": hmode["n_e"][sl], "y": hmode["T"][sl], "label": "H-mode (burning)",
         "color": "tab:red"},
        {"x": disrupt["n_e"][sl], "y": disrupt["T"][sl],
         "label": "over-fuel -> disruption", "color": "0.35"},
    ]
    out = anim.animate_operating_space(
        tracks, lmode["t"][sl], path=f"{OUT}/operating_modes.gif",
        xlabel=r"$n_e$ [m$^{-3}$]", ylabel="T [keV]", title="Tokamak operating modes",
        vlines=[(_N_G, r"Greenwald limit $n_G$")], band=(10.0, 25.0),
        xlim=(0, 1.25 * _N_G), ylim=(0, 30), fps=20, dpi=90)
    print(f"  wrote {out}")


def kink_eigenmode():
    """B1 (cylindrical MHD): the m=1/n=1 internal kink — the sawtooth trigger. With
    q(0)<1 a q=1 surface exists and the core inside it shifts rigidly sideways. Left:
    the radial eigenfunction xi_r(r) and q(r). Right: the poloidal cross-section, the
    hot core displaced into the characteristic crescent, the displacement growing and
    rotating. Validates that the kink is unstable exactly when q(0)<1."""
    q0, nu = 0.85, 1.0
    r1 = cm.rational_surface(1, 1, q0, nu)
    print(f"  [kink_eigenmode] q(0)={q0} < 1 -> internal kink UNSTABLE, "
          f"q=1 surface at r1={r1:.3f} (kink_unstable={cm.internal_kink_unstable(q0)})")

    nr, nth, nframes = 60, 120, 90
    rg = np.linspace(1e-3, 1.0, nr)
    th = np.linspace(0, 2 * np.pi, nth)
    RR, TH = np.meshgrid(rg, th, indexing="ij")
    xi_r = cm.internal_kink_xi(rg, r1)              # ideal m=1 top-hat displacement
    XI = cm.internal_kink_xi(RR, r1)
    Fcore = 1.0 - RR ** 2                            # a "temperature" peaked on axis
    q_of_r = cm.screw_pinch_q(rg, q0, nu)
    amp = 0.33 * (np.linspace(0, 1, nframes) ** 1.5)   # growing mode amplitude
    phase = np.linspace(0, 1.5 * np.pi, nframes)       # helical rotation

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 5.2))
    axL.plot(rg, xi_r, color="crimson", lw=2, label=r"$\xi_r(r)$ (kink)")
    axL.plot(rg, q_of_r, color="navy", lw=1.5, label="q(r)")
    axL.axhline(1.0, color="0.6", ls=":", lw=0.9)
    axL.axvline(r1, color="k", ls="--", lw=0.9)
    axL.text(r1, 1.7, "q=1", rotation=90, va="bottom", ha="right", fontsize=8)
    axL.set(xlabel="r/a", ylim=(0, 2.0), title="m=1 internal kink: eigenfunction & q(r)")
    axL.legend(loc="center right", fontsize=8)
    axR.set_aspect("equal")

    def draw(i):
        axR.clear(); axR.set_aspect("equal"); axR.set_xticks([]); axR.set_yticks([])
        A, ph = amp[i], phase[i]
        X = RR * np.cos(TH) + A * XI * np.cos(ph)
        Y = RR * np.sin(TH) + A * XI * np.sin(ph)
        axR.contourf(X, Y, Fcore, levels=40, cmap="inferno")
        # the q=1 surface (undisplaced reference circle)
        axR.plot(r1 * np.cos(th), r1 * np.sin(th), color="cyan", lw=1.0, ls="--")
        axR.set_xlim(-1.15, 1.15); axR.set_ylim(-1.15, 1.15)
        axR.set_title(f"core displaced (m=1), amp={A:.2f}")

    an = FuncAnimation(fig, draw, frames=nframes, blit=False)
    out = f"{OUT}/kink_eigenmode.gif"
    an.save(out, writer=PillowWriter(fps=18), dpi=90)
    plt.close(fig)
    print(f"  wrote {out}")


def tearing_island_saturation():
    """B2 (reduced MHD): a tearing mode reconnects the Harris sheet into a magnetic
    island that GROWS then SATURATES (the Rutherford regime). Left: the island width
    W(t), rising then bending over as dW/dt turns down. Right: the flux contours, the
    neutral line tearing open into the island and settling. Validated by the dW/dt
    turnover test; the linear phase obeys the FKR S^-3/5 law."""
    k = 0.5
    sim = rm.ReducedMHD(k, S=100.0, Pm=0.0, nx=160, ny=64, Lx=4.0).seed(1e-3)
    dt, t_end, nframes = 0.012, 300.0, 90
    stride = max(1, int(t_end / dt) // nframes)
    yext = np.concatenate([sim.y, [sim.Ly]])
    frames, times, Wt = [], [], []
    for i in range(int(t_end / dt)):
        sim.step(dt)
        if i % stride == 0:
            psi = sim.flux_function()
            frames.append(np.concatenate([psi, psi[:, :1]], axis=1))
            times.append(sim.t); Wt.append(sim.island_width())
    times, Wt = np.array(times), np.array(Wt)
    dWdt = np.gradient(Wt, times)
    i_peak = int(np.argmax(dWdt))
    print(f"  [tearing_island_saturation] k={k}: W grows to {Wt[-1]:.2f} sheet widths; "
          f"dW/dt peaks {dWdt[i_peak]:.2e} at t={times[i_peak]:.0f} then falls to "
          f"{dWdt[-1]:.2e} (saturating, ratio {dWdt[-1]/dWdt[i_peak]:.2f})")

    fig, (axW, axF) = plt.subplots(1, 2, figsize=(11, 5.0),
                                   gridspec_kw={"width_ratios": [1, 1.1]})
    mask = np.abs(sim.x) <= 2.5

    def draw(j):
        axW.clear(); axF.clear()
        axW.plot(times[:j + 1], Wt[:j + 1], color="crimson", lw=2)
        axW.scatter([times[j]], [Wt[j]], color="crimson", zorder=3)
        axW.set(xlim=(0, times[-1]), ylim=(0, Wt.max() * 1.1),
                xlabel=r"$t / \tau_A$", ylabel="island width W",
                title="Rutherford saturation")
        axF.contourf(yext, sim.x[mask], frames[j][mask], levels=40, cmap="RdBu_r")
        axF.contour(yext, sim.x[mask], frames[j][mask], levels=30, colors="k",
                    linewidths=0.5)
        axF.axhline(0.0, color="0.4", ls=":", lw=0.8)
        axF.set(xlabel="y", ylabel="x", title=f"flux contours   t = {times[j]:.0f}")

    an = FuncAnimation(fig, draw, frames=len(frames), blit=False)
    out = f"{OUT}/tearing_island_saturation.gif"
    an.save(out, writer=PillowWriter(fps=16), dpi=90)
    plt.close(fig)
    print(f"  wrote {out}")


def tokamak_discharge_full():
    """Track C (the integrated discharge): the F2 transport burn (seconds) coupled to
    m=1 SAWTOOTH crashes (instantaneous on that scale) — ignition -> burning H-mode
    with periodic sawteeth -> pellet fuel injection -> settling. Left: the poloidal
    cross-section T(rho,t). Right: the core T0(t) sawtoothing and q(0) crossing 1.
    The headline two-timescale 'flight simulator' (staged coupling, named as such)."""
    a = 1.0
    sim = tr.Transport1D(a, n_grid=129, chi=0.10, D=0.04, T_edge=0.1, n_edge=2e19,
                         B=5.3, beta_limit=0.04, beta_stiffness=40.0)
    sim.set_state(T=2.0, n=1e20)
    aux = tr.gaussian_deposition(sim.rho, 0.0, 0.35)
    hold = tr.gaussian_deposition(sim.rho, 0.0, 0.40)
    pellet = tr.gaussian_deposition(sim.rho, 0.35, 0.12)
    dt, t_end = 2e-3, 22.0
    stride = max(1, int(t_end / dt) // 110)
    ts, T0, q0, Tfr, crmark, n_saw = [], [], [], [], [], 0
    for k in range(int(t_end / dt)):
        t = sim.t
        p_aux = 6e5 * (0.3 + 0.7 * min(t / 4.0, 1.0)) if t < 4.0 else 0.0
        ft, fp = (1e20 / 6.0, hold)
        if 14.0 <= t < 14.2:
            ft, fp = 1e20 / 6.0 + 3e20, pellet
        sim.step(dt, p_aux_total=p_aux, aux_profile=aux, fuel_total=ft, fuel_profile=fp)
        n2, T2, crashed = sw.sawtooth_event(sim.rho, sim.n, sim.T, q_edge=2.2)
        if crashed:
            sim.n, sim.T = n2, T2; n_saw += 1
        if k % stride == 0:
            ts.append(t); T0.append(sim.T[0]); Tfr.append(sim.T.copy())
            q0.append(sw.q_from_temperature(sim.rho, sim.T, 2.2)[0])
            crmark.append(crashed)
    ts, T0, q0, Tfr = np.array(ts), np.array(T0), np.array(q0), np.array(Tfr)
    print(f"  [tokamak_discharge_full] {n_saw} sawteeth; q(0) min {q0.min():.2f}; "
          f"core T0 sawtooths {T0[(ts > 5) & (ts < 14)].min():.0f}-"
          f"{T0[(ts > 5) & (ts < 14)].max():.0f} keV")

    theta = np.linspace(0, 2 * np.pi, 160)
    RR, TT = np.meshgrid(sim.rho, theta)
    X, Y = RR * np.cos(TT), RR * np.sin(TT)
    vmax = float(Tfr.max())
    fig, (axc, axt) = plt.subplots(1, 2, figsize=(11, 5.0),
                                   gridspec_kw={"width_ratios": [1, 1.25]})
    axc.set_aspect("equal")

    def draw(i):
        axc.clear(); axc.set_aspect("equal"); axc.set_xticks([]); axc.set_yticks([])
        axc.contourf(X, Y, np.broadcast_to(Tfr[i], RR.shape), levels=40,
                     cmap="inferno", vmin=0, vmax=vmax)
        axc.set_title(f"cross-section  t = {ts[i]:.1f} s")
        axt.clear()
        axt.axvspan(0, 4, color="gold", alpha=0.15)
        axt.axvline(14.0, color="purple", ls="--", lw=1.0)
        axt.plot(ts[:i + 1], T0[:i + 1], color="crimson", lw=1.0, label="core T0 [keV]")
        axt.plot(ts[:i + 1], 10 * q0[:i + 1], color="navy", lw=0.9, label="10*q(0)")
        axt.axhline(10.0, color="navy", ls=":", lw=0.8)
        axt.set(xlim=(0, t_end), ylim=(0, max(35, vmax * 1.1)), xlabel="t [s]",
                title="ignition | burning H-mode + sawteeth | pellet")
        axt.legend(loc="upper right", fontsize=8)

    an = FuncAnimation(fig, draw, frames=len(ts), blit=False)
    out = f"{OUT}/tokamak_discharge_full.gif"
    an.save(out, writer=PillowWriter(fps=14), dpi=90)
    plt.close(fig)
    print(f"  wrote {out}")


def tokamak_3d_discharge():
    """C2 (the 3-D showcase): the event-coupled discharge rendered on the 3-D torus —
    nested flux surfaces colored by their temperature, rotating, with the burn heating
    the core and the sawtooth crashes visibly flattening it. Reuses the validated
    coupled-discharge data (Track C) and `animate.animate_torus_nested`."""
    sim = tr.Transport1D(1.0, n_grid=129, chi=0.10, D=0.04, T_edge=0.1, n_edge=2e19,
                         B=5.3, beta_limit=0.04, beta_stiffness=40.0)
    sim.set_state(T=2.0, n=1e20)
    aux = tr.gaussian_deposition(sim.rho, 0.0, 0.35)
    hold = tr.gaussian_deposition(sim.rho, 0.0, 0.40)
    pellet = tr.gaussian_deposition(sim.rho, 0.35, 0.12)
    rho_levels = np.array([0.12, 0.3, 0.5, 0.72, 1.0])
    dt, t_end = 2e-3, 22.0
    stride = max(1, int(t_end / dt) // 100)
    times, T_rt, n_saw = [], [], 0
    for k in range(int(t_end / dt)):
        t = sim.t
        p_aux = 6e5 * (0.3 + 0.7 * min(t / 4.0, 1.0)) if t < 4.0 else 0.0
        ft, fp = (1e20 / 6.0, hold)
        if 14.0 <= t < 14.2:
            ft, fp = 1e20 / 6.0 + 3e20, pellet
        sim.step(dt, p_aux_total=p_aux, aux_profile=aux, fuel_total=ft, fuel_profile=fp)
        n2, T2, crashed = sw.sawtooth_event(sim.rho, sim.n, sim.T, q_edge=2.2)
        if crashed:
            sim.n, sim.T = n2, T2; n_saw += 1
        if k % stride == 0:
            times.append(t)
            T_rt.append(np.interp(rho_levels, sim.rho, sim.T))
    times, T_rt = np.array(times), np.array(T_rt)
    print(f"  [tokamak_3d_discharge] {n_saw} sawteeth; core(rho=0.12) T "
          f"{T_rt[:, 0].min():.0f}-{T_rt[:, 0].max():.0f} keV over the discharge")
    out = anim.animate_torus_nested(
        rho_levels, T_rt, times, path=f"{OUT}/tokamak_3d_discharge.gif",
        R0=3.0, a=1.0, title="3-D discharge: ignition -> burn + sawteeth -> pellet",
        vmax=float(T_rt.max()), fps=14, dpi=90)
    print(f"  wrote {out}")


GALLERY = {
    "smoke_diffusion": smoke_diffusion,
    "burn_0d_ignition": burn_0d_ignition,
    "burn_1d_two_temperature": burn_1d_two_temperature,
    "burn_dshaped_cross_section": burn_dshaped_cross_section,
    "operating_modes": operating_modes,
    "kink_eigenmode": kink_eigenmode,
    "tearing_island_saturation": tearing_island_saturation,
    "tokamak_discharge_full": tokamak_discharge_full,
    "tokamak_3d_discharge": tokamak_3d_discharge,
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
